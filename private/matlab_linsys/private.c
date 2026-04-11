#include "private.h"
#include <string.h>

const char *scs_get_lin_sys_method(void) {
  return "sparse-direct-matlab-ldl";
}

/* Convert upper-triangular ScsMatrix (CSC) to full symmetric MATLAB sparse.
 * For each off-diagonal entry (i,j) with i < j, we store both (i,j) and (j,i).
 * This avoids calling MATLAB functions for the symmetrization.
 *
 * NOTE: MATLAB's documentation for sparse ldl says "only the upper triangle
 * is referenced", which suggests passing just the upper triangle should work.
 * However, testing shows that MATLAB's sparse ldl (which uses MA57 from HSL)
 * requires the input to be *structurally* symmetric — i.e., if entry (i,j)
 * exists then (j,i) must also exist. Passing only the upper triangle of an
 * asymmetric sparsity pattern (like a KKT matrix where A' occupies the upper-
 * right and A occupies the lower-left) causes ldl to silently produce wrong
 * factors with no error or warning. See test/ldl_diag.m for a demonstration.
 *
 * The "upper triangle only" documentation likely means that for a structurally
 * symmetric input, ldl only *reads* the upper triangle for efficiency. It does
 * NOT mean you can omit the lower triangle entries from the sparsity pattern. */
static mxArray *scs_to_mxsparse_symmetric(const ScsMatrix *M) {
  scs_int n = M->n;
  scs_int j, k, nnz_upper, nnz_full;
  scs_int *col_counts;
  mxArray *mx;
  double *pr;
  mwIndex *ir, *jc;
  scs_int *write_pos;

  nnz_upper = M->p[n];

  /* Count off-diagonal entries to compute full nnz */
  nnz_full = nnz_upper; /* start with upper triangle entries */
  for (j = 0; j < n; j++) {
    for (k = M->p[j]; k < M->p[j + 1]; k++) {
      if (M->i[k] != j) {
        nnz_full++; /* each off-diagonal entry appears twice */
      }
    }
  }

  /* Count entries per column in the full symmetric matrix */
  col_counts = (scs_int *)scs_calloc(n, sizeof(scs_int));
  if (!col_counts) return SCS_NULL;

  for (j = 0; j < n; j++) {
    for (k = M->p[j]; k < M->p[j + 1]; k++) {
      scs_int i = M->i[k];
      col_counts[j]++; /* upper triangle entry (i,j) goes in column j */
      if (i != j) {
        col_counts[i]++; /* mirror entry (j,i) goes in column i */
      }
    }
  }

  /* Create MATLAB sparse matrix and build column pointers */
  mx = mxCreateSparse((mwSize)M->m, (mwSize)n, (mwSize)nnz_full, mxREAL);
  if (!mx) {
    scs_free(col_counts);
    return SCS_NULL;
  }

  pr = mxGetPr(mx);
  ir = mxGetIr(mx);
  jc = mxGetJc(mx);

  jc[0] = 0;
  for (j = 0; j < n; j++) {
    jc[j + 1] = jc[j] + (mwIndex)col_counts[j];
  }

  /* Fill entries: write_pos[j] tracks next write position for column j */
  write_pos = (scs_int *)scs_calloc(n, sizeof(scs_int));
  if (!write_pos) {
    scs_free(col_counts);
    mxDestroyArray(mx);
    return SCS_NULL;
  }

  for (j = 0; j < n; j++) {
    write_pos[j] = (scs_int)jc[j];
  }

  for (j = 0; j < n; j++) {
    for (k = M->p[j]; k < M->p[j + 1]; k++) {
      scs_int i = M->i[k];
      scs_int pos;
      /* Original entry (i,j) in column j, i <= j (upper triangle) */
      pos = write_pos[j]++;
      ir[pos] = (mwIndex)i;
      pr[pos] = (double)M->x[k];
      if (i != j) {
        /* Mirror entry (j,i) in column i */
        pos = write_pos[i]++;
        ir[pos] = (mwIndex)j;
        pr[pos] = (double)M->x[k];
      }
    }
  }

  scs_free(col_counts);
  scs_free(write_pos);

  /* Row indices are already sorted by construction: for each column C,
   * upper-triangle entries have row indices <= C (from the input CSC order)
   * and mirror entries have row indices > C (from increasing outer loop j).
   * Concatenating these two sequences produces a sorted array. */

  return mx;
}

/* Extract L factor from MATLAB sparse mxArray into C ScsMatrix.
 * MATLAB's ldl returns L with unit diagonal stored; we strip the diagonal
 * to match QDLDL convention (only strictly lower-triangular entries). */
static scs_int extract_L(ScsLinSysWork *p, const mxArray *L_mx) {
  scs_int n_plus_m = p->n + p->m;
  mwIndex *jc = mxGetJc(L_mx);
  mwIndex *ir = mxGetIr(L_mx);
  double *pr = mxGetPr(L_mx);
  scs_int j, nnz_nodiag;
  mwIndex k;

  /* Count non-diagonal entries */
  nnz_nodiag = 0;
  for (j = 0; j < n_plus_m; j++) {
    for (k = jc[j]; k < jc[j + 1]; k++) {
      if ((scs_int)ir[k] != j) {
        nnz_nodiag++;
      }
    }
  }

  /* Free old L if present */
  if (p->L) {
    SCS(cs_spfree)(p->L);
  }

  /* Allocate new L */
  p->L = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  if (!p->L) {
    return -1;
  }
  p->L->m = n_plus_m;
  p->L->n = n_plus_m;
  p->L->p = (scs_int *)scs_calloc(n_plus_m + 1, sizeof(scs_int));
  if (nnz_nodiag > 0) {
    p->L->i = (scs_int *)scs_calloc(nnz_nodiag, sizeof(scs_int));
    p->L->x = (scs_float *)scs_calloc(nnz_nodiag, sizeof(scs_float));
  }
  if (!p->L->p || (nnz_nodiag > 0 && (!p->L->i || !p->L->x))) {
    SCS(cs_spfree)(p->L);
    p->L = SCS_NULL;
    return -1;
  }

  /* Fill L, skipping diagonal entries */
  {
    scs_int write_idx = 0;
    for (j = 0; j < n_plus_m; j++) {
      p->L->p[j] = write_idx;
      for (k = jc[j]; k < jc[j + 1]; k++) {
        if ((scs_int)ir[k] != j) {
          p->L->i[write_idx] = (scs_int)ir[k];
          p->L->x[write_idx] = (scs_float)pr[k];
          write_idx++;
        }
      }
    }
    p->L->p[n_plus_m] = write_idx;
  }
  return 0;
}

/* Extract D diagonal and sub-diagonal from MATLAB sparse mxArray.
 * D from ldl() is block diagonal with 1x1 and 2x2 blocks (tridiagonal). */
static scs_int extract_D(ScsLinSysWork *p, const mxArray *D_mx) {
  scs_int n_plus_m = p->n + p->m;
  mwIndex *jc = mxGetJc(D_mx);
  mwIndex *ir = mxGetIr(D_mx);
  double *pr = mxGetPr(D_mx);
  scs_int j;
  mwIndex k;

  /* Zero out arrays */
  memset(p->D_diag, 0, n_plus_m * sizeof(scs_float));
  memset(p->D_sub, 0, (n_plus_m - 1) * sizeof(scs_float));

  /* Extract entries from sparse D */
  for (j = 0; j < n_plus_m; j++) {
    for (k = jc[j]; k < jc[j + 1]; k++) {
      scs_int row = (scs_int)ir[k];
      if (row == j) {
        p->D_diag[j] = (scs_float)pr[k];
      } else if (row == j + 1) {
        /* Sub-diagonal entry D(j+1, j) */
        p->D_sub[j] = (scs_float)pr[k];
      }
    }
  }
  return 0;
}

/* Extract permutation vector from MATLAB double array to 0-indexed C array. */
static void extract_perm(ScsLinSysWork *p, const mxArray *perm_mx) {
  scs_int n_plus_m = p->n + p->m;
  double *pr = mxGetPr(perm_mx);
  scs_int i;

  for (i = 0; i < n_plus_m; i++) {
    p->perm[i] = (scs_int)(pr[i] - 1); /* MATLAB 1-indexed to C 0-indexed */
  }
}

/* Factorize KKT matrix using MATLAB's ldl().
 * Calls [L, D, perm] = ldl(K_sym, 'vector') and extracts factors into C.
 * The KKT is stored as upper triangular in C; we build the full symmetric
 * MATLAB sparse matrix in C before calling ldl (see scs_to_mxsparse_symmetric
 * for why the full symmetric matrix is required).
 *
 * This computes a full factorization including a new fill-reducing (AMD)
 * permutation every time it is called. Ideally we would separate the symbolic
 * analysis (which depends only on sparsity pattern) from the numeric
 * factorization (which depends on values). Unfortunately MATLAB's ldl() does
 * not expose separate symbolic/numeric phases — and while MATLAB does provide
 * amd() for computing orderings independently, there is no way to pass a
 * pre-computed permutation into ldl(). Under the hood ldl() uses HSL's MA57,
 * which does have separate analyze/factorize phases, but MATLAB's wrapper
 * bundles them into a single call. */
static scs_int matlab_ldl_factor(ScsLinSysWork *p) {
  mxArray *K_sym, *rhs[2], *lhs[3];

  /* Build full symmetric MATLAB sparse from upper-triangular C CSC */
  K_sym = scs_to_mxsparse_symmetric(p->kkt);

  /* [L, D, perm] = ldl(K_sym, 'vector') */
  rhs[0] = K_sym;
  rhs[1] = mxCreateString("vector");

  {
    mxArray *err = mexCallMATLABWithTrap(3, lhs, 2, rhs, "ldl");
    if (err != NULL) {
      scs_printf("Error in MATLAB ldl() factorization.\n");
      mxDestroyArray(K_sym);
      mxDestroyArray(rhs[1]);
      mxDestroyArray(err);
      return -1;
    }
  }

  /* Extract factors into C arrays */
  {
    scs_int status = 0;
    if (extract_L(p, lhs[0]) < 0) {
      status = -1;
    } else {
      extract_D(p, lhs[1]);
      extract_perm(p, lhs[2]);
    }

    /* Clean up all MATLAB temporaries */
    mxDestroyArray(K_sym);
    mxDestroyArray(rhs[1]);
    mxDestroyArray(lhs[0]);
    mxDestroyArray(lhs[1]);
    mxDestroyArray(lhs[2]);

    if (status < 0) {
      return -1;
    }
  }

  p->factorizations++;
  return 0;
}

/* Forward solve: (L + I) * x = b, where L is strictly lower-triangular.
 * Overwrites b with the solution. Same as QDLDL_Lsolve. */
static void forward_solve(scs_int n, const scs_int *Lp, const scs_int *Li,
                          const scs_float *Lx, scs_float *x) {
  scs_int i, j;
  for (i = 0; i < n; i++) {
    scs_float val = x[i];
    for (j = Lp[i]; j < Lp[i + 1]; j++) {
      x[Li[j]] -= Lx[j] * val;
    }
  }
}

/* Backward solve: (L + I)' * x = b, where L is strictly lower-triangular.
 * Overwrites b with the solution. Same as QDLDL_Ltsolve. */
static void backward_solve(scs_int n, const scs_int *Lp, const scs_int *Li,
                           const scs_float *Lx, scs_float *x) {
  scs_int i, j;
  for (i = n - 1; i >= 0; i--) {
    scs_float val = x[i];
    for (j = Lp[i]; j < Lp[i + 1]; j++) {
      val -= Lx[j] * x[Li[j]];
    }
    x[i] = val;
  }
}

/* Block-diagonal solve: D * z = y, where D is tridiagonal with 1x1 and 2x2
 * Bunch-Kaufman blocks. Overwrites y with the solution. */
static void diag_solve(scs_int n, const scs_float *D_diag,
                       const scs_float *D_sub, scs_float *x) {
  scs_int i = 0;
  while (i < n) {
    if (i < n - 1 && D_sub[i] != 0.0) {
      /* 2x2 block: [a b; b d] */
      scs_float a = D_diag[i], b = D_sub[i], d = D_diag[i + 1];
      scs_float det = a * d - b * b;
      scs_float y0 = x[i], y1 = x[i + 1];
      x[i] = (d * y0 - b * y1) / det;
      x[i + 1] = (a * y1 - b * y0) / det;
      i += 2;
    } else {
      /* 1x1 block */
      x[i] /= D_diag[i];
      i += 1;
    }
  }
}

ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r) {
  scs_int n_plus_m = A->n + A->m;
  ScsLinSysWork *p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  if (!p) {
    return SCS_NULL;
  }

  p->n = A->n;
  p->m = A->m;
  p->diag_p = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  p->diag_r_idxs = (scs_int *)scs_calloc(n_plus_m, sizeof(scs_int));
  p->L = SCS_NULL;
  p->D_diag = (scs_float *)scs_calloc(n_plus_m, sizeof(scs_float));
  p->D_sub = (scs_float *)scs_calloc(n_plus_m > 0 ? n_plus_m - 1 : 0,
                                      sizeof(scs_float));
  p->perm = (scs_int *)scs_calloc(n_plus_m, sizeof(scs_int));
  p->bp = (scs_float *)scs_calloc(n_plus_m, sizeof(scs_float));
  p->factorizations = 0;

  if (!p->diag_p || !p->diag_r_idxs || !p->D_diag || (n_plus_m > 1 && !p->D_sub) || !p->perm || !p->bp) {
    scs_printf("Error allocating memory for linear system workspace.\n");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* Form upper-triangular KKT matrix */
  p->kkt = SCS(form_kkt)(A, P, p->diag_p, diag_r, p->diag_r_idxs, 1);
  if (!p->kkt) {
    scs_printf("Error forming KKT matrix.\n");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* Factorize via MATLAB's ldl and cache factors in C */
  if (matlab_ldl_factor(p) < 0) {
    scs_printf("Error in initial LDL factorization.\n");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  return p;
}

/* Solve the KKT system using cached LDL factors.
 * K(p,p) = L*D*L' => K = P'*L*D*L'*P.
 * Solution overwrites b. */
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                          scs_float tol) {
  scs_int n_plus_m = p->n + p->m;
  scs_int i;

  /* Permute: bp = b(perm) */
  for (i = 0; i < n_plus_m; i++) {
    p->bp[i] = b[p->perm[i]];
  }

  /* Forward solve: L * y = bp */
  forward_solve(n_plus_m, p->L->p, p->L->i, p->L->x, p->bp);

  /* Block-diagonal solve: D * z = y */
  diag_solve(n_plus_m, p->D_diag, p->D_sub, p->bp);

  /* Backward solve: L' * x = z */
  backward_solve(n_plus_m, p->L->p, p->L->i, p->L->x, p->bp);

  /* Inverse permute: b(perm) = bp */
  for (i = 0; i < n_plus_m; i++) {
    b[p->perm[i]] = p->bp[i];
  }

  return 0;
}

/* Update diagonal of R in the KKT matrix and refactorize.
 * Only the diagonal values change — the sparsity pattern is unchanged.
 * Despite this, matlab_ldl_factor recomputes the AMD permutation because
 * MATLAB's ldl() does not accept a pre-computed permutation (see comment
 * in matlab_ldl_factor for alternatives). */
scs_int scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
  scs_int i;

  for (i = 0; i < p->n; ++i) {
    /* top left: R_x + P */
    p->kkt->x[p->diag_r_idxs[i]] = p->diag_p[i] + diag_r[i];
  }
  for (i = p->n; i < p->n + p->m; ++i) {
    /* bottom right: -R_y */
    p->kkt->x[p->diag_r_idxs[i]] = -diag_r[i];
  }

  if (matlab_ldl_factor(p) < 0) {
    scs_printf("Error in LDL refactorization.\n");
    return -1;
  }
  return 0;
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    SCS(cs_spfree)(p->L);
    SCS(cs_spfree)(p->kkt);
    scs_free(p->D_diag);
    scs_free(p->D_sub);
    scs_free(p->perm);
    scs_free(p->bp);
    scs_free(p->diag_r_idxs);
    scs_free(p->diag_p);
    scs_free(p);
  }
}
