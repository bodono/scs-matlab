#include "private.h"

const char *scs_get_lin_sys_method(void) {
  return "sparse-direct-matlab-ldl";
}

/* Convert ScsMatrix (CSC, upper triangular) to MATLAB sparse mxArray */
static mxArray *scs_to_mxsparse(const ScsMatrix *M) {
  scs_int j, nnz;
  mxArray *mx;
  double *pr;
  mwIndex *ir, *jc;

  nnz = M->p[M->n];
  mx = mxCreateSparse((mwSize)M->m, (mwSize)M->n, (mwSize)nnz, mxREAL);
  pr = mxGetPr(mx);
  ir = mxGetIr(mx);
  jc = mxGetJc(mx);

  for (j = 0; j < nnz; j++) {
    pr[j] = (double)M->x[j];
    ir[j] = (mwIndex)M->i[j];
  }
  for (j = 0; j <= M->n; j++) {
    jc[j] = (mwIndex)M->p[j];
  }
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
  p->L->i = (scs_int *)scs_calloc(nnz_nodiag, sizeof(scs_int));
  p->L->x = (scs_float *)scs_calloc(nnz_nodiag, sizeof(scs_float));

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

/* Extract inverse of diagonal D from MATLAB sparse mxArray.
 * With ldl(K, 0, 'vector'), D is guaranteed to be diagonal (no 2x2 blocks). */
static scs_int extract_Dinv(ScsLinSysWork *p, const mxArray *D_mx) {
  scs_int n_plus_m = p->n + p->m;
  mwIndex *jc = mxGetJc(D_mx);
  mwIndex *ir = mxGetIr(D_mx);
  double *pr = mxGetPr(D_mx);
  scs_int j;
  mwIndex k;

  for (j = 0; j < n_plus_m; j++) {
    p->Dinv[j] = 0.0;
    for (k = jc[j]; k < jc[j + 1]; k++) {
      if ((scs_int)ir[k] == j) {
        p->Dinv[j] = (scs_float)(1.0 / pr[k]);
        break;
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
 * Calls [L, D, perm] = ldl(K, 0, 'vector') and extracts factors into C.
 * The upper-triangular KKT is passed directly (ldl only reads upper triangle).
 * thresh=0 disables Bunch-Kaufman pivoting so D is purely diagonal. */
static scs_int matlab_ldl_factor(ScsLinSysWork *p) {
  mxArray *K_mx, *rhs[3], *lhs[3];

  /* Convert C KKT (upper triangular CSC) to MATLAB sparse */
  K_mx = scs_to_mxsparse(p->kkt);

  /* [L, D, perm] = ldl(K, 0, 'vector') */
  rhs[0] = K_mx;
  rhs[1] = mxCreateDoubleScalar(0.0);
  rhs[2] = mxCreateString("vector");

  if (mexCallMATLAB(3, lhs, 3, rhs, "ldl") != 0) {
    scs_printf("Error in MATLAB ldl() factorization.\n");
    mxDestroyArray(K_mx);
    mxDestroyArray(rhs[1]);
    mxDestroyArray(rhs[2]);
    return -1;
  }

  /* Extract factors into C arrays */
  extract_L(p, lhs[0]);
  extract_Dinv(p, lhs[1]);
  extract_perm(p, lhs[2]);

  /* Clean up all MATLAB temporaries */
  mxDestroyArray(K_mx);
  mxDestroyArray(rhs[1]);
  mxDestroyArray(rhs[2]);
  mxDestroyArray(lhs[0]);
  mxDestroyArray(lhs[1]);
  mxDestroyArray(lhs[2]);

  p->factorizations++;
  return 0;
}

/* Forward solve: (L + I) * x = b, where L is strictly lower-triangular.
 * Overwrites b with the solution. */
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
 * Overwrites b with the solution. */
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
  p->Dinv = (scs_float *)scs_calloc(n_plus_m, sizeof(scs_float));
  p->perm = (scs_int *)scs_calloc(n_plus_m, sizeof(scs_int));
  p->bp = (scs_float *)scs_calloc(n_plus_m, sizeof(scs_float));
  p->factorizations = 0;

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
 * Solves P'*L*D*L'*P * x = b where P is the fill-reducing permutation.
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

  /* Diagonal solve: D * z = y */
  for (i = 0; i < n_plus_m; i++) {
    p->bp[i] *= p->Dinv[i];
  }

  /* Backward solve: L' * x = z */
  backward_solve(n_plus_m, p->L->p, p->L->i, p->L->x, p->bp);

  /* Inverse permute: b(perm) = bp */
  for (i = 0; i < n_plus_m; i++) {
    b[p->perm[i]] = p->bp[i];
  }

  return 0;
}

/* Update diagonal of R in the KKT matrix and refactorize. */
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
    scs_free(p->Dinv);
    scs_free(p->perm);
    scs_free(p->bp);
    scs_free(p->diag_r_idxs);
    scs_free(p->diag_p);
    scs_free(p);
  }
}
