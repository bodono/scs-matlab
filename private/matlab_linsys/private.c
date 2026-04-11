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

/* Symmetrize an upper-triangular sparse matrix: K_sym = K + triu(K,1)'.
 * Returns a new mxArray (caller must destroy). Destroys K_upper. */
static mxArray *symmetrize_upper(mxArray *K_upper) {
  mxArray *triu_rhs[2], *K_strict, *K_strict_t, *plus_rhs[2], *K_sym;

  /* K_strict = triu(K_upper, 1)  — strict upper triangle */
  triu_rhs[0] = K_upper;
  triu_rhs[1] = mxCreateDoubleScalar(1.0);
  mexCallMATLAB(1, &K_strict, 2, triu_rhs, "triu");
  mxDestroyArray(triu_rhs[1]);

  /* K_strict_t = K_strict' */
  mexCallMATLAB(1, &K_strict_t, 1, &K_strict, "ctranspose");
  mxDestroyArray(K_strict);

  /* K_sym = K_upper + K_strict_t */
  plus_rhs[0] = K_upper;
  plus_rhs[1] = K_strict_t;
  mexCallMATLAB(1, &K_sym, 2, plus_rhs, "plus");
  mxDestroyArray(K_strict_t);
  mxDestroyArray(K_upper);

  return K_sym;
}

/* Build full symmetric KKT matrix and store as persistent mxArray.
 * Destroys any prior stored matrix. */
static scs_int matlab_build_kkt(ScsLinSysWork *p) {
  mxArray *K_upper;

  if (p->K_sym) {
    mxDestroyArray(p->K_sym);
    p->K_sym = SCS_NULL;
  }

  /* Convert C-level KKT (upper triangle CSC) to MATLAB sparse */
  K_upper = scs_to_mxsparse(p->kkt);

  /* Symmetrize: K_sym = K + triu(K,1)' */
  p->K_sym = symmetrize_upper(K_upper);
  /* K_upper is destroyed inside symmetrize_upper */

  mexMakeArrayPersistent(p->K_sym);
  p->factorizations++;
  return 0;
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
  p->K_sym = SCS_NULL;
  p->factorizations = 0;

  /* Form upper-triangular KKT matrix */
  p->kkt = SCS(form_kkt)(A, P, p->diag_p, diag_r, p->diag_r_idxs, 1);
  if (!p->kkt) {
    scs_printf("Error forming KKT matrix.\n");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* Build symmetric MATLAB sparse matrix */
  if (matlab_build_kkt(p) < 0) {
    scs_printf("Error building MATLAB KKT matrix.\n");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  return p;
}

/* Solves the KKT system using MATLAB's backslash: x = K_sym \ b.
 * Solution overwrites b. */
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                          scs_float tol) {
  scs_int n_plus_m = p->n + p->m;
  mxArray *b_mx, *rhs[2], *x_mx;
  double *b_pr, *x_pr;
  scs_int i;

  /* Create MATLAB dense vector from b */
  b_mx = mxCreateDoubleMatrix((mwSize)n_plus_m, 1, mxREAL);
  b_pr = mxGetPr(b_mx);
  for (i = 0; i < n_plus_m; i++) {
    b_pr[i] = (double)b[i];
  }

  /* x = K_sym \ b */
  rhs[0] = p->K_sym;
  rhs[1] = b_mx;
  if (mexCallMATLAB(1, &x_mx, 2, rhs, "mldivide") != 0) {
    scs_printf("Error in K \\ b.\n");
    mxDestroyArray(b_mx);
    return -1;
  }
  mxDestroyArray(b_mx);

  /* Copy result back into b */
  x_pr = mxGetPr(x_mx);
  for (i = 0; i < n_plus_m; i++) {
    b[i] = (scs_float)x_pr[i];
  }
  mxDestroyArray(x_mx);

  return 0;
}

/* Update diagonal of R in the KKT matrix and rebuild MATLAB matrix */
scs_int scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
  scs_int i;

  for (i = 0; i < p->n; ++i) {
    /* top left is R_x + P */
    p->kkt->x[p->diag_r_idxs[i]] = p->diag_p[i] + diag_r[i];
  }
  for (i = p->n; i < p->n + p->m; ++i) {
    /* bottom right is -R_y */
    p->kkt->x[p->diag_r_idxs[i]] = -diag_r[i];
  }

  if (matlab_build_kkt(p) < 0) {
    scs_printf("Error rebuilding MATLAB KKT matrix.\n");
    return -1;
  }
  return 0;
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    if (p->K_sym) {
      mxDestroyArray(p->K_sym);
    }
    SCS(cs_spfree)(p->kkt);
    scs_free(p->diag_r_idxs);
    scs_free(p->diag_p);
    scs_free(p);
  }
}
