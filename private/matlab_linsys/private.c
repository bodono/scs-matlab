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

/* Perform LDL factorization via MATLAB's ldl(K, 'vector').
 * Stores L, Lt, D and permutation in the workspace.
 * Destroys any prior factors. */
static scs_int matlab_ldl_factor(ScsLinSysWork *p) {
  mxArray *K_mx, *ldl_rhs[2], *ldl_lhs[3];
  mxArray *ct_lhs[1];
  double *perm_pr;
  scs_int n_plus_m = p->n + p->m;
  scs_int i;

  /* Destroy old factors if they exist */
  if (p->L) {
    mxDestroyArray(p->L);
    p->L = SCS_NULL;
  }
  if (p->Lt) {
    mxDestroyArray(p->Lt);
    p->Lt = SCS_NULL;
  }
  if (p->D) {
    mxDestroyArray(p->D);
    p->D = SCS_NULL;
  }

  /* Convert KKT to MATLAB sparse (ldl uses upper triangle only) */
  K_mx = scs_to_mxsparse(p->kkt);

  /* Call [L, D, perm_vec] = ldl(K, 'vector') */
  ldl_rhs[0] = K_mx;
  ldl_rhs[1] = mxCreateString("vector");
  if (mexCallMATLAB(3, ldl_lhs, 2, ldl_rhs, "ldl") != 0) {
    scs_printf("Error calling MATLAB ldl().\n");
    mxDestroyArray(K_mx);
    mxDestroyArray(ldl_rhs[1]);
    return -1;
  }

  /* Store L and D as persistent mxArrays */
  p->L = ldl_lhs[0];
  p->D = ldl_lhs[1];
  mexMakeArrayPersistent(p->L);
  mexMakeArrayPersistent(p->D);

  /* Compute and store Lt = L' (persistent) */
  if (mexCallMATLAB(1, ct_lhs, 1, &p->L, "ctranspose") != 0) {
    scs_printf("Error computing L transpose.\n");
    mxDestroyArray(K_mx);
    mxDestroyArray(ldl_rhs[1]);
    mxDestroyArray(ldl_lhs[2]);
    return -1;
  }
  p->Lt = ct_lhs[0];
  mexMakeArrayPersistent(p->Lt);

  /* Extract permutation vector and convert to 0-indexed C array */
  perm_pr = mxGetPr(ldl_lhs[2]);
  for (i = 0; i < n_plus_m; i++) {
    p->perm[i] = (scs_int)(perm_pr[i] - 1); /* MATLAB 1-indexed -> C 0-indexed */
  }
  mxDestroyArray(ldl_lhs[2]);

  /* Cleanup temporaries */
  mxDestroyArray(K_mx);
  mxDestroyArray(ldl_rhs[1]);

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
  p->perm = (scs_int *)scs_calloc(n_plus_m, sizeof(scs_int));
  p->L = SCS_NULL;
  p->Lt = SCS_NULL;
  p->D = SCS_NULL;
  p->factorizations = 0;

  /* Form upper-triangular KKT matrix */
  p->kkt = SCS(form_kkt)(A, P, p->diag_p, diag_r, p->diag_r_idxs, 1);
  if (!p->kkt) {
    scs_printf("Error forming KKT matrix.\n");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* Factorize via MATLAB's ldl */
  if (matlab_ldl_factor(p) < 0) {
    scs_printf("Error in MATLAB LDL initial factorization.\n");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  return p;
}

/* Solves the KKT system using the stored LDL factors.
 * Performs: permute b, then y = L\bp, z = D\y, w = Lt\z, unpermute.
 * Solution overwrites b. */
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                          scs_float tol) {
  scs_int n_plus_m = p->n + p->m;
  mxArray *bp_mx, *rhs[2], *y_mx, *z_mx, *w_mx;
  double *bp_pr, *w_pr;
  scs_int i;

  /* Create permuted b: bp[i] = b[perm[i]] */
  bp_mx = mxCreateDoubleMatrix((mwSize)n_plus_m, 1, mxREAL);
  bp_pr = mxGetPr(bp_mx);
  for (i = 0; i < n_plus_m; i++) {
    bp_pr[i] = (double)b[p->perm[i]];
  }

  /* y = L \ bp */
  rhs[0] = p->L;
  rhs[1] = bp_mx;
  if (mexCallMATLAB(1, &y_mx, 2, rhs, "mldivide") != 0) {
    scs_printf("Error in L \\ bp.\n");
    mxDestroyArray(bp_mx);
    return -1;
  }
  mxDestroyArray(bp_mx);

  /* z = D \ y */
  rhs[0] = p->D;
  rhs[1] = y_mx;
  if (mexCallMATLAB(1, &z_mx, 2, rhs, "mldivide") != 0) {
    scs_printf("Error in D \\ y.\n");
    mxDestroyArray(y_mx);
    return -1;
  }
  mxDestroyArray(y_mx);

  /* w = Lt \ z  (Lt = L') */
  rhs[0] = p->Lt;
  rhs[1] = z_mx;
  if (mexCallMATLAB(1, &w_mx, 2, rhs, "mldivide") != 0) {
    scs_printf("Error in Lt \\ z.\n");
    mxDestroyArray(z_mx);
    return -1;
  }
  mxDestroyArray(z_mx);

  /* Inverse-permute w back into b: b[perm[i]] = w[i] */
  w_pr = mxGetPr(w_mx);
  for (i = 0; i < n_plus_m; i++) {
    b[p->perm[i]] = (scs_float)w_pr[i];
  }
  mxDestroyArray(w_mx);

  return 0;
}

/* Update diagonal of R in the KKT matrix and re-factorize */
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

  if (matlab_ldl_factor(p) < 0) {
    scs_printf("Error in MATLAB LDL factorization when updating.\n");
    return -1;
  }
  return 0;
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    if (p->L) {
      mxDestroyArray(p->L);
    }
    if (p->Lt) {
      mxDestroyArray(p->Lt);
    }
    if (p->D) {
      mxDestroyArray(p->D);
    }
    SCS(cs_spfree)(p->kkt);
    scs_free(p->diag_r_idxs);
    scs_free(p->diag_p);
    scs_free(p->perm);
    scs_free(p);
  }
}
