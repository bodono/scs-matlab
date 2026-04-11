#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "csparse.h"
#include "glbopts.h"
#include "linsys.h"
#include "scs_matrix.h"

#include "mex.h"
#include "matrix.h"

struct SCS_LIN_SYS_WORK {
  scs_int m, n;
  ScsMatrix *kkt;        /* KKT matrix in CSC format (upper triangular) */
  scs_int *diag_r_idxs;  /* indices of R diagonal entries in kkt->x */
  scs_float *diag_p;     /* diagonal of P (objective matrix) */

  /* Cached LDL factors from MATLAB's ldl() */
  ScsMatrix *L;          /* Lower triangular L factor (unit diagonal NOT stored) */
  scs_float *Dinv;       /* Inverse of diagonal D */
  scs_int *perm;         /* Fill-reducing permutation (0-indexed) */
  scs_float *bp;         /* Workspace for permuted RHS */

  scs_int factorizations;
};

#ifdef __cplusplus
}
#endif
#endif
