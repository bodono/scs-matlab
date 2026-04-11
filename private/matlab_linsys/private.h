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
  /* MATLAB mxArray pointers for LDL factors: A(p,p) = L*D*L' */
  mxArray *L;            /* unit lower triangular factor */
  mxArray *D;            /* block diagonal factor (1x1 and 2x2 blocks) */
  mxArray *perm;         /* permutation vector (1-indexed, MATLAB convention) */
  scs_int factorizations;
};

#ifdef __cplusplus
}
#endif
#endif
