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
  mxArray *K_sym;        /* full symmetric KKT matrix as persistent mxArray */
  scs_int factorizations;
};

#ifdef __cplusplus
}
#endif
#endif
