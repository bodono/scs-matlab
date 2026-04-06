#include "glbopts.h"
#include "linalg.h"
#include "matrix.h"
#include "mex.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

#include <string.h>

void free_mex(ScsData *d, ScsCone *k, ScsSettings *stgs);

/* ======================== Workspace state ======================== */
static ScsWork *ws_work = SCS_NULL;
static scs_int ws_n = 0;
static scs_int ws_m = 0;

static void ws_cleanup(void) {
  if (ws_work) {
    scs_finish(ws_work);
    ws_work = SCS_NULL;
    ws_n = 0;
    ws_m = 0;
  }
}

/* ======================== Helper functions ======================== */

scs_int parse_warm_start(const mxArray *p_mex, scs_float **p, scs_int l) {
  *p = (scs_float *)scs_calloc(
      l, sizeof(scs_float)); /* this allocates memory used for ScsSolution */
  if (p_mex == SCS_NULL) {
    return 0;
  } else if (mxIsSparse(p_mex) ||
             (scs_int)mxGetNumberOfElements(p_mex) != l) {
    scs_printf("Error parsing warm start input (make sure vectors are not "
               "sparse and of correct size), running without full "
               "warm-start");
    return 0;
  } else {
    memcpy(*p, mxGetPr(p_mex), l * sizeof(scs_float));
    return 1;
  }
}

#if !(DLONG > 0)
/* this memory must be freed */
scs_int *cast_to_scs_int_arr(mwIndex *arr, scs_int len) {
  scs_int i;
  scs_int *arr_out = (scs_int *)scs_malloc(sizeof(scs_int) * len);
  for (i = 0; i < len; i++) {
    arr_out[i] = (scs_int)arr[i];
  }
  return arr_out;
}
#endif

#if SFLOAT > 0
/* this memory must be freed */
scs_float *cast_to_scs_float_arr(double *arr, scs_int len) {
  scs_int i;
  scs_float *arr_out = (scs_float *)scs_malloc(sizeof(scs_float) * len);
  for (i = 0; i < len; i++) {
    arr_out[i] = (scs_float)arr[i];
  }
  return arr_out;
}

double *cast_to_double_arr(scs_float *arr, scs_int len) {
  scs_int i;
  double *arr_out = (double *)scs_malloc(sizeof(double) * len);
  for (i = 0; i < len; i++) {
    arr_out[i] = (double)arr[i];
  }
  return arr_out;
}
#endif

void set_output_field(mxArray **pout, scs_float *out, scs_int len) {
  *pout = mxCreateDoubleMatrix(0, 0, mxREAL);
#if SFLOAT > 0
  mxSetPr(*pout, cast_to_double_arr(out, len));
  scs_free(out);
#else
  mxSetPr(*pout, out);
#endif
  mxSetM(*pout, len);
  mxSetN(*pout, 1);
}

/* Get length of a MATLAB array (handles row and column vectors) */
static scs_int get_mex_length(const mxArray *arr) {
  scs_int ndim = (scs_int)mxGetNumberOfDimensions(arr);
  const size_t *dims = mxGetDimensions(arr);
  scs_int len = (scs_int)dims[0];
  if (ndim > 1 && dims[0] == 1) {
    len = (scs_int)dims[1];
  }
  return len;
}

/* Parse data struct (A, P, b, c) into ScsData.
 * Caller must free d, d->A, d->P (if non-NULL) via free_mex(d, NULL, NULL). */
static void parse_data(const mxArray *data_mex, ScsData **d_out) {
  ScsData *d;
  ScsMatrix *A;
  ScsMatrix *P = SCS_NULL;
  const mxArray *A_mex, *P_mex, *b_mex, *c_mex;

  d = (ScsData *)scs_malloc(sizeof(ScsData));

  A_mex = (mxArray *)mxGetField(data_mex, 0, "A");
  if (A_mex == SCS_NULL) {
    scs_free(d);
    mexErrMsgTxt("ScsData struct must contain a `A` entry.");
  }
  if (!mxIsSparse(A_mex)) {
    scs_free(d);
    mexErrMsgTxt("Input matrix A must be in sparse format (pass in sparse(A))");
  }
  P_mex = (mxArray *)mxGetField(data_mex, 0, "P"); /* can be SCS_NULL */
  if (P_mex && !mxIsSparse(P_mex)) {
    scs_free(d);
    mexErrMsgTxt("Input matrix P must be in sparse format (pass in sparse(P))");
  }
  b_mex = (mxArray *)mxGetField(data_mex, 0, "b");
  if (b_mex == SCS_NULL) {
    scs_free(d);
    mexErrMsgTxt("ScsData struct must contain a `b` entry.");
  }
  if (mxIsSparse(b_mex)) {
    scs_free(d);
    mexErrMsgTxt("Input vector b must be in dense format (pass in full(b))");
  }
  c_mex = (mxArray *)mxGetField(data_mex, 0, "c");
  if (c_mex == SCS_NULL) {
    scs_free(d);
    mexErrMsgTxt("ScsData struct must contain a `c` entry.");
  }
  if (mxIsSparse(c_mex)) {
    scs_free(d);
    mexErrMsgTxt("Input vector c must be in dense format (pass in full(c))");
  }

  d->n = (scs_int) * (mxGetDimensions(c_mex));
  d->m = (scs_int) * (mxGetDimensions(b_mex));
#if SFLOAT > 0
  d->b = cast_to_scs_float_arr(mxGetPr(b_mex), d->m);
  d->c = cast_to_scs_float_arr(mxGetPr(c_mex), d->n);
#else
  d->b = (scs_float *)mxGetPr(b_mex);
  d->c = (scs_float *)mxGetPr(c_mex);
#endif

  A = (ScsMatrix *)scs_malloc(sizeof(ScsMatrix));
  A->n = d->n;
  A->m = d->m;

  if (P_mex) {
    P = (ScsMatrix *)scs_malloc(sizeof(ScsMatrix));
    P->n = d->n;
    P->m = d->n;
  }

#if DLONG > 0
  A->p = (scs_int *)mxGetJc(A_mex);
  A->i = (scs_int *)mxGetIr(A_mex);
  if (P_mex) {
    P->p = (scs_int *)mxGetJc(P_mex);
    P->i = (scs_int *)mxGetIr(P_mex);
  }
#else
  A->p = cast_to_scs_int_arr(mxGetJc(A_mex), A->n + 1);
  A->i = cast_to_scs_int_arr(mxGetIr(A_mex), A->p[A->n]);
  if (P_mex) {
    P->p = cast_to_scs_int_arr(mxGetJc(P_mex), P->n + 1);
    P->i = cast_to_scs_int_arr(mxGetIr(P_mex), P->p[P->n]);
  }
#endif
#if SFLOAT > 0
  A->x = cast_to_scs_float_arr(mxGetPr(A_mex), A->p[A->n]);
  if (P_mex) {
    P->x = cast_to_scs_float_arr(mxGetPr(P_mex), P->p[P->n]);
  }
#else
  A->x = (scs_float *)mxGetPr(A_mex);
  if (P_mex) {
    P->x = (scs_float *)mxGetPr(P_mex);
  }
#endif
  d->A = A;
  d->P = P;

  *d_out = d;
}

/* Parse settings struct into ScsSettings.
 * Caller must free via free_mex(NULL, NULL, stgs). */
static void parse_settings(const mxArray *settings_mex, ScsSettings **stgs_out) {
  ScsSettings *stgs;
  mxArray *tmp;

  stgs = (ScsSettings *)scs_malloc(sizeof(ScsSettings));
  scs_set_default_settings(stgs);

  tmp = mxGetField(settings_mex, 0, "alpha");
  if (tmp != SCS_NULL) {
    stgs->alpha = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "rho_x");
  if (tmp != SCS_NULL) {
    stgs->rho_x = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "max_iters");
  if (tmp != SCS_NULL) {
    stgs->max_iters = (scs_int)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "scale");
  if (tmp != SCS_NULL) {
    stgs->scale = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "eps_abs");
  if (tmp != SCS_NULL) {
    stgs->eps_abs = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "eps_rel");
  if (tmp != SCS_NULL) {
    stgs->eps_rel = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "eps_infeas");
  if (tmp != SCS_NULL) {
    stgs->eps_infeas = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "verbose");
  if (tmp != SCS_NULL) {
    stgs->verbose = (scs_int)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "normalize");
  if (tmp != SCS_NULL) {
    stgs->normalize = (scs_int)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "acceleration_lookback");
  if (tmp != SCS_NULL) {
    stgs->acceleration_lookback = (scs_int)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "acceleration_interval");
  if (tmp != SCS_NULL) {
    stgs->acceleration_interval = (scs_int)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "adaptive_scale");
  if (tmp != SCS_NULL) {
    stgs->adaptive_scale = (scs_int)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "time_limit_secs");
  if (tmp != SCS_NULL) {
    stgs->time_limit_secs = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "write_data_filename");
  if (tmp != SCS_NULL) {
    /* need to free this later */
    stgs->write_data_filename = mxArrayToString(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "log_csv_filename");
  if (tmp != SCS_NULL) {
    /* need to free this later */
    stgs->log_csv_filename = mxArrayToString(tmp);
  }

  *stgs_out = stgs;
}

/* Parse cone struct into ScsCone.
 * Caller must free via free_mex(NULL, k, NULL). */
static void parse_cones(const mxArray *cone_mex, ScsCone **k_out) {
  ScsCone *k;
  scs_int i, ns, ncs, nbl, nbu, blen;
  const mxArray *kf, *kz, *kl, *kbl, *kbu, *kq, *ks, *kcs, *kep, *ked, *kp;
  const double *q_mex, *bl_mex, *bu_mex, *s_mex, *cs_mex, *p_mex;
  const size_t *bl_dims, *bu_dims, *q_dims, *s_dims, *cs_dims, *p_dims;

  k = (ScsCone *)scs_malloc(sizeof(ScsCone));

  /* TODO rm this */
  kf = mxGetField(cone_mex, 0, "f");
  if (kf && !mxIsEmpty(kf)) {
    scs_printf("SCS deprecation warning: The 'f' field in the cone struct \n"
               "has been replaced by 'z' to better reflect the Zero cone. \n"
               "Please replace usage of 'f' with 'z'. If both 'f' and 'z' \n"
               "are set then we sum the two fields to get the final zero \n"
               "cone size.\n");
    k->z = (scs_int)*mxGetPr(kf);
  } else {
    k->z = 0;
  }

  kz = mxGetField(cone_mex, 0, "z");
  if (kz && !mxIsEmpty(kz)) {
    k->z += (scs_int)*mxGetPr(kz); /* TODO rm this */
  }

  kl = mxGetField(cone_mex, 0, "l");
  if (kl && !mxIsEmpty(kl)) {
    k->l = (scs_int)*mxGetPr(kl);
  } else {
    k->l = 0;
  }

  kep = mxGetField(cone_mex, 0, "ep");
  if (kep && !mxIsEmpty(kep)) {
    k->ep = (scs_int)*mxGetPr(kep);
  } else {
    k->ep = 0;
  }

  ked = mxGetField(cone_mex, 0, "ed");
  if (ked && !mxIsEmpty(ked)) {
    k->ed = (scs_int)*mxGetPr(ked);
  } else {
    k->ed = 0;
  }

  kbu = mxGetField(cone_mex, 0, "bu");
  kbl = mxGetField(cone_mex, 0, "bl");
  if (kbl && kbu && !mxIsEmpty(kbl) && !mxIsEmpty(kbu)) {
    bl_mex = mxGetPr(kbl);
    bu_mex = mxGetPr(kbu);

    nbl = (scs_int)mxGetNumberOfDimensions(kbl);
    nbu = (scs_int)mxGetNumberOfDimensions(kbu);
    if (nbl != nbu) {
      mexErrMsgTxt("bl,bu cone entries not the same size.");
    }

    bl_dims = mxGetDimensions(kbl);
    bu_dims = mxGetDimensions(kbu);
    for (i = 0; i < nbu; i++) {
      if (bl_dims[i] != bu_dims[i]) {
        mexErrMsgTxt("bl,bu cone entries not the same size.");
      }
    }
    blen = (scs_int)bu_dims[0];
    if (nbl > 1 && bu_dims[0] == 1) {
      blen = (scs_int)bu_dims[1];
    }
    k->bu = (scs_float *)scs_malloc(sizeof(scs_float) * blen);
    k->bl = (scs_float *)scs_malloc(sizeof(scs_float) * blen);

    for (i = 0; i < blen; i++) {
      k->bl[i] = (scs_float)bl_mex[i];
      k->bu[i] = (scs_float)bu_mex[i];
    }
    k->bsize = blen + 1; /* bsize is total length of cone (t,s) */
  } else {
    k->bsize = 0;
    k->bl = SCS_NULL;
    k->bu = SCS_NULL;
  }

  kq = mxGetField(cone_mex, 0, "q");
  if (kq && !mxIsEmpty(kq)) {
    q_mex = mxGetPr(kq);
    k->qsize = get_mex_length(kq);
    k->q = (scs_int *)scs_malloc(sizeof(scs_int) * k->qsize);
    for (i = 0; i < k->qsize; i++) {
      k->q[i] = (scs_int)q_mex[i];
    }
  } else {
    k->qsize = 0;
    k->q = SCS_NULL;
  }

  ks = mxGetField(cone_mex, 0, "s");
  if (ks && !mxIsEmpty(ks)) {
    s_mex = mxGetPr(ks);
    k->ssize = get_mex_length(ks);
    k->s = (scs_int *)scs_malloc(sizeof(scs_int) * k->ssize);
    for (i = 0; i < k->ssize; i++) {
      k->s[i] = (scs_int)s_mex[i];
    }
  } else {
    k->ssize = 0;
    k->s = SCS_NULL;
  }

  kcs = mxGetField(cone_mex, 0, "cs");
  if (kcs && !mxIsEmpty(kcs)) {
    cs_mex = mxGetPr(kcs);
    k->cssize = get_mex_length(kcs);
    k->cs = (scs_int *)scs_malloc(sizeof(scs_int) * k->cssize);
    for (i = 0; i < k->cssize; i++) {
      k->cs[i] = (scs_int)cs_mex[i];
    }
  } else {
    k->cssize = 0;
    k->cs = SCS_NULL;
  }

  kp = mxGetField(cone_mex, 0, "p");
  if (kp && !mxIsEmpty(kp)) {
    p_mex = mxGetPr(kp);
    k->psize = get_mex_length(kp);
    k->p = (scs_float *)scs_malloc(sizeof(scs_float) * k->psize);
    for (i = 0; i < k->psize; i++) {
      k->p[i] = (scs_float)p_mex[i];
    }
  } else {
    k->psize = 0;
    k->p = SCS_NULL;
  }

#ifdef USE_SPECTRAL_CONES
  {
    const mxArray *kd, *knuc_m, *knuc_n, *kell1, *ksl_n, *ksl_k;
    const double *d_mex_arr, *nuc_m_mex, *nuc_n_mex, *ell1_mex;
    const double *sl_n_mex, *sl_k_mex;

    kd = mxGetField(cone_mex, 0, "d");
    if (kd && !mxIsEmpty(kd)) {
      d_mex_arr = mxGetPr(kd);
      k->dsize = get_mex_length(kd);
      k->d = (scs_int *)scs_malloc(sizeof(scs_int) * k->dsize);
      for (i = 0; i < k->dsize; i++) {
        k->d[i] = (scs_int)d_mex_arr[i];
      }
    } else {
      k->dsize = 0;
      k->d = SCS_NULL;
    }

    knuc_m = mxGetField(cone_mex, 0, "nuc_m");
    knuc_n = mxGetField(cone_mex, 0, "nuc_n");
    if (knuc_m && knuc_n && !mxIsEmpty(knuc_m) && !mxIsEmpty(knuc_n)) {
      nuc_m_mex = mxGetPr(knuc_m);
      nuc_n_mex = mxGetPr(knuc_n);
      k->nucsize = get_mex_length(knuc_m);
      k->nuc_m = (scs_int *)scs_malloc(sizeof(scs_int) * k->nucsize);
      k->nuc_n = (scs_int *)scs_malloc(sizeof(scs_int) * k->nucsize);
      for (i = 0; i < k->nucsize; i++) {
        k->nuc_m[i] = (scs_int)nuc_m_mex[i];
        k->nuc_n[i] = (scs_int)nuc_n_mex[i];
      }
    } else {
      k->nucsize = 0;
      k->nuc_m = SCS_NULL;
      k->nuc_n = SCS_NULL;
    }

    kell1 = mxGetField(cone_mex, 0, "ell1");
    if (kell1 && !mxIsEmpty(kell1)) {
      ell1_mex = mxGetPr(kell1);
      k->ell1_size = get_mex_length(kell1);
      k->ell1 = (scs_int *)scs_malloc(sizeof(scs_int) * k->ell1_size);
      for (i = 0; i < k->ell1_size; i++) {
        k->ell1[i] = (scs_int)ell1_mex[i];
      }
    } else {
      k->ell1_size = 0;
      k->ell1 = SCS_NULL;
    }

    ksl_n = mxGetField(cone_mex, 0, "sl_n");
    ksl_k = mxGetField(cone_mex, 0, "sl_k");
    if (ksl_n && ksl_k && !mxIsEmpty(ksl_n) && !mxIsEmpty(ksl_k)) {
      sl_n_mex = mxGetPr(ksl_n);
      sl_k_mex = mxGetPr(ksl_k);
      k->sl_size = get_mex_length(ksl_n);
      k->sl_n = (scs_int *)scs_malloc(sizeof(scs_int) * k->sl_size);
      k->sl_k = (scs_int *)scs_malloc(sizeof(scs_int) * k->sl_size);
      for (i = 0; i < k->sl_size; i++) {
        k->sl_n[i] = (scs_int)sl_n_mex[i];
        k->sl_k[i] = (scs_int)sl_k_mex[i];
      }
    } else {
      k->sl_size = 0;
      k->sl_n = SCS_NULL;
      k->sl_k = SCS_NULL;
    }
  }
#endif

  *k_out = k;
}

/* Write ScsInfo to a MATLAB struct and assign to plhs[3] */
static void write_info(mxArray **plhs3, const ScsInfo *info) {
  const mwSize one[1] = {1};
  const int num_info_fields = 22;
  const char *info_fields[] = {
      "iter",       "status",         "pobj",          "dobj",
      "res_pri",    "res_dual",       "res_infeas",    "res_unbdd_a",
      "scale",      "status_val",     "res_unbdd_p",   "gap",
      "setup_time", "solve_time",     "scale_updates", "comp_slack",
      "lin_sys_solver", "rejected_accel_steps", "accepted_accel_steps",
      "lin_sys_time",   "cone_time",            "accel_time"};
  mxArray *tmp;

  *plhs3 = mxCreateStructArray(1, one, num_info_fields, info_fields);

  mxSetField(*plhs3, 0, "status", mxCreateString(info->status));

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "iter", tmp);
  *mxGetPr(tmp) = (scs_float)info->iter;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "scale_updates", tmp);
  *mxGetPr(tmp) = (scs_float)info->scale_updates;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "status_val", tmp);
  *mxGetPr(tmp) = (scs_float)info->status_val;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "pobj", tmp);
  *mxGetPr(tmp) = info->pobj;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "dobj", tmp);
  *mxGetPr(tmp) = info->dobj;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "res_pri", tmp);
  *mxGetPr(tmp) = info->res_pri;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "res_dual", tmp);
  *mxGetPr(tmp) = info->res_dual;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "res_infeas", tmp);
  *mxGetPr(tmp) = info->res_infeas;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "res_unbdd_a", tmp);
  *mxGetPr(tmp) = info->res_unbdd_a;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "res_unbdd_p", tmp);
  *mxGetPr(tmp) = info->res_unbdd_p;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "comp_slack", tmp);
  *mxGetPr(tmp) = info->comp_slack;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "gap", tmp);
  *mxGetPr(tmp) = info->gap;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "scale", tmp);
  *mxGetPr(tmp) = info->scale;

  /*info.time is millisecs - return value in secs */
  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "setup_time", tmp);
  *mxGetPr(tmp) = info->setup_time;

  /*info.time is millisecs - return value in secs */
  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "solve_time", tmp);
  *mxGetPr(tmp) = info->solve_time;

  mxSetField(*plhs3, 0, "lin_sys_solver",
             mxCreateString(info->lin_sys_solver));

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "rejected_accel_steps", tmp);
  *mxGetPr(tmp) = (scs_float)info->rejected_accel_steps;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "accepted_accel_steps", tmp);
  *mxGetPr(tmp) = (scs_float)info->accepted_accel_steps;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "lin_sys_time", tmp);
  *mxGetPr(tmp) = info->lin_sys_time;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "cone_time", tmp);
  *mxGetPr(tmp) = info->cone_time;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(*plhs3, 0, "accel_time", tmp);
  *mxGetPr(tmp) = info->accel_time;
}

/* ======================== MEX entry point ======================== */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
#if VERBOSITY > 0
  scs_printf("SIZE OF mwSize = %i\n", (int)sizeof(mwSize));
  scs_printf("SIZE OF mwIndex = %i\n", (int)sizeof(mwIndex));
#endif

  /* ---- Workspace command dispatch ---- */
  if (nrhs >= 1 && mxIsChar(prhs[0])) {
    char *cmd = mxArrayToString(prhs[0]);

    if (strcmp(cmd, "init") == 0) {
      /* scs_xxx('init', data, cone, settings) */
      ScsData *d;
      ScsCone *k;
      ScsSettings *stgs;
      if (nrhs != 4) {
        scs_free(cmd);
        mexErrMsgTxt("Usage: scs_xxx('init', data, cone, settings)");
      }
      ws_cleanup(); /* free any existing workspace */
      parse_data(prhs[1], &d);
      parse_cones(prhs[2], &k);
      parse_settings(prhs[3], &stgs);

      ws_n = d->n;
      ws_m = d->m;
      ws_work = scs_init(d, k, stgs);

      free_mex(d, k, stgs);

      if (!ws_work) {
        ws_n = 0;
        ws_m = 0;
        mexErrMsgTxt("SCS init failed.");
      }
      mexAtExit(ws_cleanup);
      scs_free(cmd);
      return;
    }

    if (strcmp(cmd, "update") == 0) {
      /* scs_xxx('update', b_new, c_new)
       * Either argument can be [] to leave unchanged. */
      scs_float *b_new = SCS_NULL;
      scs_float *c_new = SCS_NULL;
      if (!ws_work) {
        scs_free(cmd);
        mexErrMsgTxt("No workspace. Call scs_init first.");
      }
      if (nrhs >= 2 && !mxIsEmpty(prhs[1])) {
#if SFLOAT > 0
        b_new = cast_to_scs_float_arr(mxGetPr(prhs[1]), ws_m);
#else
        b_new = (scs_float *)mxGetPr(prhs[1]);
#endif
      }
      if (nrhs >= 3 && !mxIsEmpty(prhs[2])) {
#if SFLOAT > 0
        c_new = cast_to_scs_float_arr(mxGetPr(prhs[2]), ws_n);
#else
        c_new = (scs_float *)mxGetPr(prhs[2]);
#endif
      }
      scs_update(ws_work, b_new, c_new);
#if SFLOAT > 0
      if (b_new) scs_free(b_new);
      if (c_new) scs_free(c_new);
#endif
      scs_free(cmd);
      return;
    }

    if (strcmp(cmd, "solve") == 0) {
      /* [x,y,s,info] = scs_xxx('solve')
       * [x,y,s,info] = scs_xxx('solve', warm_start_struct) */
      ScsSolution sol = {0};
      ScsInfo info;
      scs_int warm_start = 0;
      if (!ws_work) {
        scs_free(cmd);
        mexErrMsgTxt("No workspace. Call scs_init first.");
      }
      if (nrhs >= 2 && !mxIsEmpty(prhs[1])) {
        const mxArray *ws_data = prhs[1];
        warm_start =
            parse_warm_start(mxGetField(ws_data, 0, "x"), &(sol.x), ws_n);
        warm_start |=
            parse_warm_start(mxGetField(ws_data, 0, "y"), &(sol.y), ws_m);
        warm_start |=
            parse_warm_start(mxGetField(ws_data, 0, "s"), &(sol.s), ws_m);
      }
      if (!sol.x)
        sol.x = (scs_float *)scs_calloc(ws_n, sizeof(scs_float));
      if (!sol.y)
        sol.y = (scs_float *)scs_calloc(ws_m, sizeof(scs_float));
      if (!sol.s)
        sol.s = (scs_float *)scs_calloc(ws_m, sizeof(scs_float));

      scs_solve(ws_work, &sol, &info, warm_start);

      set_output_field(&plhs[0], sol.x, ws_n);
      set_output_field(&plhs[1], sol.y, ws_m);
      set_output_field(&plhs[2], sol.s, ws_m);
      write_info(&plhs[3], &info);

      scs_free(cmd);
      return;
    }

    if (strcmp(cmd, "finish") == 0) {
      ws_cleanup();
      scs_free(cmd);
      return;
    }

    scs_free(cmd);
    mexErrMsgTxt("Unknown command. Use 'init', 'solve', 'update', or "
                 "'finish'.");
    return;
  }

  /* ---- One-shot solve: [x,y,s,info] = scs(data,cone,settings) ---- */
  {
    ScsData *d;
    ScsCone *k;
    ScsSettings *stgs;
    ScsSolution sol = {0};
    ScsInfo info;

    if (nrhs != 3) {
      mexErrMsgTxt("Three arguments are required in this order: data struct, "
                   "cone struct, settings struct");
    }
    if (nlhs > 4) {
      mexErrMsgTxt("scs returns up to 4 output arguments only.");
    }

    parse_data(prhs[0], &d);
    parse_cones(prhs[1], &k);
    parse_settings(prhs[2], &stgs);

    /* warm-start */
    stgs->warm_start =
        parse_warm_start(mxGetField(prhs[0], 0, "x"), &(sol.x), d->n);
    stgs->warm_start |=
        parse_warm_start(mxGetField(prhs[0], 0, "y"), &(sol.y), d->m);
    stgs->warm_start |=
        parse_warm_start(mxGetField(prhs[0], 0, "s"), &(sol.s), d->m);

    scs(d, k, stgs, &sol, &info);

    set_output_field(&plhs[0], sol.x, d->n);
    set_output_field(&plhs[1], sol.y, d->m);
    set_output_field(&plhs[2], sol.s, d->m);
    write_info(&plhs[3], &info);

    free_mex(d, k, stgs);
  }
}

void free_mex(ScsData *d, ScsCone *k, ScsSettings *stgs) {
  if (k) {
    if (k->q) {
      scs_free(k->q);
    }
    if (k->bl) {
      scs_free(k->bl);
    }
    if (k->bu) {
      scs_free(k->bu);
    }
    if (k->s) {
      scs_free(k->s);
    }
    if (k->cs) {
      scs_free(k->cs);
    }
    if (k->p) {
      scs_free(k->p);
    }
#ifdef USE_SPECTRAL_CONES
    if (k->d) {
      scs_free(k->d);
    }
    if (k->nuc_m) {
      scs_free(k->nuc_m);
    }
    if (k->nuc_n) {
      scs_free(k->nuc_n);
    }
    if (k->ell1) {
      scs_free(k->ell1);
    }
    if (k->sl_n) {
      scs_free(k->sl_n);
    }
    if (k->sl_k) {
      scs_free(k->sl_k);
    }
#endif
    scs_free(k);
  }
  if (stgs) {
    if (stgs->write_data_filename) {
      scs_free((void *)stgs->write_data_filename);
    }
    if (stgs->log_csv_filename) {
      scs_free((void *)stgs->log_csv_filename);
    }
    scs_free(stgs);
  }
  if (d) {
#if SFLOAT > 0 /* only free if copies, which is only when flags set */
    if (d->b) {
      scs_free(d->b);
    }
    if (d->c) {
      scs_free(d->c);
    }
#endif
    if (d->A) {
#if !(DLONG > 0) /* only free if copies, which is only when flags set */
      if (d->A->p) {
        scs_free(d->A->p);
      }
      if (d->A->i) {
        scs_free(d->A->i);
      }
#endif
#if SFLOAT > 0 /* only free if copies, which is only when flags set */
      if (d->A->x) {
        scs_free(d->A->x);
      }
#endif
      scs_free(d->A);
    }
    if (d->P) {
#if !(DLONG > 0) /* only free if copies, which is only when flags set */
      if (d->P->p) {
        scs_free(d->P->p);
      }
      if (d->P->i) {
        scs_free(d->P->i);
      }
#endif
#if SFLOAT > 0 /* only free if copies, which is only when flags set */
      if (d->P->x) {
        scs_free(d->P->x);
      }
#endif
      scs_free(d->P);
    }
    scs_free(d);
  }
}
