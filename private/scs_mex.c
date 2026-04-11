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
  if (!*p) return 0;
  if (p_mex == SCS_NULL) {
    return 0;
  } else if (mxIsSparse(p_mex) ||
             (scs_int)mxGetNumberOfElements(p_mex) != l) {
    scs_printf("Error parsing warm start input (make sure vectors are not "
               "sparse and of correct size), running without full "
               "warm-start\n");
    return 0;
  } else {
    memcpy(*p, mxGetPr(p_mex), l * sizeof(scs_float));
    return 1;
  }
}

#ifndef DLONG
/* this memory must be freed */
scs_int *cast_to_scs_int_arr(mwIndex *arr, scs_int len) {
  scs_int i;
  scs_int *arr_out = (scs_int *)scs_malloc(sizeof(scs_int) * len);
  if (!arr_out) return SCS_NULL;
  for (i = 0; i < len; i++) {
    arr_out[i] = (scs_int)arr[i];
  }
  return arr_out;
}
#endif

#ifdef SFLOAT
/* this memory must be freed */
scs_float *cast_to_scs_float_arr(double *arr, scs_int len) {
  scs_int i;
  scs_float *arr_out = (scs_float *)scs_malloc(sizeof(scs_float) * len);
  if (!arr_out) return SCS_NULL;
  for (i = 0; i < len; i++) {
    arr_out[i] = (scs_float)arr[i];
  }
  return arr_out;
}

double *cast_to_double_arr(scs_float *arr, scs_int len) {
  scs_int i;
  double *arr_out = (double *)scs_malloc(sizeof(double) * len);
  if (!arr_out) return SCS_NULL;
  for (i = 0; i < len; i++) {
    arr_out[i] = (double)arr[i];
  }
  return arr_out;
}
#endif

void set_output_field(mxArray **pout, scs_float *out, scs_int len) {
  scs_int i;
  *pout = mxCreateDoubleMatrix(len, 1, mxREAL);
  if (out != SCS_NULL) {
#ifdef SFLOAT
    double *pr = mxGetPr(*pout);
    for (i = 0; i < len; i++) {
      pr[i] = (double)out[i];
    }
#else
    memcpy(mxGetPr(*pout), out, len * sizeof(scs_float));
#endif
    scs_free(out);
  }
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
static scs_int parse_data(const mxArray *data_mex, ScsData **d_out) {
  ScsData *d;
  ScsMatrix *A;
  ScsMatrix *P = SCS_NULL;
  const mxArray *A_mex, *P_mex, *b_mex, *c_mex;

  d = (ScsData *)scs_calloc(1, sizeof(ScsData));

  A_mex = (mxArray *)mxGetField(data_mex, 0, "A");
  if (A_mex == SCS_NULL) {
    scs_free(d);
    scs_printf("ScsData struct must contain a `A` entry.\n");
    return -1;
  }
  if (!mxIsSparse(A_mex)) {
    scs_free(d);
    scs_printf("Input matrix A must be in sparse format (pass in sparse(A))\n");
    return -1;
  }
  P_mex = (mxArray *)mxGetField(data_mex, 0, "P"); /* can be SCS_NULL */
  if (P_mex && !mxIsSparse(P_mex)) {
    scs_free(d);
    scs_printf("Input matrix P must be in sparse format (pass in sparse(P))\n");
    return -1;
  }
  b_mex = (mxArray *)mxGetField(data_mex, 0, "b");
  if (b_mex == SCS_NULL) {
    scs_free(d);
    scs_printf("ScsData struct must contain a `b` entry.\n");
    return -1;
  }
  if (mxIsSparse(b_mex)) {
    scs_free(d);
    scs_printf("Input vector b must be in dense format (pass in full(b))\n");
    return -1;
  }
  c_mex = (mxArray *)mxGetField(data_mex, 0, "c");
  if (c_mex == SCS_NULL) {
    scs_free(d);
    scs_printf("ScsData struct must contain a `c` entry.\n");
    return -1;
  }
  if (mxIsSparse(c_mex)) {
    scs_free(d);
    scs_printf("Input vector c must be in dense format (pass in full(c))\n");
    return -1;
  }

  d->n = (scs_int) * (mxGetDimensions(c_mex));
  d->m = (scs_int) * (mxGetDimensions(b_mex));
#ifdef SFLOAT
  d->b = cast_to_scs_float_arr(mxGetPr(b_mex), d->m);
  d->c = cast_to_scs_float_arr(mxGetPr(c_mex), d->n);
  if (!d->b || !d->c) {
    free_mex(d, SCS_NULL, SCS_NULL);
    scs_printf("Memory allocation failed for vectors b or c.\n");
    return -1;
  }
#else
  d->b = (scs_float *)mxGetPr(b_mex);
  d->c = (scs_float *)mxGetPr(c_mex);
#endif

  A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  A->n = d->n;
  A->m = d->m;
  d->A = A;

  if (P_mex) {
    P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
    P->n = d->n;
    P->m = d->n;
    d->P = P;
  }

#ifdef DLONG
  A->p = (scs_int *)mxGetJc(A_mex);
  A->i = (scs_int *)mxGetIr(A_mex);
  if (P_mex) {
    P->p = (scs_int *)mxGetJc(P_mex);
    P->i = (scs_int *)mxGetIr(P_mex);
  }
#else
  A->p = cast_to_scs_int_arr(mxGetJc(A_mex), A->n + 1);
  if (!A->p) {
    free_mex(d, SCS_NULL, SCS_NULL);
    scs_printf("Memory allocation failed for A->p.\n");
    return -1;
  }
  A->i = cast_to_scs_int_arr(mxGetIr(A_mex), A->p[A->n]);
  if (!A->i) {
    free_mex(d, SCS_NULL, SCS_NULL);
    scs_printf("Memory allocation failed for A->i.\n");
    return -1;
  }
  if (P_mex) {
    P->p = cast_to_scs_int_arr(mxGetJc(P_mex), P->n + 1);
    if (!P->p) {
      free_mex(d, SCS_NULL, SCS_NULL);
      scs_printf("Memory allocation failed for P->p.\n");
      return -1;
    }
    P->i = cast_to_scs_int_arr(mxGetIr(P_mex), P->p[P->n]);
    if (!P->i) {
      free_mex(d, SCS_NULL, SCS_NULL);
      scs_printf("Memory allocation failed for P->i.\n");
      return -1;
    }
  }
#endif
#ifdef SFLOAT
  A->x = cast_to_scs_float_arr(mxGetPr(A_mex), A->p[A->n]);
  if (!A->x) {
    free_mex(d, SCS_NULL, SCS_NULL);
    scs_printf("Memory allocation failed for A->x.\n");
    return -1;
  }
  if (P_mex) {
    P->x = cast_to_scs_float_arr(mxGetPr(P_mex), P->p[P->n]);
    if (!P->x) {
      free_mex(d, SCS_NULL, SCS_NULL);
      scs_printf("Memory allocation failed for P->x.\n");
      return -1;
    }
  }
#else
  A->x = (scs_float *)mxGetPr(A_mex);
  if (P_mex) {
    P->x = (scs_float *)mxGetPr(P_mex);
  }
#endif

  *d_out = d;
  return 0;
}

/* Parse settings struct into ScsSettings.
 * Caller must free via free_mex(NULL, NULL, stgs). */
static void parse_settings(const mxArray *settings_mex, ScsSettings **stgs_out) {
  ScsSettings *stgs;
  mxArray *tmp;

  stgs = (ScsSettings *)scs_malloc(sizeof(ScsSettings));
  scs_set_default_settings(stgs);

#define GET_SETTING_FLOAT(field)                                               \
  tmp = mxGetField(settings_mex, 0, #field);                                   \
  if (tmp != SCS_NULL && !mxIsEmpty(tmp))                                      \
  stgs->field = (scs_float)*mxGetPr(tmp)

#define GET_SETTING_INT(field)                                                 \
  tmp = mxGetField(settings_mex, 0, #field);                                   \
  if (tmp != SCS_NULL && !mxIsEmpty(tmp))                                      \
  stgs->field = (scs_int)*mxGetPr(tmp)

  GET_SETTING_FLOAT(alpha);
  GET_SETTING_FLOAT(rho_x);
  GET_SETTING_INT(max_iters);
  GET_SETTING_FLOAT(scale);
  GET_SETTING_FLOAT(eps_abs);
  GET_SETTING_FLOAT(eps_rel);
  GET_SETTING_FLOAT(eps_infeas);
  GET_SETTING_INT(verbose);
  GET_SETTING_INT(normalize);
  GET_SETTING_INT(acceleration_lookback);
  GET_SETTING_INT(acceleration_interval);
  GET_SETTING_INT(adaptive_scale);
  GET_SETTING_FLOAT(time_limit_secs);

#undef GET_SETTING_FLOAT
#undef GET_SETTING_INT

  tmp = mxGetField(settings_mex, 0, "write_data_filename");
  if (tmp != SCS_NULL && !mxIsEmpty(tmp)) {
    /* need to free this later */
    stgs->write_data_filename = mxArrayToString(tmp);
  }

  tmp = mxGetField(settings_mex, 0, "log_csv_filename");
  if (tmp != SCS_NULL && !mxIsEmpty(tmp)) {
    /* need to free this later */
    stgs->log_csv_filename = mxArrayToString(tmp);
  }

  *stgs_out = stgs;
}

/* Parse cone struct into ScsCone.
 * Caller must free via free_mex(NULL, k, NULL). */
static scs_int parse_cones(const mxArray *cone_mex, ScsCone **k_out) {
  ScsCone *k;
  scs_int i, ns, ncs, nbl, nbu, blen;
  mxArray *tmp;
  const double *tmp_mex;

  k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  if (!k) {
    return -1;
  }

#define GET_CONE_INT(field)                                                    \
  tmp = mxGetField(cone_mex, 0, #field);                                       \
  if (tmp && !mxIsEmpty(tmp))                                                  \
  k->field += (scs_int)*mxGetPr(tmp)

  /* TODO rm this */
  tmp = mxGetField(cone_mex, 0, "f");
  if (tmp && !mxIsEmpty(tmp)) {
    scs_printf("SCS deprecation warning: The 'f' field in the cone struct \n"
               "has been replaced by 'z' to better reflect the Zero cone. \n"
               "Please replace usage of 'f' with 'z'. If both 'f' and 'z' \n"
               "are set then we sum the two fields to get the final zero \n"
               "cone size.\n");
    k->z += (scs_int)*mxGetPr(tmp);
  }

  GET_CONE_INT(z);
  GET_CONE_INT(l);
  GET_CONE_INT(ep);
  GET_CONE_INT(ed);

#undef GET_CONE_INT

#define GET_CONE_ARR(field, size_field, type)                                  \
  tmp = mxGetField(cone_mex, 0, #field);                                       \
  if (tmp && !mxIsEmpty(tmp)) {                                                \
    if (!mxIsDouble(tmp)) {                                                    \
      scs_printf("Cone field `" #field "` must be a double array.\n");         \
      free_mex(SCS_NULL, k, SCS_NULL);                                         \
      return -1;                                                               \
    }                                                                          \
    tmp_mex = mxGetPr(tmp);                                                    \
    k->size_field = get_mex_length(tmp);                                       \
    k->field = (type *)scs_calloc(k->size_field, sizeof(type));                \
    if (!k->field) {                                                           \
      free_mex(SCS_NULL, k, SCS_NULL);                                         \
      return -1;                                                               \
    }                                                                          \
    for (i = 0; i < k->size_field; i++) {                                      \
      k->field[i] = (type)tmp_mex[i];                                          \
    }                                                                          \
  }

  GET_CONE_ARR(q, qsize, scs_int);
  GET_CONE_ARR(s, ssize, scs_int);
  GET_CONE_ARR(cs, cssize, scs_int);
  GET_CONE_ARR(p, psize, scs_float);

  {
    mxArray *kbl = mxGetField(cone_mex, 0, "bl");
    mxArray *kbu = mxGetField(cone_mex, 0, "bu");
    if (kbl && kbu && !mxIsEmpty(kbl) && !mxIsEmpty(kbu)) {
      if (!mxIsDouble(kbl) || !mxIsDouble(kbu)) {
        scs_printf("bl,bu cone entries must be double arrays.\n");
        free_mex(SCS_NULL, k, SCS_NULL);
        return -1;
      }
      nbl = (scs_int)mxGetNumberOfDimensions(kbl);
      nbu = (scs_int)mxGetNumberOfDimensions(kbu);
      if (nbl != nbu) {
        scs_printf("bl,bu cone entries not the same size.\n");
        free_mex(SCS_NULL, k, SCS_NULL);
        return -1;
      }
      const size_t *bl_dims = mxGetDimensions(kbl);
      const size_t *bu_dims = mxGetDimensions(kbu);
      for (i = 0; i < nbu; i++) {
        if (bl_dims[i] != bu_dims[i]) {
          scs_printf("bl,bu cone entries not the same size.\n");
          free_mex(SCS_NULL, k, SCS_NULL);
          return -1;
        }
      }
      blen = (scs_int)bu_dims[0];
      if (nbl > 1 && bu_dims[0] == 1) {
        blen = (scs_int)bu_dims[1];
      }
      k->bu = (scs_float *)scs_calloc(blen, sizeof(scs_float));
      k->bl = (scs_float *)scs_calloc(blen, sizeof(scs_float));
      if (!k->bu || !k->bl) {
        free_mex(SCS_NULL, k, SCS_NULL);
        return -1;
      }
      const double *bl_mex = mxGetPr(kbl);
      const double *bu_mex = mxGetPr(kbu);
      for (i = 0; i < blen; i++) {
        k->bl[i] = (scs_float)bl_mex[i];
        k->bu[i] = (scs_float)bu_mex[i];
      }
      k->bsize = blen + 1;
    }
  }

#ifdef USE_SPECTRAL_CONES
  GET_CONE_ARR(d, dsize, scs_int);
  GET_CONE_ARR(ell1, ell1_size, scs_int);

  {
    mxArray *knuc_m = mxGetField(cone_mex, 0, "nuc_m");
    mxArray *knuc_n = mxGetField(cone_mex, 0, "nuc_n");
    if (knuc_m && knuc_n && !mxIsEmpty(knuc_m) && !mxIsEmpty(knuc_n)) {
      if (!mxIsDouble(knuc_m) || !mxIsDouble(knuc_n)) {
        scs_printf("nuc_m, nuc_n cone entries must be double arrays.\n");
        free_mex(SCS_NULL, k, SCS_NULL);
        return -1;
      }
      k->nucsize = get_mex_length(knuc_m);
      k->nuc_m = (scs_int *)scs_calloc(k->nucsize, sizeof(scs_int));
      k->nuc_n = (scs_int *)scs_calloc(k->nucsize, sizeof(scs_int));
      if (!k->nuc_m || !k->nuc_n) {
        free_mex(SCS_NULL, k, SCS_NULL);
        return -1;
      }
      const double *nuc_m_mex = mxGetPr(knuc_m);
      const double *nuc_n_mex = mxGetPr(knuc_n);
      for (i = 0; i < k->nucsize; i++) {
        k->nuc_m[i] = (scs_int)nuc_m_mex[i];
        k->nuc_n[i] = (scs_int)nuc_n_mex[i];
      }
    }
  }

  {
    mxArray *ksl_n = mxGetField(cone_mex, 0, "sl_n");
    mxArray *ksl_k = mxGetField(cone_mex, 0, "sl_k");
    if (ksl_n && ksl_k && !mxIsEmpty(ksl_n) && !mxIsEmpty(ksl_k)) {
      if (!mxIsDouble(ksl_n) || !mxIsDouble(ksl_k)) {
        scs_printf("sl_n, sl_k cone entries must be double arrays.\n");
        free_mex(SCS_NULL, k, SCS_NULL);
        return -1;
      }
      k->sl_size = get_mex_length(ksl_n);
      k->sl_n = (scs_int *)scs_calloc(k->sl_size, sizeof(scs_int));
      k->sl_k = (scs_int *)scs_calloc(k->sl_size, sizeof(scs_int));
      if (!k->sl_n || !k->sl_k) {
        free_mex(SCS_NULL, k, SCS_NULL);
        return -1;
      }
      const double *sl_n_mex = mxGetPr(ksl_n);
      const double *sl_k_mex = mxGetPr(ksl_k);
      for (i = 0; i < k->sl_size; i++) {
        k->sl_n[i] = (scs_int)sl_n_mex[i];
        k->sl_k[i] = (scs_int)sl_k_mex[i];
      }
    }
  }
#endif

#undef GET_CONE_ARR

  *k_out = k;
  return 0;
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
  mxSetField(*plhs3, 0, "lin_sys_solver", mxCreateString(info->lin_sys_solver));

#define SET_INFO_FIELD(field)                                                  \
  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);                                    \
  *mxGetPr(tmp) = (double)info->field;                                         \
  mxSetField(*plhs3, 0, #field, tmp)

  SET_INFO_FIELD(iter);
  SET_INFO_FIELD(scale_updates);
  SET_INFO_FIELD(status_val);
  SET_INFO_FIELD(pobj);
  SET_INFO_FIELD(dobj);
  SET_INFO_FIELD(res_pri);
  SET_INFO_FIELD(res_dual);
  SET_INFO_FIELD(res_infeas);
  SET_INFO_FIELD(res_unbdd_a);
  SET_INFO_FIELD(res_unbdd_p);
  SET_INFO_FIELD(comp_slack);
  SET_INFO_FIELD(gap);
  SET_INFO_FIELD(scale);
  SET_INFO_FIELD(setup_time);
  SET_INFO_FIELD(solve_time);
  SET_INFO_FIELD(rejected_accel_steps);
  SET_INFO_FIELD(accepted_accel_steps);
  SET_INFO_FIELD(lin_sys_time);
  SET_INFO_FIELD(cone_time);
  SET_INFO_FIELD(accel_time);

#undef SET_INFO_FIELD
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
      if (!mxIsStruct(prhs[1]) || !mxIsStruct(prhs[2])) {
        scs_free(cmd);
        mexErrMsgTxt("Input arguments 2 and 3 must be structs.");
      }
      if (!mxIsEmpty(prhs[3]) && !mxIsStruct(prhs[3])) {
        scs_free(cmd);
        mexErrMsgTxt("Input argument 4 (settings) must be a struct.");
      }
      ws_cleanup(); /* free any existing workspace */
      if (parse_data(prhs[1], &d) < 0) {
        scs_free(cmd);
        mexErrMsgTxt("Error parsing data.");
      }
      if (parse_cones(prhs[2], &k) < 0) {
        free_mex(d, SCS_NULL, SCS_NULL);
        scs_free(cmd);
        mexErrMsgTxt("Error parsing cones.");
      }
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
#ifdef SFLOAT
        b_new = cast_to_scs_float_arr(mxGetPr(prhs[1]), ws_m);
        if (!b_new) {
          scs_free(cmd);
          mexErrMsgTxt("Memory allocation failed for b_new.");
        }
#else
        b_new = (scs_float *)mxGetPr(prhs[1]);
#endif
      }
      if (nrhs >= 3 && !mxIsEmpty(prhs[2])) {
#ifdef SFLOAT
        c_new = cast_to_scs_float_arr(mxGetPr(prhs[2]), ws_n);
        if (!c_new) {
          if (b_new) scs_free(b_new);
          scs_free(cmd);
          mexErrMsgTxt("Memory allocation failed for c_new.");
        }
#else
        c_new = (scs_float *)mxGetPr(prhs[2]);
#endif
      }
      scs_update(ws_work, b_new, c_new);
#ifdef SFLOAT
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
        if (!mxIsStruct(ws_data)) {
          scs_free(cmd);
          mexErrMsgTxt("Warm start argument must be a struct.");
        }
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

      if (!sol.x || !sol.y || !sol.s) {
        if (sol.x) scs_free(sol.x);
        if (sol.y) scs_free(sol.y);
        if (sol.s) scs_free(sol.s);
        scs_free(cmd);
        mexErrMsgTxt("Memory allocation failed for solution vectors.");
      }

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
    if (!mxIsStruct(prhs[0]) || !mxIsStruct(prhs[1])) {
      mexErrMsgTxt("First two input arguments must be structs.");
    }
    if (!mxIsEmpty(prhs[2]) && !mxIsStruct(prhs[2])) {
      mexErrMsgTxt("Third input argument (settings) must be a struct.");
    }

    if (parse_data(prhs[0], &d) < 0) {
      mexErrMsgTxt("Error parsing data.");
    }
    if (parse_cones(prhs[1], &k) < 0) {
      free_mex(d, SCS_NULL, SCS_NULL);
      mexErrMsgTxt("Error parsing cones.");
    }
    parse_settings(prhs[2], &stgs);

    /* warm-start */
    stgs->warm_start =
        parse_warm_start(mxGetField(prhs[0], 0, "x"), &(sol.x), d->n);
    stgs->warm_start |=
        parse_warm_start(mxGetField(prhs[0], 0, "y"), &(sol.y), d->m);
    stgs->warm_start |=
        parse_warm_start(mxGetField(prhs[0], 0, "s"), &(sol.s), d->m);

    if (!sol.x || !sol.y || !sol.s) {
      if (sol.x) scs_free(sol.x);
      if (sol.y) scs_free(sol.y);
      if (sol.s) scs_free(sol.s);
      free_mex(d, k, stgs);
      mexErrMsgTxt("Memory allocation failed for solution vectors.");
    }

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
#ifdef SFLOAT /* only free if copies, which is only when flags set */
    if (d->b) {
      scs_free(d->b);
    }
    if (d->c) {
      scs_free(d->c);
    }
#endif
    if (d->A) {
#ifndef DLONG /* only free if copies, which is only when flags set */
      if (d->A->p) {
        scs_free(d->A->p);
      }
      if (d->A->i) {
        scs_free(d->A->i);
      }
#endif
#ifdef SFLOAT /* only free if copies, which is only when flags set */
      if (d->A->x) {
        scs_free(d->A->x);
      }
#endif
      scs_free(d->A);
    }
    if (d->P) {
#ifndef DLONG /* only free if copies, which is only when flags set */
      if (d->P->p) {
        scs_free(d->P->p);
      }
      if (d->P->i) {
        scs_free(d->P->i);
      }
#endif
#ifdef SFLOAT /* only free if copies, which is only when flags set */
      if (d->P->x) {
        scs_free(d->P->x);
      }
#endif
      scs_free(d->P);
    }
    scs_free(d);
  }
}
