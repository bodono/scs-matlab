#include "glbopts.h"
#include "linalg.h"
#include "matrix.h"
#include "mex.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

void free_mex(ScsData *d, ScsCone *k, ScsSettings *stgs);

scs_int parse_warm_start(const mxArray *p_mex, scs_float **p, scs_int l) {
  *p = (scs_float *)scs_calloc(
      l, sizeof(scs_float)); /* this allocates memory used for ScsSolution */
  if (p_mex == SCS_NULL) {
    return 0;
  } else if (mxIsSparse(p_mex) || (scs_int)*mxGetDimensions(p_mex) != l) {
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* matlab usage: [x,y,s,info] = scs(data,cone,settings); */
  scs_int i, ns, status, nbl, nbu, blen, buflen;
  ScsData *d;
  ScsCone *k;
  ScsSettings *stgs;
  ScsSolution sol = {0};
  ScsInfo info;
  ScsMatrix *A;
  ScsMatrix *P = SCS_NULL;
  char *buf;

  const mxArray *data;
  const mxArray *A_mex;
  const mxArray *P_mex;
  const mxArray *b_mex;
  const mxArray *c_mex;

  const mxArray *kf;
  const mxArray *kz;
  const mxArray *kl;
  const mxArray *kbl;
  const mxArray *kbu;
  const mxArray *kq;
  const mxArray *ks;
  const mxArray *kep;
  const mxArray *ked;
  const mxArray *kp;
  const double *q_mex;
  const double *bl_mex;
  const double *bu_mex;
  const double *s_mex;
  const double *p_mex;
  const size_t *bl_dims;
  const size_t *bu_dims;
  const size_t *q_dims;
  const size_t *s_dims;
  const size_t *p_dims;

  const mxArray *cone;
  const mxArray *settings;

  const mwSize one[1] = {1};
  const int num_info_fields = 16;
  const char *info_fields[] = {
      "iter",       "status",     "pobj",          "dobj",
      "res_pri",    "res_dual",   "res_infeas",    "res_unbdd_a",
      "scale",      "status_val", "res_unbdd_p",   "gap",
      "setup_time", "solve_time", "scale_updates", "comp_slack"};
  mxArray *tmp;
#if VERBOSITY > 0
  scs_printf("SIZE OF mwSize = %i\n", (int)sizeof(mwSize));
  scs_printf("SIZE OF mwIndex = %i\n", (int)sizeof(mwIndex));
#endif

  if (nrhs != 3) {
    mexErrMsgTxt("Three arguments are required in this order: data struct, "
                 "cone struct, settings struct");
  }
  if (nlhs > 4) {
    mexErrMsgTxt("scs returns up to 4 output arguments only.");
  }
  d = (ScsData *)mxMalloc(sizeof(ScsData));
  stgs = (ScsSettings *)mxMalloc(sizeof(ScsSettings));
  k = (ScsCone *)mxMalloc(sizeof(ScsCone));
  data = prhs[0];

  A_mex = (mxArray *)mxGetField(data, 0, "A");
  if (A_mex == SCS_NULL) {
    scs_free(d);
    scs_free(k);
    mexErrMsgTxt("ScsData struct must contain a `A` entry.");
  }
  if (!mxIsSparse(A_mex)) {
    scs_free(d);
    scs_free(k);
    mexErrMsgTxt("Input matrix A must be in sparse format (pass in sparse(A))");
  }
  P_mex = (mxArray *)mxGetField(data, 0, "P"); /* can be SCS_NULL */
  if (P_mex && !mxIsSparse(P_mex)) {
    scs_free(d);
    scs_free(k);
    mexErrMsgTxt("Input matrix P must be in sparse format (pass in sparse(A))");
  }
  b_mex = (mxArray *)mxGetField(data, 0, "b");
  if (b_mex == SCS_NULL) {
    scs_free(d);
    scs_free(k);
    mexErrMsgTxt("ScsData struct must contain a `b` entry.");
  }
  if (mxIsSparse(b_mex)) {
    scs_free(d);
    scs_free(k);
    mexErrMsgTxt("Input vector b must be in dense format (pass in full(b))");
  }
  c_mex = (mxArray *)mxGetField(data, 0, "c");
  if (c_mex == SCS_NULL) {
    scs_free(d);
    scs_free(k);
    mexErrMsgTxt("ScsData struct must contain a `c` entry.");
  }
  if (mxIsSparse(c_mex)) {
    scs_free(d);
    scs_free(k);
    mexErrMsgTxt("Input vector c must be in dense format (pass in full(c))");
  }
  cone = prhs[1];
  settings = prhs[2];
  d->n = (scs_int) * (mxGetDimensions(c_mex));
  d->m = (scs_int) * (mxGetDimensions(b_mex));
#if SFLOAT > 0
  d->b = cast_to_scs_float_arr(mxGetPr(b_mex), d->m);
  d->c = cast_to_scs_float_arr(mxGetPr(c_mex), d->n);
#else
  d->b = (scs_float *)mxGetPr(b_mex);
  d->c = (scs_float *)mxGetPr(c_mex);
#endif
  scs_set_default_settings(stgs);

  /* settings */
  tmp = mxGetField(settings, 0, "alpha");
  if (tmp != SCS_NULL) {
    stgs->alpha = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "rho_x");
  if (tmp != SCS_NULL) {
    stgs->rho_x = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "max_iters");
  if (tmp != SCS_NULL) {
    stgs->max_iters = (scs_int)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "scale");
  if (tmp != SCS_NULL) {
    stgs->scale = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "eps_abs");
  if (tmp != SCS_NULL) {
    stgs->eps_abs = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "eps_rel");
  if (tmp != SCS_NULL) {
    stgs->eps_rel = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "eps_infeas");
  if (tmp != SCS_NULL) {
    stgs->eps_infeas = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "verbose");
  if (tmp != SCS_NULL) {
    stgs->verbose = (scs_int)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "normalize");
  if (tmp != SCS_NULL) {
    stgs->normalize = (scs_int)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "acceleration_lookback");
  if (tmp != SCS_NULL) {
    stgs->acceleration_lookback = (scs_int)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "acceleration_interval");
  if (tmp != SCS_NULL) {
    stgs->acceleration_interval = (scs_int)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "adaptive_scale");
  if (tmp != SCS_NULL) {
    stgs->adaptive_scale = (scs_int)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "time_limit_secs");
  if (tmp != SCS_NULL) {
    stgs->time_limit_secs = (scs_float)*mxGetPr(tmp);
  }

  tmp = mxGetField(settings, 0, "write_data_filename");
  if (tmp != SCS_NULL) {
    buflen = mxGetNumberOfElements(tmp) + 1;
    buf = (char *)scs_calloc(buflen, sizeof(char));
    /* Copy the string data from tmp and place it into buf. */
    if (mxGetString(tmp, buf, buflen) != 0) {
      mexErrMsgIdAndTxt("MATLAB:explore:invalidStringArray",
                        "Could not convert string data.");
      mxFree(buf);
    } else {
      stgs->write_data_filename = buf;
    }
  }

  tmp = mxGetField(settings, 0, "log_csv_filename");
  if (tmp != SCS_NULL) {
    buflen = mxGetNumberOfElements(tmp) + 1;
    buf = (char *)scs_calloc(buflen, sizeof(char));
    /* Copy the string data from tmp and place it into buf. */
    if (mxGetString(tmp, buf, buflen) != 0) {
      mexErrMsgIdAndTxt("MATLAB:explore:invalidStringArray",
                        "Could not convert string data.");
      mxFree(buf);
    } else {
      stgs->log_csv_filename = buf;
    }
  }

  /* cones */

  /* TODO rm this */
  kf = mxGetField(cone, 0, "f");
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

  kz = mxGetField(cone, 0, "z");
  if (kz && !mxIsEmpty(kz)) {
    k->z += (scs_int)*mxGetPr(kz); /* TODO rm this */
  }

  kl = mxGetField(cone, 0, "l");
  if (kl && !mxIsEmpty(kl)) {
    k->l = (scs_int)*mxGetPr(kl);
  } else {
    k->l = 0;
  }

  kep = mxGetField(cone, 0, "ep");
  if (kep && !mxIsEmpty(kep)) {
    k->ep = (scs_int)*mxGetPr(kep);
  } else {
    k->ep = 0;
  }

  ked = mxGetField(cone, 0, "ed");
  if (ked && !mxIsEmpty(ked)) {
    k->ed = (scs_int)*mxGetPr(ked);
  } else {
    k->ed = 0;
  }

  kbu = mxGetField(cone, 0, "bu");
  kbl = mxGetField(cone, 0, "bl");
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
    k->bu = (scs_float *)mxMalloc(sizeof(scs_int) * blen);
    k->bl = (scs_float *)mxMalloc(sizeof(scs_int) * blen);

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

  kq = mxGetField(cone, 0, "q");
  if (kq && !mxIsEmpty(kq)) {
    q_mex = mxGetPr(kq);
    ns = (scs_int)mxGetNumberOfDimensions(kq);
    q_dims = mxGetDimensions(kq);
    k->qsize = (scs_int)q_dims[0];
    if (ns > 1 && q_dims[0] == 1) {
      k->qsize = (scs_int)q_dims[1];
    }
    k->q = (scs_int *)mxMalloc(sizeof(scs_int) * k->qsize);
    for (i = 0; i < k->qsize; i++) {
      k->q[i] = (scs_int)q_mex[i];
    }
  } else {
    k->qsize = 0;
    k->q = SCS_NULL;
  }

  ks = mxGetField(cone, 0, "s");
  if (ks && !mxIsEmpty(ks)) {
    s_mex = mxGetPr(ks);
    ns = (scs_int)mxGetNumberOfDimensions(ks);
    s_dims = mxGetDimensions(ks);
    k->ssize = (scs_int)s_dims[0];
    if (ns > 1 && s_dims[0] == 1) {
      k->ssize = (scs_int)s_dims[1];
    }
    k->s = (scs_int *)mxMalloc(sizeof(scs_int) * k->ssize);
    for (i = 0; i < k->ssize; i++) {
      k->s[i] = (scs_int)s_mex[i];
    }
  } else {
    k->ssize = 0;
    k->s = SCS_NULL;
  }

  kp = mxGetField(cone, 0, "p");
  if (kp && !mxIsEmpty(kp)) {
    p_mex = mxGetPr(kp);
    ns = (scs_int)mxGetNumberOfDimensions(kp);
    p_dims = mxGetDimensions(kp);
    k->psize = (scs_int)p_dims[0];
    if (ns > 1 && p_dims[0] == 1) {
      k->psize = (scs_int)p_dims[1];
    }
    k->p = (scs_float *)mxMalloc(sizeof(scs_float) * k->psize);
    for (i = 0; i < k->psize; i++) {
      k->p[i] = (scs_float)p_mex[i];
    }
  } else {
    k->psize = 0;
    k->p = SCS_NULL;
  }

  A = (ScsMatrix *)scs_malloc(sizeof(ScsMatrix));
  A->n = d->n;
  A->m = d->m;

  if (P_mex) {
    P = (ScsMatrix *)scs_malloc(sizeof(ScsMatrix));
    P->n = d->n;
    P->m = d->n;
  }

  /* TODO:
   * these return (mwIndex *), equivalent to (size_t *)
   * casting as (scs_int *), when scs_int = long seems to work
   * although maybe not on all machines:
   *
   * If scs_int is not long, then we explictly cast the entire
   * array to get the correct width
   */
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
  /* warm-start inputs, allocates sol->x, ->y, ->s even if warm start not used
   */
  stgs->warm_start =
      parse_warm_start((mxArray *)mxGetField(data, 0, "x"), &(sol.x), d->n);
  stgs->warm_start |=
      parse_warm_start((mxArray *)mxGetField(data, 0, "y"), &(sol.y), d->m);
  stgs->warm_start |=
      parse_warm_start((mxArray *)mxGetField(data, 0, "s"), &(sol.s), d->m);

  status = scs(d, k, stgs, &sol, &info);

  set_output_field(&plhs[0], sol.x, d->n);
  set_output_field(&plhs[1], sol.y, d->m);
  set_output_field(&plhs[2], sol.s, d->m);

  plhs[3] = mxCreateStructArray(1, one, num_info_fields, info_fields);

  /* if you add/remove fields here update the info_fields above */
  mxSetField(plhs[3], 0, "status", mxCreateString(info.status));

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "iter", tmp);
  *mxGetPr(tmp) = (scs_float)info.iter;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "scale_updates", tmp);
  *mxGetPr(tmp) = (scs_float)info.scale_updates;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "status_val", tmp);
  *mxGetPr(tmp) = (scs_float)info.status_val;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "pobj", tmp);
  *mxGetPr(tmp) = info.pobj;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "dobj", tmp);
  *mxGetPr(tmp) = info.dobj;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "res_pri", tmp);
  *mxGetPr(tmp) = info.res_pri;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "res_dual", tmp);
  *mxGetPr(tmp) = info.res_dual;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "res_infeas", tmp);
  *mxGetPr(tmp) = info.res_infeas;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "res_unbdd_a", tmp);
  *mxGetPr(tmp) = info.res_unbdd_a;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "res_unbdd_p", tmp);
  *mxGetPr(tmp) = info.res_unbdd_p;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "comp_slack", tmp);
  *mxGetPr(tmp) = info.comp_slack;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "gap", tmp);
  *mxGetPr(tmp) = info.gap;

  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "scale", tmp);
  *mxGetPr(tmp) = info.scale;

  /*info.time is millisecs - return value in secs */
  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "setup_time", tmp);
  *mxGetPr(tmp) = info.setup_time;

  /*info.time is millisecs - return value in secs */
  tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "solve_time", tmp);
  *mxGetPr(tmp) = info.solve_time;

  free_mex(d, k, stgs);
  return;
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
    if (k->p) {
      scs_free(k->p);
    }
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
