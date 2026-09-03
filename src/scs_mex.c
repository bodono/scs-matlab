#include "glbopts.h"
#include "linalg.h"
#include "matrix.h"
#include "mex.h"
#include <math.h>
#include <stdint.h>
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

#include <string.h>

void free_mex(ScsData *d, ScsCone *k, ScsSettings *stgs);

/* ======================== Workspace state ======================== */
typedef struct ScsMexWorkspace {
  uint64_t id;
  ScsWork *work;
  scs_int n;
  scs_int m;
  struct ScsMexWorkspace *next;
} ScsMexWorkspace;

static ScsMexWorkspace *ws_head = SCS_NULL;
static uint64_t ws_next_id = 1;

static void ws_cleanup_all(void) {
  ScsMexWorkspace *entry = ws_head;
  while (entry) {
    ScsMexWorkspace *next = entry->next;
    scs_finish(entry->work);
    scs_free(entry);
    mexUnlock();
    entry = next;
  }
  ws_head = SCS_NULL;
}

static ScsMexWorkspace *ws_find(uint64_t id) {
  ScsMexWorkspace *entry;
  for (entry = ws_head; entry; entry = entry->next) {
    if (entry->id == id) return entry;
  }
  return SCS_NULL;
}

static ScsMexWorkspace *ws_register(ScsWork *work, scs_int n, scs_int m) {
  ScsMexWorkspace *entry = (ScsMexWorkspace *)scs_malloc(sizeof(ScsMexWorkspace));
  if (!entry) return SCS_NULL;
  /* Handles count up within one load of the MEX. A live workspace holds
   * mexLock, so the MEX cannot be cleared underneath it; a handle kept
   * across an explicit clear of the MEX is invalid, like any other stale
   * handle. */
  entry->id = ws_next_id++;
  entry->work = work;
  entry->n = n;
  entry->m = m;
  entry->next = ws_head;
  ws_head = entry;
  /* Keep static registry state alive until the matching workspace is freed. */
  mexLock();
  return entry;
}

static void ws_remove(ScsMexWorkspace *target) {
  ScsMexWorkspace **link = &ws_head;
  while (*link && *link != target) link = &((*link)->next);
  if (*link) {
    *link = target->next;
    scs_finish(target->work);
    scs_free(target);
    mexUnlock();
  }
}

/* The workspace a handle argument names; raises (after releasing the
 * command string) if it is not a live handle. */
static ScsMexWorkspace *ws_lookup(const mxArray *value, char *cmd) {
  ScsMexWorkspace *entry = SCS_NULL;
  if (value && mxIsUint64(value) && mxGetNumberOfElements(value) == 1) {
    entry = ws_find(*((const uint64_t *)mxGetData(value)));
  }
  if (!entry) {
    scs_free(cmd);
    mexErrMsgIdAndTxt("scs:invalidWorkspace",
                      "Invalid or expired SCS workspace handle.");
  }
  return entry;
}

static mxArray *workspace_id_to_mex(uint64_t id) {
  mxArray *value = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *((uint64_t *)mxGetData(value)) = id;
  return value;
}

/* ======================== Helper functions ======================== */

scs_int parse_warm_start(const mxArray *p_mex, scs_float **p, scs_int l) {
  if (p_mex == SCS_NULL) {
    *p = (scs_float *)scs_calloc(l, sizeof(scs_float));
    return 0;
  } else if (mxIsSparse(p_mex) || mxIsComplex(p_mex) || !mxIsDouble(p_mex) ||
             mxGetNumberOfElements(p_mex) != (size_t)l) {
    /* class checks matter, not just length: mxGetPr returns NULL for
     * non-double inputs (logical/single/char), which the native copy
     * below would dereference */
    *p = (scs_float *)scs_calloc(l, sizeof(scs_float));
    scs_printf("Error parsing warm start input (vectors must be dense real "
               "double arrays of the correct size), running without full "
               "warm-start\n");
    return 0;
  } else {
#ifdef SFLOAT
    *p = cast_to_scs_float_arr(mxGetPr(p_mex), l);
#else
    *p = (scs_float *)scs_calloc(l, sizeof(scs_float));
    if (*p) {
      memcpy(*p, mxGetPr(p_mex), l * sizeof(scs_float));
    }
#endif
    return (*p != SCS_NULL);
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

/* Reject sparse/complex/non-double arrays and non-vector shapes. MATLAB
 * sparse doubles pass mxIsDouble but expose only stored entries through
 * mxGetPr while reporting full numel (out-of-bounds reads), and mxGetPr
 * returns NULL for non-double types (NULL dereference). Entries may be
 * infinite: box cone bounds bl/bu legitimately use +/-Inf. */
static scs_int validate_dense_real_vector(const mxArray *arr,
                                          const char *name) {
  if (mxIsSparse(arr) || mxIsComplex(arr) || !mxIsDouble(arr)) {
    scs_printf("Field `%s` must be a dense real double array.\n", name);
    return -1;
  }
  if (mxGetNumberOfDimensions(arr) > 2 ||
      (mxGetM(arr) > 1 && mxGetN(arr) > 1)) {
    scs_printf("Field `%s` must be a vector.\n", name);
    return -1;
  }
  return 0;
}

static scs_int cast_scs_int_checked(double v, const char *name,
                                    scs_int *out);

/* Cone arrays parsed by GET_CONE_ARR additionally require finite entries;
 * integer dimension lists (`integer` set) must also be integer-valued and
 * in scs_int range, so 2.5 is rejected instead of truncated. Power-cone
 * exponents are floats. */
static scs_int validate_cone_array(const mxArray *arr, const char *name,
                                   scs_int integer) {
  size_t i, n;
  const double *p;
  scs_int unused;
  if (validate_dense_real_vector(arr, name) < 0) {
    return -1;
  }
  n = mxGetNumberOfElements(arr);
  p = mxGetPr(arr);
  for (i = 0; i < n; i++) {
    if (integer) {
      if (cast_scs_int_checked(p[i], name, &unused) < 0) {
        return -1;
      }
    } else if (!isfinite(p[i])) {
      scs_printf("Cone field `%s` entries must be finite (entry %ld "
                 "is %g).\n",
                 name, (long)i, p[i]);
      return -1;
    }
  }
  return 0;
}

/* Read a scalar struct field as a double. Accepts real double scalars and
 * logical scalars (verbose=true is idiomatic MATLAB); rejects everything
 * else -- mxGetPr returns NULL for non-double arrays, so the old
 * unconditional *mxGetPr(tmp) dereferenced NULL for e.g. logical or char
 * inputs, and a sparse scalar with no stored entry reads out of bounds. */
static scs_int get_scalar_field(const mxArray *arr, const char *name,
                                double *out) {
  if (mxIsLogicalScalar(arr)) {
    *out = mxIsLogicalScalarTrue(arr) ? 1.0 : 0.0;
    return 0;
  }
  if (mxIsSparse(arr) || mxIsComplex(arr) || !mxIsDouble(arr) ||
      mxGetNumberOfElements(arr) != 1) {
    scs_printf("Field `%s` must be a real double (or logical) scalar.\n",
               name);
    return -1;
  }
  *out = *mxGetPr(arr);
  return 0;
}

/* Checked double -> scs_int conversion: the value must be finite,
 * integral, and strictly inside the active scs_int range (int32 for
 * forced-int32 builds such as cuDSS, int64 otherwise). 2^(w-1) is
 * exactly representable as a double for both widths, so v < lim also
 * rejects +/-Inf and NaN. */
static scs_int cast_scs_int_checked(double v, const char *name,
                                    scs_int *out) {
  double lim = ldexp(1.0, (int)(8 * sizeof(scs_int)) - 1);
  if (!(v > -lim && v < lim) || v != floor(v)) {
    scs_printf("Field `%s` must be a finite integer representable as "
               "scs_int (got %g).\n",
               name, v);
    return -1;
  }
  *out = (scs_int)v;
  return 0;
}

static scs_int get_cone_count_field(const mxArray *arr, const char *name,
                                    scs_int *out) {
  double v;
  if (get_scalar_field(arr, name, &v) < 0) {
    return -1;
  }
  if (!(v >= 0)) {
    scs_printf("Cone field `%s` must be nonnegative (got %g).\n", name, v);
    return -1;
  }
  return cast_scs_int_checked(v, name, out);
}

#ifdef USE_SPECTRAL_CONES
/* Two cone arrays that describe one cone type together (nuc_m/nuc_n,
 * sl_n/sl_k): both present or both absent, equal length, finite entries.
 * The copy loop reads both with the same length, so this is what keeps it
 * in bounds. */
static scs_int parse_paired_int_arrays(const mxArray *cone_mex,
                                       const char *name_a, const char *name_b,
                                       scs_int **out_a, scs_int **out_b,
                                       scs_int *size) {
  mxArray *a = mxGetField(cone_mex, 0, name_a);
  mxArray *b = mxGetField(cone_mex, 0, name_b);
  scs_int has_a = a && !mxIsEmpty(a), has_b = b && !mxIsEmpty(b), i;
  const double *pa, *pb;
  if (has_a != has_b) {
    scs_printf("%s and %s must both be provided and nonempty.\n", name_a,
               name_b);
    return -1;
  }
  if (!has_a) {
    return 0;
  }
  if (validate_cone_array(a, name_a, 1) < 0 ||
      validate_cone_array(b, name_b, 1) < 0) {
    return -1;
  }
  if (mxGetNumberOfElements(a) != mxGetNumberOfElements(b)) {
    scs_printf("%s and %s must have the same number of entries.\n", name_a,
               name_b);
    return -1;
  }
  *size = (scs_int)mxGetNumberOfElements(a);
  *out_a = (scs_int *)scs_calloc(*size, sizeof(scs_int));
  *out_b = (scs_int *)scs_calloc(*size, sizeof(scs_int));
  if (!*out_a || !*out_b) {
    return -1;
  }
  pa = mxGetPr(a);
  pb = mxGetPr(b);
  for (i = 0; i < *size; i++) {
    (*out_a)[i] = (scs_int)pa[i];
    (*out_b)[i] = (scs_int)pb[i];
  }
  return 0;
}
#endif

/* Parse data struct (A, P, b, c) into ScsData.
 * Caller must free d, d->A, d->P (if non-NULL) via free_mex(d, NULL, NULL). */
static scs_int parse_data(const mxArray *data_mex, ScsData **d_out) {
  ScsData *d;
  ScsMatrix *A;
  ScsMatrix *P = SCS_NULL;
  const mxArray *A_mex, *P_mex, *b_mex, *c_mex;

  d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  if (!d) {
    return -1;
  }

  A_mex = (mxArray *)mxGetField(data_mex, 0, "A");
  if (A_mex == SCS_NULL) {
    scs_free(d);
    scs_printf("ScsData struct must contain a `A` entry.\n");
    return -1;
  }
  if (!mxIsSparse(A_mex) || mxIsComplex(A_mex) || !mxIsDouble(A_mex)) {
    scs_free(d);
    scs_printf("Input matrix A must be a sparse real double matrix "
               "(pass in sparse(A))\n");
    return -1;
  }
  P_mex = (mxArray *)mxGetField(data_mex, 0, "P"); /* can be SCS_NULL */
  if (P_mex && (!mxIsSparse(P_mex) || mxIsComplex(P_mex) ||
                !mxIsDouble(P_mex))) {
    scs_free(d);
    scs_printf("Input matrix P must be a sparse real double matrix "
               "(pass in sparse(P))\n");
    return -1;
  }
  b_mex = (mxArray *)mxGetField(data_mex, 0, "b");
  if (b_mex == SCS_NULL) {
    scs_free(d);
    scs_printf("ScsData struct must contain a `b` entry.\n");
    return -1;
  }
  if (validate_dense_real_vector(b_mex, "b") < 0) {
    scs_free(d);
    return -1;
  }
  c_mex = (mxArray *)mxGetField(data_mex, 0, "c");
  if (c_mex == SCS_NULL) {
    scs_free(d);
    scs_printf("ScsData struct must contain a `c` entry.\n");
    return -1;
  }
  if (validate_dense_real_vector(c_mex, "c") < 0) {
    scs_free(d);
    return -1;
  }

  /* dimensions from numel (a row vector's leading dimension is 1 and would
   * silently drop elements), and A's shape must agree, or the CSC copy
   * below reads past mxGetJc/mxGetIr */
  d->n = (scs_int)mxGetNumberOfElements(c_mex);
  d->m = (scs_int)mxGetNumberOfElements(b_mex);
  if (mxGetM(A_mex) != (size_t)d->m || mxGetN(A_mex) != (size_t)d->n) {
    scs_free(d);
    scs_printf("A must be %ld x %ld to match numel(b) x numel(c), got "
               "%ld x %ld.\n",
               (long)d->m, (long)d->n, (long)mxGetM(A_mex),
               (long)mxGetN(A_mex));
    return -1;
  }
  if (P_mex && (mxGetM(P_mex) != (size_t)d->n ||
                mxGetN(P_mex) != (size_t)d->n)) {
    scs_free(d);
    scs_printf("P must be %ld x %ld (square, matching numel(c)), got "
               "%ld x %ld.\n",
               (long)d->n, (long)d->n, (long)mxGetM(P_mex),
               (long)mxGetN(P_mex));
    return -1;
  }
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
  if (!A) {
    free_mex(d, SCS_NULL, SCS_NULL);
    return -1;
  }
  A->n = d->n;
  A->m = d->m;
  d->A = A;

  if (P_mex) {
    P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
    if (!P) {
      free_mex(d, SCS_NULL, SCS_NULL);
      return -1;
    }
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
    scs_printf("Could not convert A->p (allocation or index range).\n");
    return -1;
  }
  A->i = cast_to_scs_int_arr(mxGetIr(A_mex), A->p[A->n]);
  if (!A->i) {
    free_mex(d, SCS_NULL, SCS_NULL);
    scs_printf("Could not convert A->i (allocation or index range).\n");
    return -1;
  }
  if (P_mex) {
    P->p = cast_to_scs_int_arr(mxGetJc(P_mex), P->n + 1);
    if (!P->p) {
      free_mex(d, SCS_NULL, SCS_NULL);
      scs_printf("Could not convert P->p (allocation or index range).\n");
      return -1;
    }
    P->i = cast_to_scs_int_arr(mxGetIr(P_mex), P->p[P->n]);
    if (!P->i) {
      free_mex(d, SCS_NULL, SCS_NULL);
      scs_printf("Could not convert P->i (allocation or index range).\n");
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
static scs_int parse_settings(const mxArray *settings_mex,
                              ScsSettings **stgs_out) {
  ScsSettings *stgs;
  mxArray *tmp;

  stgs = (ScsSettings *)scs_malloc(sizeof(ScsSettings));
  if (!stgs) {
    return -1;
  }
  scs_set_default_settings(stgs);

/* Non-finite/fractional domain checks belong to the core's validate();
 * the MEX layer only guarantees a type-safe read. */
#define GET_SETTING_FLOAT(field)                                               \
  tmp = mxGetField(settings_mex, 0, #field);                                   \
  if (tmp != SCS_NULL && !mxIsEmpty(tmp)) {                                    \
    double v_;                                                                 \
    if (get_scalar_field(tmp, #field, &v_) < 0) {                              \
      scs_free(stgs);                                                          \
      return -1;                                                               \
    }                                                                          \
    stgs->field = (scs_float)v_;                                               \
  }

#define GET_SETTING_INT(field)                                                 \
  tmp = mxGetField(settings_mex, 0, #field);                                   \
  if (tmp != SCS_NULL && !mxIsEmpty(tmp)) {                                    \
    double v_;                                                                 \
    if (get_scalar_field(tmp, #field, &v_) < 0) {                              \
      scs_free(stgs);                                                          \
      return -1;                                                               \
    }                                                                          \
    if (cast_scs_int_checked(v_, #field, &stgs->field) < 0) {                  \
      scs_free(stgs);                                                          \
      return -1;                                                               \
    }                                                                          \
  }

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
  GET_SETTING_INT(acceleration_type_1);
  GET_SETTING_FLOAT(acceleration_regularization);
  GET_SETTING_FLOAT(acceleration_relaxation);
  GET_SETTING_INT(adaptive_scale);
  GET_SETTING_INT(adaptive_diag_scale);
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
  return 0;
}

/* Parse cone struct into ScsCone.
 * Caller must free via free_mex(NULL, k, NULL). */
static scs_int parse_cones(const mxArray *cone_mex, ScsCone **k_out) {
  ScsCone *k;
  scs_int i, blen;
  mxArray *tmp;
  const double *tmp_mex;

  k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  if (!k) {
    return -1;
  }

#define GET_CONE_INT(field)                                                    \
  tmp = mxGetField(cone_mex, 0, #field);                                       \
  if (tmp && !mxIsEmpty(tmp)) {                                                \
    scs_int v_;                                                                \
    if (get_cone_count_field(tmp, #field, &v_) < 0) {                          \
      free_mex(SCS_NULL, k, SCS_NULL);                                         \
      return -1;                                                               \
    }                                                                          \
    k->field += v_;                                                            \
  }

  /* TODO rm this */
  tmp = mxGetField(cone_mex, 0, "f");
  if (tmp && !mxIsEmpty(tmp)) {
    scs_printf("SCS deprecation warning: The 'f' field in the cone struct \n"
               "has been replaced by 'z' to better reflect the Zero cone. \n"
               "Please replace usage of 'f' with 'z'. If both 'f' and 'z' \n"
               "are set then we sum the two fields to get the final zero \n"
               "cone size.\n");
    {
      scs_int v_;
      if (get_cone_count_field(tmp, "f", &v_) < 0) {
        free_mex(SCS_NULL, k, SCS_NULL);
        return -1;
      }
      k->z += v_;
    }
  }

  GET_CONE_INT(z);
  GET_CONE_INT(l);
  GET_CONE_INT(ep);
  GET_CONE_INT(ed);

#undef GET_CONE_INT

#define GET_CONE_ARR(field, size_field, type, integer)                         \
  tmp = mxGetField(cone_mex, 0, #field);                                       \
  if (tmp && !mxIsEmpty(tmp)) {                                                \
    if (validate_cone_array(tmp, #field, integer) < 0) {                       \
      free_mex(SCS_NULL, k, SCS_NULL);                                         \
      return -1;                                                               \
    }                                                                          \
    tmp_mex = mxGetPr(tmp);                                                    \
    k->size_field = (scs_int)mxGetNumberOfElements(tmp);                       \
    k->field = (type *)scs_calloc(k->size_field, sizeof(type));                \
    if (!k->field) {                                                           \
      free_mex(SCS_NULL, k, SCS_NULL);                                         \
      return -1;                                                               \
    }                                                                          \
    for (i = 0; i < k->size_field; i++) {                                      \
      k->field[i] = (type)tmp_mex[i];                                          \
    }                                                                          \
  }

  GET_CONE_ARR(q, qsize, scs_int, 1);
  GET_CONE_ARR(s, ssize, scs_int, 1);
  GET_CONE_ARR(cs, cssize, scs_int, 1);
  GET_CONE_ARR(p, psize, scs_float, 0);

  {
    mxArray *kbl = mxGetField(cone_mex, 0, "bl");
    mxArray *kbu = mxGetField(cone_mex, 0, "bu");
    scs_int has_bl = kbl && !mxIsEmpty(kbl);
    scs_int has_bu = kbu && !mxIsEmpty(kbu);
    if (has_bl != has_bu) {
      /* a one-sided box would otherwise be silently dropped */
      scs_printf("bl and bu must both be provided and nonempty.\n");
      free_mex(SCS_NULL, k, SCS_NULL);
      return -1;
    }
    if (has_bl) {
      if (validate_dense_real_vector(kbl, "bl") < 0 ||
          validate_dense_real_vector(kbu, "bu") < 0) {
        free_mex(SCS_NULL, k, SCS_NULL);
        return -1;
      }
      if (mxGetNumberOfElements(kbl) != mxGetNumberOfElements(kbu)) {
        scs_printf("bl,bu cone entries not the same size.\n");
        free_mex(SCS_NULL, k, SCS_NULL);
        return -1;
      }
      blen = (scs_int)mxGetNumberOfElements(kbu);
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
  GET_CONE_ARR(d, dsize, scs_int, 1);
  GET_CONE_ARR(ell1, ell1_size, scs_int, 1);

  if (parse_paired_int_arrays(cone_mex, "nuc_m", "nuc_n", &k->nuc_m,
                              &k->nuc_n, &k->nucsize) < 0 ||
      parse_paired_int_arrays(cone_mex, "sl_n", "sl_k", &k->sl_n, &k->sl_k,
                              &k->sl_size) < 0) {
    free_mex(SCS_NULL, k, SCS_NULL);
    return -1;
  }
#endif

#undef GET_CONE_ARR

  *k_out = k;
  return 0;
}

/* Write ScsInfo to a MATLAB struct */
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

/* MATLAB reserves at least one output slot for ans. Copy only requested
 * vectors and release every solver buffer, including unrequested ones. */
static void write_outputs(int nlhs, mxArray *plhs[], const ScsSolution *sol,
                          scs_int n, scs_int m, const ScsInfo *info) {
  scs_float *vectors[3] = {sol->x, sol->y, sol->s};
  scs_int lengths[3] = {n, m, m};
  int i, want = nlhs < 1 ? 1 : nlhs;
  for (i = 0; i < 3; i++) {
    if (i < want) {
      set_output_field(&plhs[i], vectors[i], lengths[i]);
    } else {
      scs_free(vectors[i]);
    }
  }
  if (want >= 4) {
    write_info(&plhs[3], info);
  }
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
      /* id = scs_xxx('init', data, cone, settings) */
      ScsData *d;
      ScsCone *k;
      ScsSettings *stgs;
      ScsWork *work;
      ScsMexWorkspace *entry;
      if (nrhs != 4) {
        scs_free(cmd);
        mexErrMsgIdAndTxt("scs:usage",
                          "Usage: id = scs_xxx('init', data, cone, settings)");
      }
      if (!mxIsStruct(prhs[1]) || !mxIsStruct(prhs[2])) {
        scs_free(cmd);
        mexErrMsgTxt("Input arguments 2 and 3 must be structs.");
      }
      if (!mxIsEmpty(prhs[3]) && !mxIsStruct(prhs[3])) {
        scs_free(cmd);
        mexErrMsgTxt("Input argument 4 (settings) must be a struct.");
      }
      if (parse_data(prhs[1], &d) < 0) {
        scs_free(cmd);
        mexErrMsgIdAndTxt("scs:invalidData", "Error parsing data.");
      }
      if (parse_cones(prhs[2], &k) < 0) {
        free_mex(d, SCS_NULL, SCS_NULL);
        scs_free(cmd);
        mexErrMsgIdAndTxt("scs:invalidCone", "Error parsing cones.");
      }
      if (parse_settings(prhs[3], &stgs) < 0) {
        free_mex(d, k, SCS_NULL);
        scs_free(cmd);
        mexErrMsgIdAndTxt("scs:invalidSettings", "Error parsing settings.");
      }

      work = scs_init(d, k, stgs);
      entry = work ? ws_register(work, d->n, d->m) : SCS_NULL;

      free_mex(d, k, stgs);

      if (!entry) {
        if (work) scs_finish(work);
        scs_free(cmd);
        mexErrMsgTxt("SCS init failed.");
      }
      mexAtExit(ws_cleanup_all); /* idempotent: only the last registration is kept */
      plhs[0] = workspace_id_to_mex(entry->id);
      scs_free(cmd);
      return;
    }

    if (strcmp(cmd, "update") == 0) {
      /* scs_xxx('update', id, b_new, c_new)
       * Either argument can be [] to leave unchanged. */
      scs_float *b_new = SCS_NULL;
      scs_float *c_new = SCS_NULL;
      ScsMexWorkspace *entry;
      if (nrhs < 2 || nrhs > 4) {
        scs_free(cmd);
        mexErrMsgIdAndTxt("scs:usage",
                          "Usage: scs_xxx('update', id[, b_new[, c_new]])");
      }
      entry = ws_lookup(prhs[1], cmd);
      if (nrhs >= 3 && !mxIsEmpty(prhs[2])) {
        /* the native update copies entry->m doubles from the pointer */
        if (validate_dense_real_vector(prhs[2], "b") < 0 ||
            mxGetNumberOfElements(prhs[2]) != (size_t)entry->m) {
          scs_free(cmd);
          mexErrMsgIdAndTxt(
              "scs:updateInvalidB",
              "b must be a dense, real double vector with %ld elements.",
              (long)entry->m);
        }
#ifdef SFLOAT
        b_new = cast_to_scs_float_arr(mxGetPr(prhs[2]), entry->m);
        if (!b_new) {
          scs_free(cmd);
          mexErrMsgTxt("Memory allocation failed for b_new.");
        }
#else
        b_new = (scs_float *)mxGetPr(prhs[2]);
#endif
      }
      if (nrhs >= 4 && !mxIsEmpty(prhs[3])) {
        if (validate_dense_real_vector(prhs[3], "c") < 0 ||
            mxGetNumberOfElements(prhs[3]) != (size_t)entry->n) {
#ifdef SFLOAT
          if (b_new) scs_free(b_new);
#endif
          scs_free(cmd);
          mexErrMsgIdAndTxt(
              "scs:updateInvalidC",
              "c must be a dense, real double vector with %ld elements.",
              (long)entry->n);
        }
#ifdef SFLOAT
        c_new = cast_to_scs_float_arr(mxGetPr(prhs[3]), entry->n);
        if (!c_new) {
          if (b_new) scs_free(b_new);
          scs_free(cmd);
          mexErrMsgTxt("Memory allocation failed for c_new.");
        }
#else
        c_new = (scs_float *)mxGetPr(prhs[3]);
#endif
      }
      scs_update(entry->work, b_new, c_new);
#ifdef SFLOAT
      if (b_new) scs_free(b_new);
      if (c_new) scs_free(c_new);
#endif
      scs_free(cmd);
      return;
    }

    if (strcmp(cmd, "solve") == 0) {
      /* [x,y,s,info] = scs_xxx('solve', id)
       * [x,y,s,info] = scs_xxx('solve', id, warm_start_struct) */
      ScsSolution sol = {0};
      ScsInfo info;
      scs_int warm_start = 0;
      ScsMexWorkspace *entry;
      if (nrhs < 2 || nrhs > 3) {
        scs_free(cmd);
        mexErrMsgIdAndTxt(
            "scs:usage",
            "Usage: [x,y,s,info] = scs_xxx('solve', id[, warm_start])");
      }
      entry = ws_lookup(prhs[1], cmd);
      if (nrhs >= 3 && !mxIsEmpty(prhs[2])) {
        const mxArray *ws_data = prhs[2];
        if (!mxIsStruct(ws_data)) {
          scs_free(cmd);
          mexErrMsgTxt("Warm start argument must be a struct.");
        }
        warm_start =
            parse_warm_start(mxGetField(ws_data, 0, "x"), &(sol.x), entry->n);
        warm_start |=
            parse_warm_start(mxGetField(ws_data, 0, "y"), &(sol.y), entry->m);
        warm_start |=
            parse_warm_start(mxGetField(ws_data, 0, "s"), &(sol.s), entry->m);
      }
      if (!sol.x)
        sol.x = (scs_float *)scs_calloc(entry->n, sizeof(scs_float));
      if (!sol.y)
        sol.y = (scs_float *)scs_calloc(entry->m, sizeof(scs_float));
      if (!sol.s)
        sol.s = (scs_float *)scs_calloc(entry->m, sizeof(scs_float));

      if (!sol.x || !sol.y || !sol.s) {
        if (sol.x) scs_free(sol.x);
        if (sol.y) scs_free(sol.y);
        if (sol.s) scs_free(sol.s);
        scs_free(cmd);
        mexErrMsgTxt("Memory allocation failed for solution vectors.");
      }

      scs_solve(entry->work, &sol, &info, warm_start);

      write_outputs(nlhs, plhs, &sol, entry->n, entry->m, &info);

      scs_free(cmd);
      return;
    }

    if (strcmp(cmd, "finish") == 0) {
      ScsMexWorkspace *entry;
      if (nrhs != 2) {
        scs_free(cmd);
        mexErrMsgIdAndTxt("scs:usage", "Usage: scs_xxx('finish', id)");
      }
      entry = ws_lookup(prhs[1], cmd);
      ws_remove(entry);
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
      mexErrMsgIdAndTxt("scs:invalidData", "Error parsing data.");
    }
    if (parse_cones(prhs[1], &k) < 0) {
      free_mex(d, SCS_NULL, SCS_NULL);
      mexErrMsgIdAndTxt("scs:invalidCone", "Error parsing cones.");
    }
    if (parse_settings(prhs[2], &stgs) < 0) {
      free_mex(d, k, SCS_NULL);
      mexErrMsgIdAndTxt("scs:invalidSettings", "Error parsing settings.");
    }

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

    write_outputs(nlhs, plhs, &sol, d->n, d->m, &info);

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
