#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

SEXP cvr_c_simplify(SEXP x, SEXP tolerance, SEXP preserve_boundary);
SEXP cvr_c_is_valid(SEXP x, SEXP gap_width);
SEXP cvr_c_union(SEXP x);
SEXP cvr_c_n_coord(SEXP x);
SEXP cvr_c_unary_union(SEXP x);
SEXP cvr_c_simplify_each(SEXP x, SEXP tolerance);
SEXP cvr_c_area(SEXP x);
SEXP cvr_c_geos_version(void);

static const R_CallMethodDef CallEntries[] = {
  {"cvr_c_simplify",      (DL_FUNC) &cvr_c_simplify,      3},
  {"cvr_c_is_valid",      (DL_FUNC) &cvr_c_is_valid,      2},
  {"cvr_c_union",         (DL_FUNC) &cvr_c_union,         1},
  {"cvr_c_n_coord",       (DL_FUNC) &cvr_c_n_coord,       1},
  {"cvr_c_unary_union",   (DL_FUNC) &cvr_c_unary_union,   1},
  {"cvr_c_simplify_each", (DL_FUNC) &cvr_c_simplify_each, 2},
  {"cvr_c_area",          (DL_FUNC) &cvr_c_area,          1},
  {"cvr_c_geos_version",  (DL_FUNC) &cvr_c_geos_version,  0},
  {NULL, NULL, 0}
};

void attribute_visible R_init_cvr(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
