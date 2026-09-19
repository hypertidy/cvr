#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <stdio.h>
#include <geos_c.h>

/* ------------------------------------------------------------------ *
 * Error capture
 *
 * GEOS reports failures through a message handler rather than a return
 * code, and we must not longjmp out of R while GEOS resources are held.
 * So messages are buffered here, every GEOS handle is released, and
 * Rf_error() is called last.
 * ------------------------------------------------------------------ */

#define CVR_ERRBUF 1024

typedef struct {
  char msg[CVR_ERRBUF];
  int seen;
} cvr_err;

static void cvr_handler(const char *message, void *userdata) {
  cvr_err *e = (cvr_err *) userdata;
  if (!e->seen) {
    snprintf(e->msg, CVR_ERRBUF, "%s", message == NULL ? "unknown GEOS error" : message);
    e->seen = 1;
  }
}

typedef struct {
  GEOSContextHandle_t ctx;
  GEOSWKBReader *reader;
  GEOSWKBWriter *writer;
  cvr_err err;
} cvr_ctx;

static void cvr_open(cvr_ctx *s) {
  s->err.seen = 0;
  s->err.msg[0] = '\0';
  s->ctx = GEOS_init_r();
  GEOSContext_setErrorMessageHandler_r(s->ctx, cvr_handler, &s->err);
  s->reader = GEOSWKBReader_create_r(s->ctx);
  s->writer = GEOSWKBWriter_create_r(s->ctx);
  /* Fixed little-endian 2D output so results are byte-stable across
     machines and always readable by wk. */
  GEOSWKBWriter_setByteOrder_r(s->ctx, s->writer, 1);
  GEOSWKBWriter_setOutputDimension_r(s->ctx, s->writer, 2);
  GEOSWKBWriter_setIncludeSRID_r(s->ctx, s->writer, 0);
}

/* Copy the buffered message out, then release every GEOS handle. */
static void cvr_close(cvr_ctx *s, char *msg_out) {
  snprintf(msg_out, CVR_ERRBUF, "%s", s->err.seen ? s->err.msg : "unknown GEOS error");
  if (s->writer != NULL) GEOSWKBWriter_destroy_r(s->ctx, s->writer);
  if (s->reader != NULL) GEOSWKBReader_destroy_r(s->ctx, s->reader);
  GEOS_finish_r(s->ctx);
}

/* ------------------------------------------------------------------ *
 * WKB list -> GEOS GEOMETRYCOLLECTION
 *
 * The coverage operations are set-level: they take one collection and
 * return one collection, positionally parallel to the input. That is
 * the whole reason this cannot be expressed as a row-wise operator.
 *
 * GEOSGeom_createCollection_r takes ownership of its members on
 * success, so on failure we free the members we built ourselves.
 * ------------------------------------------------------------------ */

static GEOSGeometry *cvr_collection(GEOSContextHandle_t ctx, GEOSWKBReader *reader,
                                    SEXP x, R_xlen_t n, R_xlen_t *bad) {
  GEOSGeometry **members = (GEOSGeometry **) R_alloc(n, sizeof(GEOSGeometry *));
  R_xlen_t i;
  GEOSGeometry *coll;

  for (i = 0; i < n; i++) members[i] = NULL;

  for (i = 0; i < n; i++) {
    SEXP item = VECTOR_ELT(x, i);
    if (item == R_NilValue || TYPEOF(item) != RAWSXP) {
      *bad = i + 1;
      goto fail;
    }
    members[i] = GEOSWKBReader_read_r(ctx, reader, RAW(item), (size_t) Rf_xlength(item));
    if (members[i] == NULL) {
      *bad = i + 1;
      goto fail;
    }
  }

  coll = GEOSGeom_createCollection_r(ctx, GEOS_GEOMETRYCOLLECTION,
                                     members, (unsigned int) n);
  if (coll == NULL) goto fail;
  return coll;

fail:
  for (i = 0; i < n; i++) {
    if (members[i] != NULL) GEOSGeom_destroy_r(ctx, members[i]);
  }
  return NULL;
}

/* GEOS collection -> R list of raw WKB, one element per member.
   Returns an unprotected SEXP; PROTECT it at the call site before any
   further allocation. Returns R_NilValue on failure. */

static SEXP cvr_unpack(GEOSContextHandle_t ctx, GEOSWKBWriter *writer,
                       const GEOSGeometry *coll) {
  int n = GEOSGetNumGeometries_r(ctx, coll);
  int i;
  SEXP out;

  if (n < 0) return R_NilValue;

  out = PROTECT(Rf_allocVector(VECSXP, n));
  for (i = 0; i < n; i++) {
    const GEOSGeometry *part = GEOSGetGeometryN_r(ctx, coll, i);
    size_t size = 0;
    unsigned char *buf;
    SEXP raw;

    if (part == NULL) { UNPROTECT(1); return R_NilValue; }

    buf = GEOSWKBWriter_write_r(ctx, writer, part, &size);
    if (buf == NULL) { UNPROTECT(1); return R_NilValue; }

    raw = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t) size));
    memcpy(RAW(raw), buf, size);
    GEOSFree_r(ctx, buf);
    SET_VECTOR_ELT(out, i, raw);
    UNPROTECT(1);
  }

  UNPROTECT(1);
  return out;
}

/* Raise the error the operation buffered. Never returns. */
static void cvr_stop(const char *what, R_xlen_t bad, const char *msg) {
  if (bad > 0) {
    Rf_error("could not read geometry %ld as WKB: %s", (long) bad, msg);
  }
  Rf_error("GEOS coverage %s failed: %s", what, msg);
}

/* ------------------------------------------------------------------ *
 * GEOSCoverageSimplifyVW
 * ------------------------------------------------------------------ */

SEXP cvr_c_simplify(SEXP x, SEXP tolerance, SEXP preserve_boundary) {
  R_xlen_t n = Rf_xlength(x);
  double tol = Rf_asReal(tolerance);
  int preserve = (Rf_asLogical(preserve_boundary) == TRUE) ? 1 : 0;
  cvr_ctx s;
  R_xlen_t bad = 0;
  const char *what = NULL;
  char msg[CVR_ERRBUF];
  SEXP out = R_NilValue;
  int nprot = 0;
  GEOSGeometry *coll;

  if (n == 0) return Rf_allocVector(VECSXP, 0);

  cvr_open(&s);

  coll = cvr_collection(s.ctx, s.reader, x, n, &bad);
  if (coll == NULL) {
    what = "read";
  } else {
    GEOSGeometry *res = GEOSCoverageSimplifyVW_r(s.ctx, coll, tol, preserve);
    GEOSGeom_destroy_r(s.ctx, coll);
    if (res == NULL) {
      what = "simplify";
    } else {
      out = cvr_unpack(s.ctx, s.writer, res);
      PROTECT(out); nprot++;
      GEOSGeom_destroy_r(s.ctx, res);
      if (out == R_NilValue) what = "write";
    }
  }

  cvr_close(&s, msg);

  if (what != NULL) {
    if (nprot) UNPROTECT(nprot);
    cvr_stop(what, bad, msg);
  }

  UNPROTECT(nprot);
  return out;
}

/* ------------------------------------------------------------------ *
 * GEOSCoverageIsValid
 *
 * Returns list(valid = logical(1), edges = <list of raw WKB>), where
 * edges is parallel to the input: a MULTILINESTRING of offending edges
 * per member, EMPTY where that member is a valid participant.
 * ------------------------------------------------------------------ */

SEXP cvr_c_is_valid(SEXP x, SEXP gap_width) {
  R_xlen_t n = Rf_xlength(x);
  double gap = Rf_asReal(gap_width);
  cvr_ctx s;
  R_xlen_t bad = 0;
  const char *what = NULL;
  char msg[CVR_ERRBUF];
  SEXP out = R_NilValue;
  int nprot = 0;
  int verdict = 2;
  GEOSGeometry *coll;

  if (n == 0) {
    SEXP empty = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(empty, 1, Rf_allocVector(VECSXP, 0));
    SET_VECTOR_ELT(empty, 0, Rf_ScalarLogical(TRUE));
    UNPROTECT(1);
    return empty;
  }

  cvr_open(&s);

  coll = cvr_collection(s.ctx, s.reader, x, n, &bad);
  if (coll == NULL) {
    what = "read";
  } else {
    GEOSGeometry *edges = NULL;
    verdict = GEOSCoverageIsValid_r(s.ctx, coll, gap, &edges);
    GEOSGeom_destroy_r(s.ctx, coll);

    if (verdict == 2) {
      if (edges != NULL) GEOSGeom_destroy_r(s.ctx, edges);
      what = "validate";
    } else {
      SEXP edge_list;
      out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
      edge_list = (edges == NULL) ? Rf_allocVector(VECSXP, 0)
                                  : cvr_unpack(s.ctx, s.writer, edges);
      if (edges != NULL) GEOSGeom_destroy_r(s.ctx, edges);

      if (edge_list == R_NilValue) {
        what = "write";
      } else {
        /* Fill the allocating slot first: edge_list is not protected. */
        SET_VECTOR_ELT(out, 1, edge_list);
        SET_VECTOR_ELT(out, 0, Rf_ScalarLogical(verdict == 1 ? TRUE : FALSE));
      }
    }
  }

  cvr_close(&s, msg);

  if (what != NULL) {
    if (nprot) UNPROTECT(nprot);
    cvr_stop(what, bad, msg);
  }

  UNPROTECT(nprot);
  return out;
}

/* ------------------------------------------------------------------ *
 * GEOSCoverageUnion
 * ------------------------------------------------------------------ */

SEXP cvr_c_union(SEXP x) {
  R_xlen_t n = Rf_xlength(x);
  cvr_ctx s;
  R_xlen_t bad = 0;
  const char *what = NULL;
  char msg[CVR_ERRBUF];
  SEXP out = R_NilValue;
  int nprot = 0;
  GEOSGeometry *coll;

  if (n == 0) return Rf_allocVector(VECSXP, 0);

  cvr_open(&s);

  coll = cvr_collection(s.ctx, s.reader, x, n, &bad);
  if (coll == NULL) {
    what = "read";
  } else {
    GEOSGeometry *res = GEOSCoverageUnion_r(s.ctx, coll);
    GEOSGeom_destroy_r(s.ctx, coll);
    if (res == NULL) {
      what = "union";
    } else {
      size_t size = 0;
      unsigned char *buf = GEOSWKBWriter_write_r(s.ctx, s.writer, res, &size);
      if (buf == NULL) {
        what = "write";
      } else {
        SEXP raw;
        out = PROTECT(Rf_allocVector(VECSXP, 1)); nprot++;
        raw = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t) size));
        memcpy(RAW(raw), buf, size);
        SET_VECTOR_ELT(out, 0, raw);
        UNPROTECT(1);
        GEOSFree_r(s.ctx, buf);
      }
      GEOSGeom_destroy_r(s.ctx, res);
    }
  }

  cvr_close(&s, msg);

  if (what != NULL) {
    if (nprot) UNPROTECT(nprot);
    cvr_stop(what, bad, msg);
  }

  UNPROTECT(nprot);
  return out;
}

/* ------------------------------------------------------------------ *
 * Vertex count, for measuring what a tolerance actually bought.
 * ------------------------------------------------------------------ */

SEXP cvr_c_n_coord(SEXP x) {
  R_xlen_t n = Rf_xlength(x);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *pout = REAL(out);
  cvr_ctx s;
  R_xlen_t bad = 0;
  char msg[CVR_ERRBUF];
  R_xlen_t i;

  cvr_open(&s);

  for (i = 0; i < n; i++) {
    SEXP item = VECTOR_ELT(x, i);
    GEOSGeometry *g;
    int nc;

    if (item == R_NilValue || TYPEOF(item) != RAWSXP) { bad = i + 1; break; }
    g = GEOSWKBReader_read_r(s.ctx, s.reader, RAW(item), (size_t) Rf_xlength(item));
    if (g == NULL) { bad = i + 1; break; }
    nc = GEOSGetNumCoordinates_r(s.ctx, g);
    GEOSGeom_destroy_r(s.ctx, g);
    pout[i] = (nc < 0) ? NA_REAL : (double) nc;
  }

  cvr_close(&s, msg);

  UNPROTECT(1);
  if (bad > 0) cvr_stop("read", bad, msg);
  return out;
}

/* ------------------------------------------------------------------ *
 * GEOSUnaryUnion.
 *
 * Coverage union is the fast path and refuses overlapping input,
 * which is exactly the state a bad simplification leaves you in. The
 * general union still works there, so gap measurement has something
 * to stand on when the coverage assumption has already failed.
 * ------------------------------------------------------------------ */

SEXP cvr_c_unary_union(SEXP x) {
  R_xlen_t n = Rf_xlength(x);
  cvr_ctx s;
  R_xlen_t bad = 0;
  const char *what = NULL;
  char msg[CVR_ERRBUF];
  SEXP out = R_NilValue;
  int nprot = 0;
  GEOSGeometry *coll;

  if (n == 0) return Rf_allocVector(VECSXP, 0);

  cvr_open(&s);

  coll = cvr_collection(s.ctx, s.reader, x, n, &bad);
  if (coll == NULL) {
    what = "read";
  } else {
    GEOSGeometry *res = GEOSUnaryUnion_r(s.ctx, coll);
    GEOSGeom_destroy_r(s.ctx, coll);
    if (res == NULL) {
      what = "union";
    } else {
      size_t size = 0;
      unsigned char *buf = GEOSWKBWriter_write_r(s.ctx, s.writer, res, &size);
      if (buf == NULL) {
        what = "write";
      } else {
        SEXP raw;
        out = PROTECT(Rf_allocVector(VECSXP, 1)); nprot++;
        raw = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t) size));
        memcpy(RAW(raw), buf, size);
        SET_VECTOR_ELT(out, 0, raw);
        UNPROTECT(1);
        GEOSFree_r(s.ctx, buf);
      }
      GEOSGeom_destroy_r(s.ctx, res);
    }
  }

  cvr_close(&s, msg);

  if (what != NULL) {
    if (nprot) UNPROTECT(nprot);
    cvr_stop(what, bad, msg);
  }

  UNPROTECT(nprot);
  return out;
}

/* ------------------------------------------------------------------ *
 * GEOSTopologyPreserveSimplify, applied to each member independently.
 *
 * The control case. Same library, same Douglas-Peucker family, but
 * each geometry is simplified with no knowledge of its neighbours,
 * which is what opens gaps and overlaps along shared edges. Here so
 * the comparison against cvr_simplify() isolates coverage-awareness
 * as the only variable.
 * ------------------------------------------------------------------ */

SEXP cvr_c_simplify_each(SEXP x, SEXP tolerance) {
  R_xlen_t n = Rf_xlength(x);
  double tol = Rf_asReal(tolerance);
  cvr_ctx s;
  R_xlen_t bad = 0;
  const char *what = NULL;
  char msg[CVR_ERRBUF];
  SEXP out;
  R_xlen_t i;

  out = PROTECT(Rf_allocVector(VECSXP, n));
  cvr_open(&s);

  for (i = 0; i < n; i++) {
    SEXP item = VECTOR_ELT(x, i);
    GEOSGeometry *g, *simp;
    unsigned char *buf;
    size_t size = 0;
    SEXP raw;

    if (item == R_NilValue || TYPEOF(item) != RAWSXP) { bad = i + 1; break; }
    g = GEOSWKBReader_read_r(s.ctx, s.reader, RAW(item), (size_t) Rf_xlength(item));
    if (g == NULL) { bad = i + 1; break; }

    simp = GEOSTopologyPreserveSimplify_r(s.ctx, g, tol);
    GEOSGeom_destroy_r(s.ctx, g);
    if (simp == NULL) { what = "simplify"; break; }

    buf = GEOSWKBWriter_write_r(s.ctx, s.writer, simp, &size);
    GEOSGeom_destroy_r(s.ctx, simp);
    if (buf == NULL) { what = "write"; break; }

    raw = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t) size));
    memcpy(RAW(raw), buf, size);
    GEOSFree_r(s.ctx, buf);
    SET_VECTOR_ELT(out, i, raw);
    UNPROTECT(1);
  }

  cvr_close(&s, msg);

  UNPROTECT(1);
  if (bad > 0) cvr_stop("read", bad, msg);
  if (what != NULL) cvr_stop(what, 0, msg);
  return out;
}

/* ------------------------------------------------------------------ *
 * Area, so a coverage can be checked for gaps and overlaps without
 * reaching for another geometry stack.
 * ------------------------------------------------------------------ */

SEXP cvr_c_area(SEXP x) {
  R_xlen_t n = Rf_xlength(x);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *pout = REAL(out);
  cvr_ctx s;
  R_xlen_t bad = 0;
  char msg[CVR_ERRBUF];
  R_xlen_t i;

  cvr_open(&s);

  for (i = 0; i < n; i++) {
    SEXP item = VECTOR_ELT(x, i);
    GEOSGeometry *g;
    double a = 0;

    if (item == R_NilValue || TYPEOF(item) != RAWSXP) { bad = i + 1; break; }
    g = GEOSWKBReader_read_r(s.ctx, s.reader, RAW(item), (size_t) Rf_xlength(item));
    if (g == NULL) { bad = i + 1; break; }
    pout[i] = GEOSArea_r(s.ctx, g, &a) ? a : NA_REAL;
    GEOSGeom_destroy_r(s.ctx, g);
  }

  cvr_close(&s, msg);

  UNPROTECT(1);
  if (bad > 0) cvr_stop("read", bad, msg);
  return out;
}

/* GEOS version this package was built and linked against. */

SEXP cvr_c_geos_version(void) {
  return Rf_mkString(GEOSversion());
}
