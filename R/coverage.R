#' Polygonal coverage operations
#'
#' Bindings to the GEOS coverage API (GEOS >= 3.12). A polygonal
#' coverage is a set of valid polygons that are non-overlapping and
#' edge-matched: vertices along a shared edge are identical in both
#' neighbours.
#'
#' These are set-level operators, not row-wise ones. They consume the
#' whole vector as a single coverage and return a result parallel to
#' it. This is why PostGIS exposes `ST_CoverageSimplify` as a window
#' function rather than an ordinary one, and why the result of
#' `cvr_simplify()` cannot be computed feature by feature.
#'
#' `cvr_simplify()` applies a Visvalingam-Whyatt simplification to the
#' coverage edges. Shared edges are simplified once and stay identical
#' in both neighbours, so no gaps or overlaps open up. Node points
#' (inner vertices shared by three or more polygons, boundary vertices
#' shared by two or more) are never moved. Features never disappear,
#' though they may be reduced to a triangle.
#'
#' With `boundary = FALSE` only the inner, shared edges are simplified
#' and the outer boundary of the coverage is left untouched. A subset
#' simplified this way still fits exactly against the rest of the
#' coverage, which is what makes chunked or incremental generalization
#' possible.
#'
#' @param x A vector of polygons: anything [wk::as_wkb()] accepts, such
#'   as a `wk_wkb`, `wk_wkt`, `sfc` or `sf` object. All members must be
#'   POLYGON or MULTIPOLYGON.
#' @param tolerance Simplification tolerance in the units of `x`,
#'   roughly the square root of the area of the triangles removed.
#'   Length 1.
#' @param boundary Simplify the outer boundary of the coverage as well
#'   as the inner shared edges? `FALSE` simplifies inner edges only.
#' @param gap_width The maximum width of gaps to report as invalid.
#'
#' @return
#' `cvr_simplify()` a `wk_wkb` the same length as `x`, in the same
#' order, carrying the CRS of `x`.
#'
#' `cvr_is_valid()` a single `TRUE` or `FALSE` for the coverage as a
#' whole.
#'
#' `cvr_invalid_edges()` a `wk_wkb` the same length as `x`: a
#' MULTILINESTRING of the offending edges for each invalid member, and
#' an empty geometry where that member is a valid participant.
#'
#' `cvr_union()` a length-1 `wk_wkb`.
#'
#' `cvr_n_coord()` a numeric vector of vertex counts, one per member.
#'
#' `cvr_area()` a numeric vector of planar areas, one per member.
#'
#' `cvr_simplify_each()` a `wk_wkb` the same length as `x`, each
#' member simplified independently with
#' `GEOSTopologyPreserveSimplify`. The control case, not the one to
#' use: same library, same algorithm family, but no knowledge of
#' neighbours. It is here so the cost of ignoring them can be
#' measured rather than asserted.
#'
#' `cvr_gap()` a single number: the absolute difference between the
#' summed member areas and the area of their union. Zero for a
#' coverage. Grows as soon as shared edges drift apart, so it is the
#' direct test of whether a simplification stayed edge-matched.
#'
#' @examples
#' library(wk)
#'
#' # The four-polygon coverage from the PostGIS ST_CoverageSimplify docs
#' cov <- as_wkb(wkt(c(
#'   paste0("POLYGON ((160 150, 110 130, 90 100, 90 70, 60 60, 50 10, 30 30,",
#'          " 40 50, 25 40, 10 60, 30 100, 30 120, 20 170, 60 180, 90 190,",
#'          " 130 180, 130 160, 160 150),",
#'          " (40 160, 50 140, 66 125, 60 100, 80 140, 90 170, 60 160, 40 160))"),
#'   "POLYGON ((40 160, 60 160, 90 170, 80 140, 60 100, 66 125, 50 140, 40 160))",
#'   "POLYGON ((110 130, 160 50, 140 50, 120 33, 90 30, 50 10, 60 60, 90 70, 90 100, 110 130))",
#'   "POLYGON ((160 150, 150 120, 160 90, 160 50, 110 130, 160 150))"
#' )))
#'
#' cvr_is_valid(cov)
#' sum(cvr_n_coord(cov))
#'
#' simp <- cvr_simplify(cov, 30)
#' sum(cvr_n_coord(simp))
#'
#' # still a coverage: shared edges moved together
#' cvr_is_valid(simp)
#'
#' # inner-only simplification pins the outer boundary
#' inner <- cvr_simplify(cov, 30, boundary = FALSE)
#' cvr_is_valid(inner)
#'
#' @name cvr_coverage
NULL

#' @rdname cvr_coverage
#' @export
cvr_simplify <- function(x, tolerance, boundary = TRUE) {
  wkb <- cvr_as_wkb(x)
  tolerance <- cvr_scalar_double(tolerance, "tolerance")
  boundary <- cvr_scalar_logical(boundary, "boundary")

  # GEOS takes preserveBoundary, the inverse of the PostGIS and shapely
  # simplify_boundary argument. Follow the latter, which is what people
  # coming from those tools expect.
  out <- .Call(cvr_c_simplify, unclass(wkb), tolerance, !boundary)
  cvr_restore(out, wkb)
}

#' @rdname cvr_coverage
#' @export
cvr_is_valid <- function(x, gap_width = 0) {
  wkb <- cvr_as_wkb(x)
  gap_width <- cvr_scalar_double(gap_width, "gap_width")
  .Call(cvr_c_is_valid, unclass(wkb), gap_width)[[1L]]
}

#' @rdname cvr_coverage
#' @export
cvr_invalid_edges <- function(x, gap_width = 0) {
  wkb <- cvr_as_wkb(x)
  gap_width <- cvr_scalar_double(gap_width, "gap_width")
  res <- .Call(cvr_c_is_valid, unclass(wkb), gap_width)
  edges <- res[[2L]]

  # A valid coverage reports no edges at all; return one empty geometry
  # per member so the result stays parallel to the input.
  if (length(edges) == 0L && length(wkb) > 0L) {
    edges <- unclass(wk::as_wkb(wk::wkt(rep("MULTILINESTRING EMPTY", length(wkb)))))
  }
  cvr_restore(edges, wkb)
}

#' @rdname cvr_coverage
#' @export
cvr_union <- function(x) {
  wkb <- cvr_as_wkb(x)
  cvr_restore(.Call(cvr_c_union, unclass(wkb)), wkb)
}

#' @rdname cvr_coverage
#' @export
cvr_n_coord <- function(x) {
  .Call(cvr_c_n_coord, unclass(cvr_as_wkb(x)))
}

#' @rdname cvr_coverage
#' @export
cvr_simplify_each <- function(x, tolerance) {
  wkb <- cvr_as_wkb(x)
  tolerance <- cvr_scalar_double(tolerance, "tolerance")
  cvr_restore(.Call(cvr_c_simplify_each, unclass(wkb), tolerance), wkb)
}

#' @rdname cvr_coverage
#' @export
cvr_area <- function(x) {
  .Call(cvr_c_area, unclass(cvr_as_wkb(x)))
}

#' @rdname cvr_coverage
#' @export
cvr_gap <- function(x) {
  wkb <- cvr_as_wkb(x)
  # Unary union, not coverage union: this is measured precisely when
  # the coverage assumption may already have failed, and coverage
  # union refuses overlapping input.
  u <- cvr_restore(.Call(cvr_c_unary_union, unclass(wkb)), wkb)
  abs(sum(cvr_area(wkb)) - sum(cvr_area(u)))
}

#' GEOS version cvr was built against
#'
#' @return A version string.
#' @export
#' @examples
#' cvr_geos_version()
cvr_geos_version <- function() {
  .Call(cvr_c_geos_version)
}

# helpers -----------------------------------------------------------

cvr_as_wkb <- function(x) {
  wkb <- wk::as_wkb(x)
  if (anyNA(unclass(wkb))) {
    stop("`x` must not contain missing geometries: a coverage has no room for them")
  }
  wkb
}

cvr_restore <- function(lst, template) {
  wk::new_wk_wkb(lst, crs = wk::wk_crs(template))
}

cvr_scalar_double <- function(x, arg) {
  x <- as.double(x)
  if (length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single non-missing number", arg))
  }
  x
}

cvr_scalar_logical <- function(x, arg) {
  x <- as.logical(x)
  if (length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single TRUE or FALSE", arg))
  }
  x
}
