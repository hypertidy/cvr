library(cvr)
library(wk)

ok <- function(cond, what) {
  if (!isTRUE(cond)) stop("FAILED: ", what, call. = FALSE)
  cat("ok  ", what, "\n", sep = "")
}

# The coverage from the PostGIS ST_CoverageSimplify documentation.
cov <- as_wkb(wkt(c(
  paste0("POLYGON ((160 150, 110 130, 90 100, 90 70, 60 60, 50 10, 30 30,",
         " 40 50, 25 40, 10 60, 30 100, 30 120, 20 170, 60 180, 90 190,",
         " 130 180, 130 160, 160 150),",
         " (40 160, 50 140, 66 125, 60 100, 80 140, 90 170, 60 160, 40 160))"),
  "POLYGON ((40 160, 60 160, 90 170, 80 140, 60 100, 66 125, 50 140, 40 160))",
  "POLYGON ((110 130, 160 50, 140 50, 120 33, 90 30, 50 10, 60 60, 90 70, 90 100, 110 130))",
  "POLYGON ((160 150, 150 120, 160 90, 160 50, 110 130, 160 150))"
)))

## ---- structural guarantees ----------------------------------------

simp <- cvr_simplify(cov, 30)

ok(inherits(simp, "wk_wkb"), "cvr_simplify returns a wk_wkb")
ok(length(simp) == length(cov), "result is parallel to input")
ok(cvr_is_valid(cov), "input is a valid coverage")
ok(cvr_is_valid(simp), "simplified result is still a valid coverage")
ok(sum(cvr_n_coord(simp)) < sum(cvr_n_coord(cov)), "vertices were actually removed")

# GEOS documents this: features never disappear, they bottom out at a
# triangle. A ridiculous tolerance must not drop a row.
huge <- cvr_simplify(cov, 1e6)
ok(length(huge) == length(cov), "features survive an absurd tolerance")
ok(!any(vapply(unclass(huge), is.null, logical(1))), "no NULL members")

## ---- the property the whole scheme rests on -----------------------
## Simplifying a coverage must not change its total area by opening
## gaps or overlaps between neighbours. cvr_gap() is exactly that
## test: summed member areas against the area of the union.

ok(cvr_gap(cov) < 1e-6, "input: member areas sum to the union area")
ok(cvr_gap(simp) < 1e-6, "simplified: member areas still sum to the union area")

indep <- cvr_simplify_each(cov, 30)
ok(length(indep) == length(cov), "cvr_simplify_each is parallel to input")

## ---- inner-only simplification ------------------------------------
## The outer boundary must come back byte-identical, which is what
## lets a subset of a coverage be simplified on its own and still fit
## against the rest.

inner <- cvr_simplify(cov, 30, boundary = FALSE)
ok(cvr_is_valid(inner), "inner-only result is a valid coverage")

ok(identical(unclass(cvr_union(cov)), unclass(cvr_union(inner))),
   "boundary = FALSE leaves the outer boundary byte-identical")
ok(sum(cvr_n_coord(inner)) < sum(cvr_n_coord(cov)), "inner edges were still simplified")
ok(sum(cvr_n_coord(inner)) > sum(cvr_n_coord(simp)), "inner-only keeps more vertices than full")

## ---- monotone ladder ----------------------------------------------
## A level-of-detail ladder needs vertex count to fall monotonically
## with tolerance, and every rung to remain a valid coverage.

tol <- c(0, 1, 5, 10, 20, 40, 80)
nv <- vapply(tol, function(t) sum(cvr_n_coord(cvr_simplify(cov, t))), numeric(1))
valid <- vapply(tol, function(t) cvr_is_valid(cvr_simplify(cov, t)), logical(1))

ok(all(diff(nv) <= 0), "vertex count is monotone non-increasing in tolerance")
ok(all(valid), "every rung of the ladder is a valid coverage")

## ---- validation and invalid edges ---------------------------------

overlapping <- as_wkb(wkt(c(
  "POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))",
  "POLYGON ((5 0, 15 0, 15 10, 5 10, 5 0))"
)))
ok(!cvr_is_valid(overlapping), "overlaps are detected as invalid")

edges <- cvr_invalid_edges(overlapping)
ok(length(edges) == length(overlapping), "invalid edges are parallel to input")
ok(inherits(edges, "wk_wkb"), "invalid edges come back as wk_wkb")

valid_edges <- cvr_invalid_edges(cov)
ok(length(valid_edges) == length(cov), "valid coverage still gets a parallel edge result")

## ---- plumbing ------------------------------------------------------

crs_cov <- as_wkb(wkt(c("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"), crs = "EPSG:3031"))
ok(identical(wk_crs(cvr_simplify(crs_cov, 0.1)), "EPSG:3031"), "CRS is carried through")

ok(length(cvr_simplify(as_wkb(wkt(character(0))), 1)) == 0, "zero-length input is allowed")
ok(inherits(try(cvr_simplify(cov, c(1, 2)), silent = TRUE), "try-error"),
   "non-scalar tolerance is rejected")
ok(inherits(try(cvr_simplify(as_wkb(wkt(NA_character_)), 1), silent = TRUE), "try-error"),
   "missing geometry is rejected")

ok(length(cvr_union(cov)) == 1, "cvr_union returns a single geometry")
ok(abs(sum(cvr_area(cov)) - sum(cvr_area(cvr_union(cov)))) < 1e-6,
   "cvr_area agrees across members and union")
ok(!("package:sf" %in% search()), "no sf was needed")
ok(is.character(cvr_geos_version()), "GEOS version reports")

## ---- a real coverage ----------------------------------------------
## The four-polygon example above is too coarse to show the failure:
## Douglas-Peucker happens to drop the same vertices from both sides
## of every shared edge. On real data it does not.

f <- system.file("extdata/nc-32119.csv.gz", package = "cvr")
if (nzchar(f)) {
  nc <- as_wkb(wkt(sub('^"(.*)",[^,]*$', "\\1", readLines(f, warn = FALSE)[-1]),
                   crs = "EPSG:32119"))

  ok(length(nc) == 100, "NC fixture reads")
  ok(cvr_is_valid(nc), "NC fixture is a valid coverage")

  nc_cov <- cvr_simplify(nc, 1000)
  nc_each <- cvr_simplify_each(nc, 1000)

  ok(cvr_is_valid(nc_cov), "coverage simplify keeps it a coverage")
  ok(!cvr_is_valid(nc_each), "row-wise simplify destroys the coverage")
  ok(cvr_gap(nc_cov) < 1e-3, "coverage simplify leaves no area gap")
  ok(cvr_gap(nc_each) > 1e6, "row-wise simplify opens a large area gap")

  # inner-only lets chunks be simplified independently and rejoined
  half <- seq_len(50)
  re <- nc
  re[half] <- cvr_simplify(nc[half], 1000, boundary = FALSE)
  re[-half] <- cvr_simplify(nc[-half], 1000, boundary = FALSE)
  ok(cvr_is_valid(re), "independently simplified chunks rejoin as a coverage")
  ok(cvr_gap(re) < 1e-3, "rejoined chunks leave no area gap")
}

cat("\nAll tests passed against GEOS ", cvr_geos_version(), "\n", sep = "")
