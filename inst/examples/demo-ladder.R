# A level-of-detail ladder over a real coverage, and the property that
# makes chunked generalization safe.
#
# NC counties: 100 adjacent polygons in NC State Plane metres, so the
# tolerance reads as ground distance. The fixture is WKT written once
# by ogr2ogr; in real use read straight to WKB:
#
#   vapour::vapour_read_geometry(f) |> wk::new_wk_wkb(crs = crs)
#
# wk and cvr only. No sf.

library(cvr)
library(wk)

f <- system.file("extdata/nc-32119.csv.gz", package = "cvr")
cov <- as_wkb(wkt(sub('^"(.*)",[^,]*$', "\\1", readLines(f, warn = FALSE)[-1]),
                  crs = "EPSG:32119"))

cat("GEOS:", cvr_geos_version(), "\n")
cat("features:", length(cov), " vertices:", sum(cvr_n_coord(cov)), "\n")
cat("valid coverage:", cvr_is_valid(cov), "\n\n")

# --- the ladder ------------------------------------------------------
# Screen-space error is the decision variable: pick the tolerance from
# the ground distance one pixel covers at each zoom, not from a feature
# count.
#
# cvr_gap() is summed member areas minus the area of their union. Zero
# for a coverage; it grows the moment shared edges drift apart.

tol <- c(0, 25, 50, 100, 250, 500, 1000, 2500, 5000)
base_n <- sum(cvr_n_coord(cov))

cat(sprintf("%10s %10s %8s %8s %14s\n",
            "tol (m)", "vertices", "pct", "valid", "area gap (m2)"))
for (t in tol) {
  s <- cvr_simplify(cov, t)
  n <- sum(cvr_n_coord(s))
  cat(sprintf("%10.0f %10.0f %7.1f%% %8s %14.2e\n",
              t, n, 100 * n / base_n, cvr_is_valid(s), cvr_gap(s)))
}

# For contrast: the same GEOS library, the same Douglas-Peucker
# family, applied to each feature independently. The only difference
# is that it does not know neighbours exist.
cat(sprintf("\nsimplify_each(1000): area gap %.2e m2\n",
            cvr_gap(cvr_simplify_each(cov, 1000))))
cat(sprintf("cvr_simplify(1000):  area gap %.2e m2\n",
            cvr_gap(cvr_simplify(cov, 1000))))

# --- chunked reassembly ----------------------------------------------
# The failure mode in chunked tiling is that each chunk is generalized
# against only its own neighbours, so geometry that survives with full
# context collapses when processed alone, and chunk seams no longer
# line up.
#
# Inner-only simplification is the fix: pin each chunk's outer edges,
# simplify only edges interior to the chunk. Independently simplified
# chunks then reassemble into a valid coverage.

cx <- vapply(seq_along(cov), function(i) {
  b <- unclass(wk_bbox(cov[i]))
  (b$xmin + b$xmax) / 2
}, numeric(1))
west <- cx < 500000
cat(sprintf("\nchunks: %d west, %d east\n", sum(west), sum(!west)))

chunk_w <- cvr_simplify(cov[west], 1000, boundary = FALSE)
chunk_e <- cvr_simplify(cov[!west], 1000, boundary = FALSE)

rejoined <- cov
rejoined[west] <- chunk_w
rejoined[!west] <- chunk_e

cat("rejoined is a valid coverage:", cvr_is_valid(rejoined), "\n")
cat(sprintf("rejoined area gap: %.2e m2\n", cvr_gap(rejoined)))
cat(sprintf("rejoined vertices: %.0f (%.1f%% of original)\n",
            sum(cvr_n_coord(rejoined)),
            100 * sum(cvr_n_coord(rejoined)) / base_n))

# The outer boundary of each chunk comes back byte-identical, so a
# chunk can be re-simplified later without disturbing its neighbours.
cat("chunk seam preserved:",
    identical(unclass(cvr_union(cov[west])), unclass(cvr_union(chunk_w))), "\n")

# The cost of pinning: inner-only keeps the seams over-detailed.
full <- cvr_simplify(cov, 1000)
cat(sprintf("\nfull-context: %.0f vertices, inner-only chunked: %.0f (+%.1f%%)\n",
            sum(cvr_n_coord(full)), sum(cvr_n_coord(rejoined)),
            100 * (sum(cvr_n_coord(rejoined)) / sum(cvr_n_coord(full)) - 1)))

cat("\nno sf on the search path:", !("package:sf" %in% search()), "\n")
