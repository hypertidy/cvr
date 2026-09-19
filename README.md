# cvr

Polygonal coverage operations for R, bound directly to the GEOS C API
(GEOS >= 3.12). Proof of concept.

A **coverage** is a set of valid polygons that are non-overlapping and
edge-matched: vertices along a shared edge are identical in both
neighbours. Administrative boundaries, sea-ice concentration polygons,
land-cover patches: most polygon data people actually simplify is a
coverage, and simplifying it one feature at a time is wrong.

```r
cvr_simplify(x, tolerance, boundary = TRUE)
cvr_is_valid(x, gap_width = 0)
cvr_invalid_edges(x, gap_width = 0)
cvr_union(x)
cvr_n_coord(x)
```

plus `cvr_area()`, `cvr_gap()`, `cvr_n_coord()` for measuring, and
`cvr_simplify_each()` as the row-wise control case.

Input is anything `wk::as_wkb()` accepts. Output is `wk_wkb`, same
length, same order, CRS carried through. Package wk is the interop
provider. 

Read your data however you like:

```r
vapour::vapour_read_geometry(f) |> wk::new_wk_wkb(crs = crs)
gdalraster::GDALVector$new(f, layer)      # WKB per feature
```

## Why this is not a row-wise operator

The GEOS coverage functions take a `GEOMETRYCOLLECTION` and return a
`GEOMETRYCOLLECTION`. They are set-level by construction: you cannot
compute a feature's simplified boundary without knowing who it borders.
That shape does not fit Simple Features at all, which is why PostGIS
had to expose `ST_CoverageSimplify` as a *window function* over a
winset rather than an ordinary scalar function.

`cvr` keeps that shape rather than hiding it. The vector in is the
coverage; the vector out is parallel to it.

## What it buys you

100 NC counties, projected to metres, simplified at 1000 m tolerance.
`cvr_gap()` is the difference between the summed member areas and the
area of their union: if shared edges drift apart, gaps and overlaps
open and the two diverge.

The control is `cvr_simplify_each()`, which is
`GEOSTopologyPreserveSimplify` applied to each feature on its own.
Same library, same Douglas-Peucker family; the only difference is
that it does not know neighbours exist. (It reproduces
`sf::st_simplify()`'s number exactly, which is how it was validated.)

```
cvr_simplify_each(1000):  area gap 2.24e+08 m2   coverage valid: FALSE
cvr_simplify(1000):       area gap 4.58e-05 m2   coverage valid: TRUE
```

224 square kilometres of slivers versus floating-point noise. And the
row-wise output is not merely gappy - it overlaps, badly enough that
`GEOSCoverageUnion` refuses to process it at all. `cvr_gap()` uses a
general union for exactly that reason.

A level-of-detail ladder over the same data, every rung still a valid
coverage:

```
   tol (m)   vertices      pct    valid
         0       2529   100.0%     TRUE
       250       2527    99.9%     TRUE
      1000       2448    96.8%     TRUE
      2500       1647    65.1%     TRUE
      5000       1051    41.6%     TRUE
```

Vertex count is monotone in tolerance, and features never disappear:
GEOS bottoms them out at a triangle rather than dropping them. That
matters for LOD, where a feature vanishing between zoom levels reads
as a bug.

## Inner-only simplification, and why it is the interesting one

`boundary = FALSE` simplifies only the edges shared by two polygons
and leaves the outer boundary of the coverage untouched. Node points
(inner vertices shared by three or more polygons, boundary vertices
shared by two or more) are never moved in either mode.

This is the guarantee that makes chunked or parallel generalization
safe. The usual failure in chunked tiling is that each chunk is
generalized against only its own neighbours, so geometry that survives
with full context collapses when processed alone and chunk seams stop
lining up. Pin each chunk's outer edges and the problem goes away:

```
chunks: 40 west, 60 east        (simplified independently, boundary = FALSE)

rejoined is a valid coverage: TRUE
rejoined area gap:            9.16e-05 m2     (floating point noise)
chunk seam preserved:         TRUE

full-context:  2448 vertices
inner-only chunked: 2463 vertices   (+0.6%)
```

Pinning the seams costs under one percent of the vertex budget on this
data, and buys embarrassingly parallel generalization with exact
reassembly.

See `inst/examples/demo-ladder.R` for the script that produces all of
the above. It runs against `inst/extdata/nc-32119.csv.gz`, a WKT
fixture written once by ogr2ogr, so the demo has no reader dependency
of any kind.

## Install

Needs GEOS 3.12.0 or later, where the coverage API landed.

```
apt-get install libgeos-dev     # or: brew install geos, dnf install geos-devel
R CMD INSTALL cvr
```

`configure` finds GEOS via `geos-config` and refuses to build against
anything older than 3.12. Set `GEOS_CONFIG` to override discovery.

## Notes

- GEOS takes `preserveBoundary`; this package exposes `boundary`, the
  inverse, matching PostGIS `simplifyBoundary` and shapely
  `simplify_boundary`. Coming from either of those, the argument means
  what you expect.
- `cvr_is_valid()` returns one verdict for the coverage as a whole.
  `cvr_invalid_edges()` returns the detail: a MULTILINESTRING of
  offending edges per member, empty where that member is a valid
  participant. Run it before simplifying, since an invalid coverage is
  still simplified but the artifacts are silent.
- An invalid coverage produces no error, just bad output. That is
  GEOS behaviour, not a wrapper choice.

## Status and what is missing

Proof of concept. Deliberately not filed against paleolimbot/geos,
which is a larger change to a package on a slower cycle.

The gap this does not close: GEOS gives you a coverage at *one*
tolerance. A multi-resolution store wants the vertex *ranking* rather
than the thresholded result. Visvalingam is repeated minimum-area
vertex removal, so every vertex has an effective area, and with that
column in hand any tolerance is a predicate rather than a
recomputation. GEOS does not expose it. Getting at it means noding the
coverage into an arc table and ranking vertices per arc, at which
point a feature becomes an ordered list of arc references and a level
of detail becomes `WHERE rank >= f(tolerance)`.

That prototype now exists and is measured against this package: same
inputs, same ladder, same three checks per rung. Quality matches GEOS
(identical vertex counts through 2500 m, within 0.5% beyond), the
levels come out strictly nested, and the threshold runs about 11x
faster per rung at 200k vertices. See the arc-table prototype.
