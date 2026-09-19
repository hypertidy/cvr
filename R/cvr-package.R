# The @useDynLib directive below is load-bearing, and its absence
# fails in a way that looks nothing like its cause.
#
# NAMESPACE is regenerated from these sources, so if nothing here
# declares the DLL the first document() drops the useDynLib line. The
# package then builds, installs and loads perfectly cleanly, and every
# .Call() fails with:
#
#   object 'cvr_c_simplify' not found
#
# .registration = TRUE binds each routine registered by
# R_registerRoutines as an object named after the routine, which is
# what the .Call()s in this package pass. R_forceSymbols(TRUE) in
# src/init.c then refuses lookup by character name, so a missing
# directive cannot be papered over with .Call("cvr_c_simplify", ...).
#
# Keep the directive on its own line: roxygen tags run to the next
# tag, so trailing prose inside the block becomes part of the value.

#' @keywords internal
#' @useDynLib cvr, .registration = TRUE
"_PACKAGE"
