#include <R.h>
#include <Rinternals.h>

extern "C" {

SEXP MBAPoints(SEXP xyz, SEXP xyzEst, SEXP m, SEXP n, SEXP h, SEXP extend, SEXP verbose);

SEXP MBASurf(SEXP xyz, SEXP noX, SEXP noY, SEXP m, SEXP n, SEXP h, SEXP extend, SEXP hpts,
	     SEXP xMinDom, SEXP xMaxDom, SEXP yMinDom, SEXP yMaxDom);

}
