#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "xcms.h"
#include "util.h"
#include "binners.h"
/*
 * Contains binning utils.
 */


/*
 * ----------------------- R ENTRY POINTS -----------------------
 */

/*
 * Define same-sized bins on y and select the max value in x corresponding to
 * y-values within each bin.
 * Binning is defined based on the number of bins (nBin) and the range of values
 * in x that should be binned (fromX to toX).
 * This binning corresponds to seq(fromX, toX, length.out = (nBins + 1))
 */
SEXP binXonY_nBins_max(SEXP x, SEXP y, SEXP nBins, SEXP fromX, SEXP toX) {
  SEXP ans;
  int n_bin = asInteger(nBins);
  double from_x, to_x;
  from_x = REAL(fromX)[0];
  to_x = REAL(toX)[0];

  // Check that n_bin > 0

  // Create output
  PROTECT(ans = NEW_NUMERIC(n_bin + 1));
  // Calculate the breaks
  _breaks_on_nBins(from_x, to_x, n_bin, REAL(ans));
  UNPROTECT(1);
  return ans;
}

/*
 * Same as binXonY_nBins_max, but the binning is defined by the binSize.
 */
SEXP binXonY_binSize_max(SEXP x, SEXP y, SEXP binSize, SEXP fromX, SEXP toX) {
  SEXP ans;
  int n_bin;
  double bin_size, from_x, to_x;
  bin_size = REAL(binSize)[0];
  from_x = REAL(fromX)[0];
  to_x = REAL(toX)[0];
  n_bin = (int)ceil((to_x - from_x) / bin_size);
  // Create output
  PROTECT(ans = NEW_NUMERIC(n_bin + 1));
  // Calculate breaks
  _breaks_on_binSize(from_x, to_x, n_bin, bin_size, REAL(ans));
  UNPROTECT(1);
  return ans;
}

/*
 * ----------------------- INTERNAL FUNCTIONS -----------------------
 */

/*
 * Create breaks for binning: seq(from_val, to_val, length.out = (n_bin + 1))
 */
void _breaks_on_nBins(double from_val, double to_val, int n_bin,
		      double *brks) {
  int i;
  double current_val, bin_size;

  bin_size = (to_val - from_val) / (double)n_bin;
  current_val = from_val;
  for (i = 0; i <= n_bin; i++) {
    brks[i] = current_val;
    current_val += bin_size;
  }
  return;
}

/*
 * Create breaks for binning: seq(from_val, to_val, by = bin_size)
 */
void _breaks_on_binSize(double from_val, double to_val, int n_bin,
			double bin_size, double *brks) {
  // We have to make a decision here, the max break should be to_val!
  // no of breaks: ceil(from_val - to_val / bin_size)
  // We have to assume that the *brks array has already been sized to the
  // correct length (i.e. ceil((from_val - to_val) / bin_size))
  double current_val;
  current_val = from_val;
  for(int i = 0; i < n_bin; i++) {
    brks[i] = current_val;
    current_val += bin_size;
  }
  brks[n_bin] = to_val;
  return;
}



/*
 * Some simple functions to check passing of arguments.
 */
SEXP test_integer(SEXP x) {
  int x_val = asInteger(x);
  Rprintf("input asInteger(x): %d\n", x_val);

  //
  int *p_ans;
  int *p_x = INTEGER(x);
  Rprintf("getting the first value from the pointer: %d\n", p_x[0]);

  SEXP ans = allocVector(INTSXP, LENGTH(x));
  p_ans = INTEGER(ans);
  p_ans[0] = x_val;
  return ans;
}

SEXP test_real(SEXP x) {
  int x_val = asReal(x);
  Rprintf("input asReal(x): %f\n", x_val);

  //
  double *p_ans;
  double *p_x = REAL(x);
  Rprintf("getting the first value from the pointer: %f\n", p_x[0]);

  SEXP ans = allocVector(REALSXP, LENGTH(x));
  p_ans = REAL(ans);
  p_ans[0] = x_val;
  return ans;
}


