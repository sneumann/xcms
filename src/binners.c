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
 * Arguments:
 * x: numeric vector of values on which y should be binned.
 * y: numeric vector that should be binned.
 * nBins: number of bins.
 * fromX: the lowest x-value form which binning should start.
 * toX: the largest x-value to be included in the binning.
 * fromIdx, toIdx: indices in array x (0-based) that allow to specify a sub-set
 *     of x/y on which the binning should be done.
 * shiftByHalfBinSize: either 0 or 1. If 0 breaks are defined from fromX to toX,
 *     if 1 breaks are defined from fromX - (bin_size / 2) to toX + (bin_size / 2).
 *     The bin_size is calculated as (toX - fromX) / (nBins - 1), such that the
 *     number of bins corresponds to nBins. If 1 binning is performed similarly to
 *     the profBin method.
 * The function returns a list with the first element (x) being the bin
 * midpoints and the second element (y) the max y-value within each bin.
 */
SEXP binXonY_nBins_max(SEXP x, SEXP y, SEXP nBins, SEXP fromX, SEXP toX,
		       SEXP fromIdx, SEXP toIdx, SEXP shiftByHalfBinSize) {
  SEXP ans, brks, bin_mids, ans_list, names;
  int n_bin = asInteger(nBins);
  int from_idx = asInteger(fromIdx);
  int to_idx = asInteger(toIdx);
  double from_x, to_x, *p_ans, *p_brks;
  from_x = REAL(fromX)[0];
  to_x = REAL(toX)[0];

  /* Eventually shift by a half bin_size */
  int shift_by_half_bin_size = asInteger(shiftByHalfBinSize);

  /* Error checks */
  /* Check that n_bin > 0 */
  if (n_bin <= 0)
    error("'nBins' must be larger 1!");
  // Check that from_x to_x are > 0 and < length(x)
  if (from_idx < 0 | to_idx < 0)
    error("'fromIdx' and 'toIdx' have to be >= 0!");
  if (from_idx > to_idx)
    error("'fromIdx' has to be smaller than 'toIdx'!");
  if (to_idx >= LENGTH(x))
    error("'toIdx' can not be larger than length(x)!");

  /* Create output */
  PROTECT(brks = allocVector(REALSXP, n_bin + 1));
  PROTECT(ans = allocVector(REALSXP, n_bin));
  /* Calculate the breaks */
  p_brks = REAL(brks);
  _breaks_on_nBins(from_x, to_x, n_bin, p_brks, shift_by_half_bin_size);
  /* Do the binning. */
  p_ans = REAL(ans);
  for(int i = 0; i < n_bin; i++) {
    p_ans[i] = 0;
  }
  /* max*/
  _bin_y_on_x_with_breaks_max(REAL(x), REAL(y), p_brks, p_ans, n_bin,
			      from_idx, to_idx);
  /* Calculate bin mid-points */
  PROTECT(bin_mids = allocVector(REALSXP, n_bin));
  _bin_midPoint(p_brks, REAL(bin_mids), n_bin);
  /* Now create the result list. */
  PROTECT(ans_list = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ans_list, 0, bin_mids);
  SET_VECTOR_ELT(ans_list, 1, ans);
  /* Setting names */
  PROTECT(names = allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("x"));
  SET_STRING_ELT(names, 1, mkChar("y"));
  setAttrib(ans_list, R_NamesSymbol, names);
  UNPROTECT(5);
  return ans_list;
  /* UNPROTECT(2); */
  /* return ans; */
}

/*
 * Same as binXonY_nBins_max, but the binning is defined by the binSize.
 */
SEXP binXonY_binSize_max(SEXP x, SEXP y, SEXP binSize, SEXP fromX, SEXP toX,
			 SEXP fromIdx, SEXP toIdx, SEXP shiftByHalfBinSize) {
  SEXP ans, brks, bin_mids, ans_list, names;
  int n_bin;
  int from_idx = asInteger(fromIdx);
  int to_idx = asInteger(toIdx);
  double bin_size, from_x, to_x, *p_ans, *p_brks;
  from_x = REAL(fromX)[0];
  to_x = REAL(toX)[0];
  bin_size = REAL(binSize)[0];
  /* Eventually shift by a half bin_size */
  int shift_by_half_bin_size = asInteger(shiftByHalfBinSize);
  if (shift_by_half_bin_size > 0) {
    from_x = from_x - bin_size / 2;
  }

  n_bin = (int)ceil((to_x - from_x) / bin_size);
  // Input checking
  if (bin_size < 0)
    error("'binSize' has to be > 0!");
  if (from_idx < 0 | to_idx < 0)
    error("'fromIdx' and 'toIdx' have to be >= 0!");
  if (from_idx > to_idx)
    error("'fromIdx' has to be smaller than 'toIdx'!");
  if (to_idx >= LENGTH(x))
    error("'toIdx' can not be larger than length(x)!");

  // Create output
  PROTECT(brks = allocVector(REALSXP, n_bin + 1));
  PROTECT(ans = allocVector(REALSXP, n_bin));
  // Calculate breaks
  p_brks = REAL(brks);
  _breaks_on_binSize(from_x, to_x, n_bin, bin_size, REAL(brks));
  // Do the binning.
  p_ans = REAL(ans);
  for(int i = 0; i < n_bin; i++) {
    p_ans[i] = 0;
  }
  _bin_y_on_x_with_breaks_max(REAL(x), REAL(y), p_brks, p_ans, n_bin,
			      from_idx, to_idx);
  // Calculate bin mid-points
  PROTECT(bin_mids = allocVector(REALSXP, n_bin));
  _bin_midPoint(p_brks, REAL(bin_mids), n_bin);
  // Now create the result list.
  PROTECT(ans_list = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ans_list, 0, bin_mids);
  SET_VECTOR_ELT(ans_list, 1, ans);
  // Setting names
  PROTECT(names = allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("x"));
  SET_STRING_ELT(names, 1, mkChar("y"));
  setAttrib(ans_list, R_NamesSymbol, names);
  UNPROTECT(5);
  return ans_list;
  /* UNPROTECT(2); */
  /* return ans; */
}

/*
 * Same as binXonY_nBins_max, but with pre-defined breaks.
 */
SEXP binXonY_breaks_max(SEXP x, SEXP y, SEXP breaks, SEXP fromIdx,
			SEXP toIdx) {
  SEXP ans, bin_mids, ans_list, names;
  int n_bin;
  int from_idx = asInteger(fromIdx);
  int to_idx = asInteger(toIdx);
  double bin_size, from_x, to_x, *p_ans, *p_brks;
  n_bin = (LENGTH(breaks) - 1);
  p_brks = REAL(breaks);

  // Input checking
  if (n_bin < 1)
    error("Not enough breaks defined!");
  if (from_idx < 0 | to_idx < 0)
    error("'fromIdx' and 'toIdx' have to be >= 0!");
  if (from_idx > to_idx)
    error("'fromIdx' has to be smaller than 'toIdx'!");
  if (to_idx >= LENGTH(x))
    error("'toIdx' can not be larger than length(x)!");

  // Create output
  PROTECT(ans = allocVector(REALSXP, n_bin));
  // Do the binning.
  p_ans = REAL(ans);
  for(int i = 0; i < n_bin; i++) {
    p_ans[i] = 0;
  }
  _bin_y_on_x_with_breaks_max(REAL(x), REAL(y), p_brks, p_ans, n_bin,
			      from_idx, to_idx);
  // Calculate bin mid-points
  PROTECT(bin_mids = allocVector(REALSXP, n_bin));
  _bin_midPoint(p_brks, REAL(bin_mids), n_bin);
  // Now create the result list.
  PROTECT(ans_list = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ans_list, 0, bin_mids);
  SET_VECTOR_ELT(ans_list, 1, ans);
  // Setting names
  PROTECT(names = allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("x"));
  SET_STRING_ELT(names, 1, mkChar("y"));
  setAttrib(ans_list, R_NamesSymbol, names);
  UNPROTECT(4);
  return ans_list;
  /* UNPROTECT(1); */
  /* return ans; */
}


/*
 * Some other test functions.
 */

/*
 * Get breaks given from x, to x and number of bins.
 */
SEXP breaks_on_nBins(SEXP fromX, SEXP toX, SEXP nBins) {
  SEXP ans;
  int n_bin;
  double from_x, to_x;
  n_bin = asInteger(nBins);
  from_x = REAL(fromX)[0];
  to_x = REAL(toX)[0];
  PROTECT(ans = allocVector(REALSXP, n_bin + 1));
  _breaks_on_nBins(from_x, to_x, n_bin, REAL(ans), 0);
  UNPROTECT(1);
  return ans;
}

/*
 * Get breaks given from x, to x and number of bins.
 */
SEXP breaks_on_binSize(SEXP fromX, SEXP toX, SEXP binSize) {
  SEXP ans;
  int n_bin;
  double from_x, to_x, bin_size;
  bin_size = REAL(binSize)[0];
  from_x = REAL(fromX)[0];
  to_x = REAL(toX)[0];
  n_bin = (int)ceil((to_x - from_x) / bin_size);
  PROTECT(ans = allocVector(REALSXP, n_bin + 1));
  _breaks_on_binSize(from_x, to_x, n_bin, bin_size, REAL(ans));
  UNPROTECT(1);
  return ans;
}

/*
 * ----------------------- INTERNAL FUNCTIONS -----------------------
 */

/*
 * Create breaks for binning: seq(from_val, to_val, length.out = (n_bin + 1))
 * shift_by_half_bin_size: either 0 or 1, if 1 the bin mid-points will be shifted left
 * by half of the bin_size.
 */
void _breaks_on_nBins(double from_val, double to_val, int n_bin,
		      double *brks, int shift_by_half_bin_size) {
  int i;
  double current_val, bin_size;

  /* If we're going to shift the bin mid-points we have to ensure that the bin_size
   * will be slightly larger to include also the to_val (+ bin_size/2).
   */
  bin_size = (to_val - from_val) / ((double)n_bin - (double)shift_by_half_bin_size);
  if (shift_by_half_bin_size > 0) {
    current_val = from_val - (bin_size / 2);
  } else {
    current_val = from_val;
  }
  for (i = 0; i <= n_bin; i++) {
    brks[i] = current_val;
    current_val += bin_size;
  }
  return;
}

/*
 * Create breaks for binning: seq(from_val, to_val, by = bin_size).
 */
void _breaks_on_binSize(double from_val, double to_val, int n_bin,
			double bin_size, double *brks) {
  // We have to make a decision here, the max break should be to_val!
  // no of breaks: ceil(from_val - to_val / bin_size)
  // We have to assume that the *brks array has already been sized to the
  // correct length (i.e. ceil((from_val - to_val) / bin_size))
  double current_val;
  current_val = from_val;
  for (int i = 0; i < n_bin; i++) {
    brks[i] = current_val;
    current_val += bin_size;
  }
  brks[n_bin] = to_val;
  return;
}

/*
 * Bin Y on X based on defined breaks.
 * x and breaks are expected to be sorted incrementally.
 * brks is supposed to be an array with length n_bin + 1 defining the start
 * and end values for the bins.
 * ans is an array with length = n_bin which should have been initialized
 * with 0s!
 * x_start_idx and x_end_idx are the start and end indices in array x in
 * which we're looking for values to be within the bins. This allows to
 * bin on a subset of the x/y arrays. We suppose these have been checked
 * BEFORE (i.e. both being positive and x_end_idx <= length(x)).
 */
static void _bin_y_on_x_with_breaks_max(double *x, double *y, double *brks,
					double *ans, int n_bin, int x_start_idx,
					int x_end_idx) {
  int x_current_idx, last_bin_idx;
  double x_current_value;
  last_bin_idx = n_bin - 1;
  x_current_idx = x_start_idx;

  // o Loop through the bins/brks
  for (int i = 0; i < n_bin; i++) {
    // loop through the x values; assumes x sorted increasingly
    while (x_current_idx <= x_end_idx) {
      x_current_value = x[x_current_idx];
      if (x_current_value >= brks[i]) {
	/* OK, now check if the value is smaller the upper border
	 * OR if we're in the last bin, whether the value matches the upper border.
	 */
	if ((x_current_value < brks[i + 1]) | (x_current_value == brks[i + 1] &
					       i == last_bin_idx)) {
	  /* if (x_current_value == brks[i + 1] & i == last_bin_idx) */
	  /*   Rprintf(" We're in the last bin and the value is equal to the upper brks.\n"); */
	  /* Check if the corresponding y value is larger than the one we have and
	   * replace if so.
	   */
	  if (y[x_current_idx] > ans[i])
	    ans[i] = y[x_current_idx];
	} else {
	  /* Break without incrementing the x_current_idx, thus the same value will
	   * be evaluated for the next bin i.
	   */
	  break;
	}
      }
      x_current_idx++;
    }
  }
  return;
}

static void _bin_midPoint(double *brks, double *bin_mids, int n_bin) {
  for (int i = 0; i < n_bin; i++) {
    bin_mids[i] = (brks[i] + brks[i+1]) / 2;
  }
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


