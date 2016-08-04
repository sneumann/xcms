#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

// From binners.c
void binYonX_max(double *x, double *y, int *numin, double *xstart,
		 double *xend, int *numout, double *out);

/*
 * Calculate breaks from_val to to_val based on the number of bins.
 * Results are directly stored into array brks which size has to be
 * correctly defined in the calling function.
 * shift_by_half_bin_size: either 0 or 1, if 1 the bin mid-points will be shifted left
 * by half of the bin_size.
 */
void _breaks_on_nBins(double from_val, double to_val, int n_bin,
		      double *brks, int shift_by_half_bin_size);

/*
 * Calculate breaks from_val to to_val providing the bin size.
 * Results are directly stored into the array brks which size has to
 * be correctly defined in the calling function. Note that the last
 * element in the array will always be to_val (thus the size of the
 * brks array has to be ceil((to_val - from_val) / bin_size)).
 */
void _breaks_on_binSize(double from_val, double to_val, int n_bin,
			double bin_size, double *brks);

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
					int x_end_idx);

static void _bin_y_on_x_with_breaks_min(double *x, double *y, double *brks,
					double *ans, int n_bin, int x_start_idx,
					int x_end_idx);

static void _bin_y_on_x_with_breaks_sum(double *x, double *y, double *brks,
					double *ans, int n_bin, int x_start_idx,
					int x_end_idx);

static void _bin_y_on_x_with_breaks_mean(double *x, double *y, double *brks,
					 double *ans, int n_bin, int x_start_idx,
					 int x_end_idx);


/*
 * Simple function to calculate the midpoint of bins, breaks provided.
 */
static void _bin_midPoint(double *brks, double *bin_mids, int n_bin);
