// From binners.c
void binYonX_max(double *x, double *y, int *numin, double *xstart,
		 double *xend, int *numout, double *out);

/*
 * Calculate breaks from_val to to_val based on the number of bins.
 * Results are directly stored into array brks which size has to be
 * correctly defined in the calling function.
 */
void _breaks_on_nBins(double from_val, double to_val, int n_bin,
		      double *brks);

/*
 * Calculate breaks from_val to to_val providing the bin size.
 * Results are directly stored into the array brks which size has to
 * be correctly defined in the calling function. Note that the last
 * element in the array will always be to_val (thus the size of the
 * brks array has to be ceil((to_val - from_val) / bin_size)).
 */
void _breaks_on_binSize(double from_val, double to_val, int n_bin,
			double bin_size, double *brks);
