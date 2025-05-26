#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

// Forward declarations for functions from mzROI_xcms_extracted.c
SEXP findmzROI(SEXP mz, SEXP intensity, SEXP scanindex, SEXP mzrange, SEXP scanrange, SEXP lastscan, SEXP dev, SEXP minEntries, SEXP prefilter, SEXP noise);
SEXP getEIC(SEXP mz, SEXP intensity, SEXP scanindex, SEXP mzrange, SEXP scanrange, SEXP lastscan);
SEXP getWeightedMZ(SEXP mz, SEXP intensity, SEXP scanindex, SEXP mzrange, SEXP scanrange, SEXP lastscan);

// Forward declarations for functions from util_xcms_extracted.c
// Note: RectUnique is void but called via .Call, so its SEXP version for registration might differ
// or it might be intended for .C call. Assuming .Call for now as per R's modern practice.
// If it's a .C interface, the signature in CallEntries would be different.
// For the provided CallEntries, its signature implies it's being treated as if it returns SEXP,
// which is unusual for a void C function via .Call.
// However, sticking to the provided CallEntries structure.
// The same applies to DescendMin and continuousPtsAboveThreshold if they were void.
// The C functions themselves are void, but the DL_FUNC cast in CallEntries is generic.
// The number of arguments in CallEntries is what R checks.

void RectUnique(const double *m, const int *order, const int *nrow, const int *ncol, const double *xdiff, const double *ydiff, int *keep);
void DescendMin(double *yvals, int *numin, int *istart, int *ilower, int *iupper);
void continuousPtsAboveThreshold(double *x, int *istart, int *numin, double *threshold, int *num, int *n);

static const R_CallMethodDef CallEntries[] = {
    {"findmzROI",                   (DL_FUNC) &findmzROI,                   10},
    {"getEIC",                      (DL_FUNC) &getEIC,                       6},
    {"getWeightedMZ",               (DL_FUNC) &getWeightedMZ,                6},
    {"RectUnique",                  (DL_FUNC) &RectUnique,                   7},
    {"DescendMin_C",                (DL_FUNC) &DescendMin,                   5},
    {"continuousPtsAboveThreshold_C", (DL_FUNC) &continuousPtsAboveThreshold, 6},
    {NULL, NULL, 0}
};

void R_init_xcms_extracted_lib(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    // Optional: If you want to make symbols available for .C calls specifically
    // R_forceSymbols(dll, TRUE); 
}
