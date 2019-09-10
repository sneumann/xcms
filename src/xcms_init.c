#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void ColMax(void *, void *, void *, void *);
extern void continuousPtsAboveThreshold(void *, void *, void *, void *, void *, void *);
extern void continuousPtsAboveThresholdIdx(void *, void *, void *, void *, void *, void *);
extern void DescendMin(void *, void *, void *, void *, void *);
extern void DescendValue(void *, void *, void *, void *, void *, void *);
extern void DescendZero(void *, void *, void *, void *, void *);
extern void FindEqualGreater(void *, void *, void *, void *);
extern void FindEqualGreaterM(void *, void *, void *, void *, void *);
extern void FindEqualGreaterUnsorted(void *, void *, void *, void *);
extern void FindEqualLess(void *, void *, void *, void *);
extern void MedianFilter(void *, void *, void *, void *, void *, void *);
extern void NetCDFClose(void *, void *);
extern void NetCDFMSPoints(void *, void *, void *, void *, void *, void *, void *);
extern void NetCDFOpen(void *, void *, void *);
extern void NetCDFStrError(void *, void *, void *);
extern void NetCDFVarDouble(void *, void *, void *, void *);
extern void NetCDFVarID(void *, void *, void *, void *);
extern void NetCDFVarInt(void *, void *, void *, void *);
extern void NetCDFVarLen(void *, void *, void *, void *);
extern void ProfBin(void *, void *, void *, void *, void *, void *, void *);
extern void ProfBinLin(void *, void *, void *, void *, void *, void *, void *);
extern void ProfBinLinBase(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ProfBinLinBaseM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ProfBinLinM(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ProfBinM(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ProfIntLin(void *, void *, void *, void *, void *, void *, void *);
extern void ProfIntLinM(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ProfMaxIdx(void *, void *, void *, void *, void *, void *, void *);
extern void ProfMaxIdxM(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void R_mzClust_hclust(void *, void *, void *, void *, void *, void *);
extern void RectUnique(void *, void *, void *, void *, void *, void *, void *);
extern void RowMax(void *, void *, void *, void *);
extern void WhichColMax(void *, void *, void *, void *);
extern void WhichRowMax(void *, void *, void *, void *);

/* .Call calls */
extern SEXP binYonX(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP binYonX_multi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP breaks_on_binSize(SEXP, SEXP, SEXP);
extern SEXP breaks_on_nBins(SEXP, SEXP, SEXP, SEXP);
extern SEXP DoubleMatrix(SEXP, SEXP);
extern SEXP fastMatch(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP findmzROI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getEIC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getMZ(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impute_with_linear_interpolation(SEXP, SEXP);
extern SEXP impute_with_linear_interpolation_base(SEXP, SEXP, SEXP);
extern SEXP IntegerMatrix(SEXP, SEXP);
extern SEXP LogicalMatrix(SEXP, SEXP);
extern SEXP massifquant(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_set_from_xcms(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP test_integer(SEXP);
extern SEXP test_real(SEXP);

static const R_CMethodDef CEntries[] = {
    {"ColMax",                         (DL_FUNC) &ColMax,                          4},
    {"continuousPtsAboveThreshold",    (DL_FUNC) &continuousPtsAboveThreshold,     6},
    {"continuousPtsAboveThresholdIdx", (DL_FUNC) &continuousPtsAboveThresholdIdx,  6},
    {"DescendMin",                     (DL_FUNC) &DescendMin,                      5},
    {"DescendValue",                   (DL_FUNC) &DescendValue,                    6},
    {"DescendZero",                    (DL_FUNC) &DescendZero,                     5},
    {"FindEqualGreater",               (DL_FUNC) &FindEqualGreater,                4},
    {"FindEqualGreaterM",              (DL_FUNC) &FindEqualGreaterM,               5},
    {"FindEqualGreaterUnsorted",       (DL_FUNC) &FindEqualGreaterUnsorted,        4},
    {"FindEqualLess",                  (DL_FUNC) &FindEqualLess,                   4},
    {"MedianFilter",                   (DL_FUNC) &MedianFilter,                    6},
    {"NetCDFClose",                    (DL_FUNC) &NetCDFClose,                     2},
    {"NetCDFMSPoints",                 (DL_FUNC) &NetCDFMSPoints,                  7},
    {"NetCDFOpen",                     (DL_FUNC) &NetCDFOpen,                      3},
    {"NetCDFStrError",                 (DL_FUNC) &NetCDFStrError,                  3},
    {"NetCDFVarDouble",                (DL_FUNC) &NetCDFVarDouble,                 4},
    {"NetCDFVarID",                    (DL_FUNC) &NetCDFVarID,                     4},
    {"NetCDFVarInt",                   (DL_FUNC) &NetCDFVarInt,                    4},
    {"NetCDFVarLen",                   (DL_FUNC) &NetCDFVarLen,                    4},
    {"ProfBin",                        (DL_FUNC) &ProfBin,                         7},
    {"ProfBinLin",                     (DL_FUNC) &ProfBinLin,                      7},
    {"ProfBinLinBase",                 (DL_FUNC) &ProfBinLinBase,                  9},
    {"ProfBinLinBaseM",                (DL_FUNC) &ProfBinLinBaseM,                11},
    {"ProfBinLinM",                    (DL_FUNC) &ProfBinLinM,                     9},
    {"ProfBinM",                       (DL_FUNC) &ProfBinM,                        9},
    {"ProfIntLin",                     (DL_FUNC) &ProfIntLin,                      7},
    {"ProfIntLinM",                    (DL_FUNC) &ProfIntLinM,                     9},
    {"ProfMaxIdx",                     (DL_FUNC) &ProfMaxIdx,                      7},
    {"ProfMaxIdxM",                    (DL_FUNC) &ProfMaxIdxM,                     9},
    {"R_mzClust_hclust",               (DL_FUNC) &R_mzClust_hclust,                6},
    {"RectUnique",                     (DL_FUNC) &RectUnique,                      7},
    {"RowMax",                         (DL_FUNC) &RowMax,                          4},
    {"WhichColMax",                    (DL_FUNC) &WhichColMax,                     4},
    {"WhichRowMax",                    (DL_FUNC) &WhichRowMax,                     4},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"binYonX",                               (DL_FUNC) &binYonX,                               14},
    {"binYonX_multi",                         (DL_FUNC) &binYonX_multi,                         14},
    {"breaks_on_binSize",                     (DL_FUNC) &breaks_on_binSize,                      3},
    {"breaks_on_nBins",                       (DL_FUNC) &breaks_on_nBins,                        4},
    {"DoubleMatrix",                          (DL_FUNC) &DoubleMatrix,                           2},
    {"fastMatch",                             (DL_FUNC) &fastMatch,                              6},
    {"findmzROI",                             (DL_FUNC) &findmzROI,                             10},
    {"getEIC",                                (DL_FUNC) &getEIC,                                 6},
    {"getMZ",                                 (DL_FUNC) &getMZ,                                  6},
    {"impute_with_linear_interpolation",      (DL_FUNC) &impute_with_linear_interpolation,       2},
    {"impute_with_linear_interpolation_base", (DL_FUNC) &impute_with_linear_interpolation_base,  3},
    {"IntegerMatrix",                         (DL_FUNC) &IntegerMatrix,                          2},
    {"LogicalMatrix",                         (DL_FUNC) &LogicalMatrix,                          2},
    {"massifquant",                           (DL_FUNC) &massifquant,                           14},
    {"R_set_from_xcms",                       (DL_FUNC) &R_set_from_xcms,                       18},
    {"test_integer",                          (DL_FUNC) &test_integer,                           1},
    {"test_real",                             (DL_FUNC) &test_real,                              1},
    {NULL, NULL, 0}
};

void R_init_xcms(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
