#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include "util.h"

void DescendZero(double *yvals, int *numin, int *istart, 
                 int *ilower, int *iupper) {
    
    int i;
    
    for (i = *istart; i >= 0; i--)
        if (yvals[i] < 0)
            break;
    *ilower = i + 1;
    
    for (i = *istart; i < *numin; i++)
        if (yvals[i] < 0)
            break;
    *iupper = i - 1;
}

void DescendValue(const double *yvals, const int *numin, const int *istart, 
                  const double *yval, int *ilower, int *iupper) {
    
    int i;
    
    for (i = *istart; i >= 0; i--)
        if (yvals[i] < *yval)
            break;
    *ilower = i + 1;
    
    for (i = *istart; i < *numin; i++)
        if (yvals[i] < *yval)
            break;
    *iupper = i - 1;
}

void DescendMin(double *yvals, int *numin, int *istart, 
                int *ilower, int *iupper) {
    
    int i;
    
    for (i = *istart; i > 0; i--)
        if (yvals[i-1] >= yvals[i])
            break;
    *ilower = i;
    
    for (i = *istart; i < *numin-1; i++)
        if (yvals[i+1] >= yvals[i])
            break;
    *iupper = i;
}

void FindEqualGreaterM(const double *in, const int *size, const double *values, 
                       const int *valsize, int *index) {

    int i, idx = 0;
    
    for (i = 0; i < *valsize; i++) {
        while (idx < *size && in[idx] < values[i])
            idx++;
        index[i] = idx;
    }
}

void FindEqualGreater(const double *in, const int *size, const double *target, 
                      int *index) {
    
    int min = 0, max = *size-1, i = max/2;
    
    while (min < max) {
        if (in[i] < *target)
            min = i+1;
        else
            max = i;
        i = (min+max)/2;
    }
    
    *index = i;
}

void FindEqualLess(const double *in, const int *size, const double *target, 
                   int *index) {

    int min = 0, max = *size-1, i = max/2;
    
    while (min < max) {
        if (in[i] > *target)
            max = i-1;
        else
            min = i;
        i = (int) ceil((min+max)/(float)2);
    }
    
    *index = i;
}

void ColMax(const double *in, const int *n, const int *dn, double *out) {

    int i, j;
    
    for (i = 0; i < *dn; i++) {
        out[i] = in[*n*i];
        for (j = 1; j < *n; j++)
            if (in[*n*i + j] > out[i])
                out[i] = in[*n*i + j];
    }
}

void RowMax(const double *in, const int *dn, const int *p, double *out) {

    int i, j;
    
    for (i = 0; i < *dn; i++) {
        out[i] = in[i];
        for (j = 1; j < *p; j++)
            if (in[i + *dn*j] > out[i])
                out[i] = in[i + *dn*j];
    }
}

SEXP DoubleMatrix(SEXP nrow, SEXP ncol) {
    
    SEXP matrix, dim;
    int  nrowint, ncolint;
    
    nrowint = INTEGER_POINTER(nrow)[0];
    ncolint = INTEGER_POINTER(ncol)[0];
    
    PROTECT(matrix = NEW_NUMERIC(nrowint*ncolint));
    PROTECT(dim = NEW_INTEGER(2));
    INTEGER_POINTER(dim)[0] = nrowint;
    INTEGER_POINTER(dim)[1] = ncolint;
    SET_DIM(matrix, dim);
    
    UNPROTECT(2);
    
    return matrix;
}

SEXP LogicalMatrix(SEXP nrow, SEXP ncol) {
    
    SEXP matrix, dim;
    int  nrowint, ncolint;
    
    nrowint = INTEGER_POINTER(nrow)[0];
    ncolint = INTEGER_POINTER(ncol)[0];
    
    PROTECT(matrix = NEW_LOGICAL(nrowint*ncolint));
    PROTECT(dim = NEW_INTEGER(2));
    INTEGER_POINTER(dim)[0] = nrowint;
    INTEGER_POINTER(dim)[1] = ncolint;
    SET_DIM(matrix, dim);
    
    UNPROTECT(2);
    
    return matrix;
}
