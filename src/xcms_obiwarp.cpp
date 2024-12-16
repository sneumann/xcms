// STDLIB:
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include "string.h"

// Obiwarp
#include "obiwarp/vec.h"
#include "obiwarp/mat.h"
#include "obiwarp/lmat.h"
#include "obiwarp/xcms_dynprog.h"

// R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>


/********************************************/
//char * VERSION = "0.9.2";
/********************************************/

#define DEBUG (0)

extern "C" SEXP R_set_from_xcms(SEXP valscantime, SEXP scantime, SEXP mzrange, SEXP mz, SEXP intensity,
				SEXP valscantime2, SEXP scantime2, SEXP mzrange2, SEXP mz2, SEXP intensity2,
				SEXP response, SEXP score,
				SEXP gap_init, SEXP gap_extend,
				SEXP factor_diag, SEXP factor_gap,
				SEXP local_alignment, SEXP init_penalty)
{

  // Create two matrices in LMata format

    int pvalscantime, pmzrange;
    int pvalscantime2, pmzrange2;
    double *pscantime, *pmz, *pintensity;
    double *pscantime2, *pmz2, *pintensity2;
    SEXP corrected;

    PROTECT(valscantime = Rf_coerceVector(valscantime, INTSXP));
    mzrange = Rf_coerceVector(mzrange, INTSXP);
    pvalscantime = INTEGER(valscantime)[0];
    pmzrange = INTEGER(mzrange)[0];
    pscantime = REAL(scantime);
    pmz = REAL(mz);
    pintensity = REAL(intensity);

    PROTECT(valscantime2 = Rf_coerceVector(valscantime2, INTSXP));
    mzrange2 = Rf_coerceVector(mzrange2, INTSXP);
    pvalscantime2 = INTEGER(valscantime2)[0];
    pmzrange2 = INTEGER(mzrange2)[0];
    pscantime2 = REAL(scantime2);
    pmz2 = REAL(mz2);
    pintensity2 = REAL(intensity2);

    // ************************************************************
    // * READ IN FILES TO GET MAT
    // ************************************************************
    LMat lmat1;
    LMat lmat2;
    MatF smat;
    DynProg dyn;

    lmat1.set_from_xcms(pvalscantime, pscantime, pmzrange, pmz, pintensity);
    lmat2.set_from_xcms(pvalscantime2, pscantime2, pmzrange2, pmz2, pintensity2);

    // ************************************************************
    // * SCORE THE MATRICES
    // ************************************************************
    if (DEBUG) {
      std::cerr << "Scoring the mats!\n";
    }

    dyn.score(*(lmat1.mat()), *(lmat2.mat()), smat, CHAR(STRING_ELT(score, 0)));

    if (DEBUG) {
      std::cerr << "Checking scoring\n";
    }

    if (!strcmp(CHAR(STRING_ELT(score, 0)),"euc")) {
      smat *= -1; // inverting euclidean
    }


    // ************************************************************
    // * PREPARE GAP PENALTY ARRAY
    // ************************************************************

    MatF time_tester;
    MatF time_tester_trans;
    VecF mpt;
    VecF npt;
    VecF mOut_tm;
    VecF nOut_tm;

    int gp_length = smat.rows() + smat.cols();

    VecF gp_array;
    dyn.linear_less_before(*REAL(gap_extend), *REAL(gap_init), gp_length, gp_array);

    // ************************************************************
    // * DYNAMIC PROGRAM
    // ************************************************************
    int minimize = 0;
    if (DEBUG) {
        std::cerr << "Dynamic Time Warping Score Matrix!\n";
    }
    dyn.find_path(smat, gp_array, minimize,
		  *REAL(factor_diag), *REAL(factor_gap), *INTEGER(AS_INTEGER(local_alignment)), *REAL(init_penalty));

    VecI mOut;
    VecI nOut;
    dyn.warp_map(mOut, nOut, *INTEGER(AS_INTEGER(response)), minimize);

    VecF nOutF;
    VecF mOutF;
    lmat1.tm_axis_vals(mOut, mOutF);
    lmat2.tm_axis_vals(nOut, nOutF);
    lmat2.warp_tm(nOutF, mOutF);

    PROTECT(corrected = Rf_allocVector(REALSXP, Rf_length(scantime2)));
    for(int i=0; i < Rf_length(scantime2);i++){
      REAL(corrected)[i] = lmat2.tm()->back()[i];
    }

    UNPROTECT(3);

    return corrected;

}