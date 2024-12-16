// Massifquant

//STDLIB
#include <math.h>
#include <vector>
#include <list>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>

//MASSIFQUANT
#include "OpOverload.h"
#include "Tracker.h"
#include "TrMgr.h"
#include "DataKeeper.h"
#include "SegProc.h"

// R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

const int N_NAMES = 7;
using namespace std;

extern "C" SEXP massifquant(SEXP mz, SEXP intensity, SEXP scanindex,
        SEXP scantime, SEXP mzrange, SEXP scanrange, SEXP lastscan,
        SEXP minIntensity, SEXP minCentroids, SEXP consecMissedLim,
        SEXP ppm, SEXP criticalVal, SEXP segs, SEXP scanBack) {

    //the return data structure and its elemental components
    //jo SEXP peaklist,entrylist,list_names,vmz,vmzmin,vmzmax,vstcenter,vscmin,vscmax,vintensity,vintenmax, vlength;
    SEXP peaklist,entrylist,list_names,vmz,vmzmin,vmzmax,vscmin,vscmax,vintensity,vlength;
    int scanrangeTo, scanrangeFrom;
    int firstScan = 1;

    scanrangeFrom = INTEGER(scanrange)[0];
    scanrangeTo = INTEGER(scanrange)[1];

    //store data
    DataKeeper dkeep(mz, intensity, scanindex, lastscan, scantime);
    dkeep.ghostScanR();
    std::vector<double> mzScan;
    std::vector<double> intenScan;

    int totalScanNums = dkeep.getTotalScanNumbers();
    int nmz = dkeep.getTotalCentroidCount();
    double iq =  dkeep.getInitIS2();
    double mzq = dkeep.getInitMZS2();
    double mzr =  sqrt(mzq);
    double ir = dkeep.getInitIS();
    double * pscantime = REAL(scantime);

    //model breaks down otherwise
    if (mzq == 0) {
        mzq = 1e-6;
        mzr = sqrt(mzq);
    }

    if ((scanrangeFrom <  firstScan) || (scanrangeFrom > totalScanNums) || (scanrangeTo < firstScan) || (scanrangeTo > totalScanNums))
        Rf_error("Error in scanrange \n");

    //show the progress please
    Rprintf("\n Detecting Kalman ROI's ... \n percent finished: ");

    //initialize tracker manager
    TrMgr busybody(scanrangeTo, sqrt(REAL(minIntensity)[0]),
            INTEGER(minCentroids)[0], REAL(consecMissedLim)[0],
            REAL(ppm)[0], REAL(criticalVal)[0], INTEGER(scanBack)[0]);
    dkeep.getScanXcms(scanrangeTo, nmz, totalScanNums, mzScan, intenScan);
    busybody.setDataScan(mzScan, intenScan);
    busybody.initTrackers(iq, mzq, ir, mzr, scanrangeTo);
    //begin feature finding
    //Rprintf("scanrangeTo: %d\n", scanrangeTo);
    double progCount = 0;
    //jo double maxScanNums = double(scanrangeTo);
    double progThresh = 10;
    for (int k = scanrangeTo - 1; k >= scanrangeFrom; k--) {

        //progress
        double perc  = (progCount/scanrangeTo) * 100;
        //Rprintf("perc: %f\t", perc);
        if (perc > progThresh) {
            Rprintf(" %d  ", int(perc));
            progThresh += 10;
        }

        busybody.setCurrScanIdx(k);
        dkeep.getScanXcms(k, nmz, totalScanNums, mzScan, intenScan);
        busybody.predictScan(mzScan, intenScan);
        busybody.competeAct();
        busybody.manageMissed();
        busybody.manageTracked();
        busybody.initTrackers(iq, mzq, ir, mzr, k);
        progCount++;
    }

    busybody.removeOvertimers();

    //option to include segmentation correction
    //string segsStr(segs);
    //if (segsStr.compare("segOn") == 0) {
    if (INTEGER(segs)[0] == 1) {
        SegProc sproc(busybody.getPicCounts());
        sproc.groupSegments(busybody);
        sproc.collapseSubsets();
        //sproc.segsToFile(busybody);
        sproc.solderSegs(busybody);
    }

    Rprintf(" %d\n", 100);

    const char *names[N_NAMES] = {"mz", "mzmin", "mzmax", "scmin", "scmax", "length", "intensity"};
    PROTECT(list_names = Rf_allocVector(STRSXP, N_NAMES));
    for(int j = 0; j < N_NAMES; j++)
        SET_STRING_ELT(list_names, j,  Rf_mkChar(names[j]));

    PROTECT(peaklist = Rf_allocVector(VECSXP, busybody.getPicCounts()));
    for (int i=0;i<busybody.getPicCounts();i++) {

        std::vector<double> featInfo = busybody.iterOverFeatures(i, pscantime);
	//jo int scanLength = int(featInfo.at(5) - featInfo.at(4) + 1);

        PROTECT(entrylist = Rf_allocVector(VECSXP, N_NAMES));

        //allow for new vars declared to be passed out
        PROTECT(vmz = NEW_NUMERIC(1));
        PROTECT(vmzmin = NEW_NUMERIC(1));
        PROTECT(vmzmax = NEW_NUMERIC(1));

        PROTECT(vscmin = NEW_INTEGER(1));
        PROTECT(vscmax = NEW_INTEGER(1));
        PROTECT(vlength = NEW_INTEGER(1));
        PROTECT(vintensity = NEW_INTEGER(1));

        //assign a copy
        NUMERIC_POINTER(vmz)[0]  = featInfo.at(0);
        NUMERIC_POINTER(vmzmin)[0] = featInfo.at(1);
        NUMERIC_POINTER(vmzmax)[0] = featInfo.at(2);

        INTEGER_POINTER(vscmin)[0] = int(featInfo.at(4));
        INTEGER_POINTER(vscmax)[0] = int(featInfo.at(5));
        INTEGER_POINTER(vlength)[0] = int(featInfo.at(3));
        INTEGER_POINTER(vintensity)[0] = int(featInfo.at(6));

        //enter into nested list
        SET_VECTOR_ELT(entrylist, 0, vmz);
        SET_VECTOR_ELT(entrylist, 1, vmzmin);
        SET_VECTOR_ELT(entrylist, 2, vmzmax);
        SET_VECTOR_ELT(entrylist, 3, vscmin);
        SET_VECTOR_ELT(entrylist, 4, vscmax);
        SET_VECTOR_ELT(entrylist, 5, vlength);
        SET_VECTOR_ELT(entrylist, 6, vintensity);

        Rf_setAttrib(entrylist, R_NamesSymbol, list_names); //attaching the vector names
        SET_VECTOR_ELT(peaklist, i, entrylist);
        UNPROTECT(N_NAMES + 1); //entrylist + values
    }

    //busybody.writePICsToFile();
    //Rprintf("Number detected: %d\n", busybody.getPicCounts());

    UNPROTECT(2);//peaklist, list_names

    return (peaklist);
}
