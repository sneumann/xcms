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

const int N_NAMES = 8;
using namespace std;

extern "C" SEXP massifquant(SEXP mz, SEXP intensity, SEXP scanindex, 
        SEXP scantime, SEXP mzrange, SEXP scanrange, SEXP lastscan, 
        SEXP minIntensity, SEXP minCentroids, SEXP consecMissedLim, 
        SEXP ppm, SEXP criticalVal, SEXP segs, SEXP scanBack) { 
 
    //the return data structure and its elemental components
    SEXP featurelist,entrylist,list_names,vmzmean,vmzmin,vmzmax,vstcenter,vstmin,vstmax,varea,vintenmax;
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

    if ((scanrangeFrom <  firstScan) || (scanrangeFrom > totalScanNums) || (scanrangeTo < firstScan) || (scanrangeTo > totalScanNums))
        error("Error in scanrange \n");

    //show the progress please
    Rprintf("\n Detecting chromatographic peaks ... \n percent finished: ");
       
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
    double maxScanNums = double(scanrangeTo);
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
        sproc.segsToFile(busybody);
        sproc.solderSegs(busybody);
    }
    
    Rprintf(" %d\n", 100);
 

/*
    char *names[N_NAMES] = {"mz", "mz min", "mz max", "st min", "st max", "integration", "max intensity"};
    PROTECT(list_names = allocVector(STRSXP, N_NAMES));
    for(int j = 0; j < N_NAMES; j++)
        SET_STRING_ELT(list_names, j,  mkChar(names[j]));
*/

    PROTECT(featurelist = allocVector(VECSXP, busybody.getPicCounts()));
    for (int i=0;i<busybody.getPicCounts();i++) {

        std::vector<double> featInfo = busybody.iterOverFeatures(i, pscantime);
        
        PROTECT(entrylist = allocVector(VECSXP, N_NAMES));

        //allow for new vars declared to be passed out
        PROTECT(vmzmean = NEW_NUMERIC(1));
        PROTECT(vmzmin = NEW_NUMERIC(1));
        PROTECT(vmzmax = NEW_NUMERIC(1));
        PROTECT(vstcenter = NEW_NUMERIC(1));
        PROTECT(vstmin = NEW_NUMERIC(1));
        PROTECT(vstmax = NEW_NUMERIC(1));
        PROTECT(varea = NEW_NUMERIC(1));
        PROTECT(vintenmax = NEW_NUMERIC(1));

        //assign a copy
        NUMERIC_POINTER(vmzmean)[0]  = featInfo.at(0);
        NUMERIC_POINTER(vmzmin)[0] = featInfo.at(1);
        NUMERIC_POINTER(vmzmax)[0] = featInfo.at(2);
        NUMERIC_POINTER(vstcenter)[0] = featInfo.at(3);
        NUMERIC_POINTER(vstmin)[0] = featInfo.at(4);
        NUMERIC_POINTER(vstmax)[0] = featInfo.at(5);
        NUMERIC_POINTER(varea)[0] = featInfo.at(6);
        NUMERIC_POINTER(vintenmax)[0] = featInfo.at(7);

        //enter into nested list
        SET_VECTOR_ELT(entrylist, 0, vmzmean);
        SET_VECTOR_ELT(entrylist, 1, vmzmin);
        SET_VECTOR_ELT(entrylist, 2, vmzmax);
        SET_VECTOR_ELT(entrylist, 3, vstcenter);
        SET_VECTOR_ELT(entrylist, 4, vstmin);
        SET_VECTOR_ELT(entrylist, 5, vstmax);
        SET_VECTOR_ELT(entrylist, 6, varea);
        SET_VECTOR_ELT(entrylist, 7, vintenmax);

  //      setAttrib(entrylist, R_NamesSymbol, list_names); //attaching the vector names
        SET_VECTOR_ELT(featurelist, i, entrylist);
        UNPROTECT(N_NAMES + 1); //entrylist + values
    }

    //busybody.writePICsToFile();
    Rprintf("There were %d features detected\n", busybody.getPicCounts());

    UNPROTECT(1);//featurelist

    return (featurelist);  
}

