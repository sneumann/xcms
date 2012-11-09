//DataKeeper.cpp

#include "OpOverload.h"
#include "DataKeeper.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
//#include <fstream>
//#include <iostream>
#include <list>
#include <algorithm>

using namespace std;


DataKeeper::DataKeeper(SEXP mz, SEXP inten, SEXP scanindex, SEXP ls, SEXP scantime) {
 pmz = REAL(mz);
 pinten = REAL(inten);
 pscanindex = INTEGER(scanindex);
 nmz = GET_LENGTH(mz);
 lastScan = INTEGER(ls)[0];
 num_scans = lastScan;
 pscantime = REAL(scantime);
}

DataKeeper::DataKeeper(const char* dotplms1) {
    //something like "window1_5200_4900_764_795.plms1",
    strcpy(filename, dotplms1);
}

DataKeeper::~DataKeeper() {
}

/*void DataKeeper::assign_values(float64* data, uint32_t data_len, vector<double> & vec, int vec_len) {

    uint32_t i;
    //vec.resize(data_len + vec_len);
    for (i = 0; i < data_len; ++i) {
        vec.push_back(data[i]);
    }
}*/

std::vector<double> DataKeeper::getMZScan(int s) {

    int start_pos = scan_idx.at(s);
    int stop_pos = scan_idx.at(s + 1);
   // //cout << "Start" << start_pos << endl;
   // //cout << "Stop" << stop_pos << endl;
    std::vector<double> the_scan(stop_pos - start_pos, 0);
    int  copyvec_idx = 0;
    for (int i = start_pos; i < stop_pos; ++i) {
        the_scan[copyvec_idx] = mz.at(i);
//        //cout << "MZ VAL" << the_scan.at(i) << endl;

        copyvec_idx++;
    }
    return the_scan;
}

void DataKeeper::getScanMQ(int s, std::vector<double> & mzScan, std::vector<double> & intenScan) {

    mzScan.clear();
    intenScan.clear();

    int start_pos = scan_idx.at(s);
    int stop_pos = scan_idx.at(s + 1);

    mzScan = std::vector<double>(stop_pos - start_pos, 0);
    intenScan = std::vector<double>(stop_pos - start_pos, 0);

    int  copyvec_idx = 0;
    for (int i = start_pos; i < stop_pos; ++i) {
        mzScan[copyvec_idx] = mz.at(i);
        intenScan[copyvec_idx] = intensity.at(i);
        copyvec_idx++;
    }
}

void  DataKeeper::getScanXcms(int scan, int nmz, int lastScan, std::vector<double> & mzScan, std::vector<double> & intenScan) {

    //pass in as reference and changes scan to scan
    mzScan.clear();
    intenScan.clear();

    int idx, idx1, idx2, i=0, N=0;
    idx1 = pscanindex[scan - 1] + 1;

    if (scan == lastScan) {
        idx2 = nmz - 1;
    }
    else {
        idx2 = pscanindex[scan];
    }

   // cout << "idx1: " << idx1 << endl;
    N = idx2 - idx1 + 1;
//    scanBuf scbuf;
    if (N > 0) {
        mzScan = std::vector<double>(N);
        intenScan = std::vector<double>(N);
//        cout << "MZ for scan  "<< scan << endl;
        for (idx = idx1; idx <= idx2; idx++) {
//            cout << pmz[idx - 1] << " ";
            mzScan[i] = pmz[idx - 1];
            intenScan[i] = sqrt(pinten[idx - 1]);
            i++;
        }
        //cout << endl;
    }
}

double DataKeeper::getScanTime(int s) {
    return pscantime[s];
}

std::vector<double> DataKeeper::getIScan(int s) {

    int start_pos = scan_idx.at(s);
    int stop_pos = scan_idx.at(s + 1);

    std::vector<double> the_scan(stop_pos - start_pos, 0);
    int  copyvec_idx = 0;
    for (int i = start_pos; i < stop_pos; ++i) {
        the_scan[copyvec_idx] = intensity.at(i);
        copyvec_idx++;
    }
    return the_scan;
}

std::vector<double> DataKeeper::privGetMZScan(int s) {

    int start_pos = scan_idx.at(s);
    int stop_pos = scan_idx.at(s + 1);
   // //cout << "Start" << start_pos << endl;
   // //cout << "Stop" << stop_pos << endl;
    std::vector<double> the_scan(stop_pos - start_pos, 0);
    int  copyvec_idx = 0;
    for (int i = start_pos; i < stop_pos; ++i) {
        the_scan[copyvec_idx] = mz.at(i);
//        //cout << "MZ VAL" << the_scan.at(i) << endl;
        copyvec_idx++;
    }
    return the_scan;
}
std::vector<double> DataKeeper::privGetIScan(int s) {

    int start_pos = scan_idx.at(s);
    int stop_pos = scan_idx.at(s + 1);

    std::vector<double> the_scan(stop_pos - start_pos, 0);
    int  copyvec_idx = 0;
    for (int i = start_pos; i < stop_pos; ++i) {
        the_scan[copyvec_idx] = intensity.at(i);
        copyvec_idx++;
    }
    return the_scan;
}

void DataKeeper::privGetScanXcms(int scan, std::vector<double> & mzScan,
        std::vector<double> & intenScan) {

    //pass in as reference and changes scan to scan
    mzScan.clear();
    intenScan.clear();

    int idx, idx1, idx2, i=0, N=0;
    idx1 = pscanindex[scan - 1] + 1;

    if (scan == lastScan) {
        idx2 = nmz - 1;
    }
    else {
        idx2 = pscanindex[scan];
    }

   // cout << "idx1: " << idx1 << endl;
    N = idx2 - idx1 + 1;
//    scanBuf scbuf;
    if (N > 0) {
        mzScan = std::vector<double>(N);
        intenScan = std::vector<double>(N);
        //Rprintf("extract Scan: %d\n", scan);
//        cout << "Inside the method: " << endl;
        for (idx = idx1; idx <= idx2; idx++) {
//            cout << pmz[idx - 1] << " ";
            mzScan[i] = pmz[idx - 1];
            intenScan[i] = pinten[idx - 1];
            i++;
        }
//        cout << endl;
    }
}

uint32_t DataKeeper::getTotalScanNumbers() {
    return num_scans;
}

int DataKeeper::getTotalCentroidCount() {
    return nmz;
}
double DataKeeper::getInitMZS2() {
    return initMZS2;
}
double DataKeeper::getInitIS2() {
    return  initIS2;
}
double DataKeeper::getInitIS() {
    return initIS;
}

void DataKeeper::ghostScanR() {

    //Rprintf("Inside ghostScanR()\n");
    //find max intensity index
    double apexVal = *max_element(pinten, pinten + nmz);
    //Rprintf("apexVal is %f\n", apexVal);
    initIS = sqrt(apexVal);
     //cout << "R Apex Val: " << apexVal << endl;
    int apexIdx = int(max_element(pinten, pinten + nmz)
                   - pinten);
    double mzApex = pmz[apexIdx];
    //Rprintf("mzApex: %f\n", mzApex);
    //cout << "R apex mz is  " << mz.at(apexIdx) << endl;
    int low = int(lower_bound(pscanindex, pscanindex + lastScan, apexIdx)
                   - pscanindex);
    int up  = int(upper_bound(pscanindex, pscanindex + lastScan, apexIdx)
                   - pscanindex);


    //Rprintf("low: %f\n", low);
    //Rprintf("up: %f\n", up);
    //cout << "R with STL algo, low" << low << endl;
    //cout << "R with STL algo, up" << up << endl;

    int centerScanStartIdx;
    if (low == up) {
        centerScanStartIdx = low - 1;
    }
    else {
        centerScanStartIdx = low;
    }

    /*eg 0, 1, 2, ..., 7 if 3 were the middle idx*/
    std::list<int> range;
    int i;
    for (i = 3; i > 0; i--) {
        range.push_back(centerScanStartIdx - i);

        //cout << "Range is" << centerScanStartIdx - i << endl;
    }
    range.push_back(centerScanStartIdx);
    for (i = 1; i <  4; i++) {
        range.push_back(centerScanStartIdx + i);

                //cout << "Range is" << centerScanStartIdx + i << endl;
    }

    std::list<double> iMaxFeat;
    std::list<double> mMaxFeat;
    std::vector<double> mGhost;
    std::vector<double> iGhost;
    /*Now find the surrounding values*/
    std::list<int>::iterator it_r;
    for (it_r = range.begin(); it_r != range.end(); ++it_r) {

        privGetScanXcms(*it_r + 1, mGhost, iGhost);

        //mGhost = privGetMZScan(*it_r);
        //iGhost = privGetIScan(*it_r);

        double left = mzApex - 0.015;
        double right = mzApex + 0.015;
        std::vector<int> leftLogic = mGhost >= left;
        std::vector<int> rightLogic = mGhost<= right;
        std::vector<int> twos = leftLogic + rightLogic;
        std::vector<int> neighborIdx = twos == 2;

     //cout << "My First Neighbor  is " << mGhost.at(neighborIdx.at(0)) << endl;

        std::vector<double> iGhostSub;
        iGhostSub = copySubIdx(iGhost, neighborIdx);
       //cout << "iSub" << endl;
//        printVec(iGhostSub);

        if (neighborIdx.size() > 0) {
            std::vector<double>::iterator it_imax;
            it_imax = max_element(iGhostSub.begin(), iGhostSub.end());
            int imaxIdx = int(it_imax -iGhostSub.begin());

            //cout << "imaxIdx is "  << imaxIdx << endl;
            //cout << "neignboridx is  "  <<  neighborIdx.at(0) << endl;
            //cout << "val of mz is " << mz.at(neighborIdx.at(imaxIdx)) << endl;
             //cout << "val of int is " << *it_imax << endl;
            iMaxFeat.push_back(*it_imax);
            mMaxFeat.push_back(mGhost.at(neighborIdx.at(imaxIdx)));
        }
    }
   //cout << "Range is " << endl;
    //printList(range);
   //cout << "Final Picks for MZ is " << endl;
    //printList(mMaxFeat);
   //cout << " Final Picks for Intensity is " << endl;
    //printList(iMaxFeat);
    initMZS2 = computeAnySampVar(mMaxFeat);
    initIS2  = computeAnySampVar(iMaxFeat);

  // printVec(intensity);
//intensity = transformIntensity(intensity); //take the sqrt

  // transformIntensityR(); //take the sqrt
}

void DataKeeper::ghostScan() {
    //find max intensity index
    double apexVal = *max_element(intensity.begin(), intensity.end());
    initIS = sqrt(apexVal);
     //cout << "STL Apex Val: " << apexVal << endl;

    int apexIdx = int(max_element(intensity.begin(), intensity.end())
                   - intensity.begin());
    double mzApex = mz.at(apexIdx);
//cout << "the apex mz is  " << mz.at(apexIdx) << endl;
    int low = int(lower_bound(scan_idx.begin(), scan_idx.end(), apexIdx)
                   - scan_idx.begin());
    int up  = int(upper_bound(scan_idx.begin(), scan_idx.end(), apexIdx)
                   - scan_idx.begin());

    //cout << "STL with STL algo, low" << low << endl;
    //cout << "STL with STL algo, up" << up << endl;


    int centerScanStartIdx;
    if (low == up) {
        centerScanStartIdx = low - 1;
    }
    else {
        centerScanStartIdx = low;
    }

    /*eg 0, 1, 2, ..., 7 if 3 were the middle idx*/
    std::list<int> range;
    int i;
    for (i = 3; i > 0; i--) {
        range.push_back(centerScanStartIdx - i);

        //cout << "Range is" << centerScanStartIdx - i << endl;
    }
    range.push_back(centerScanStartIdx);
    for (i = 1; i <  4; i++) {
        range.push_back(centerScanStartIdx + i);

                //cout << "Range is" << centerScanStartIdx + i << endl;
    }

    std::list<double> iMaxFeat;
    std::list<double> mMaxFeat;
   /*Now find the surrounding values*/
    std::list<int>::iterator it_r;
    for (it_r = range.begin(); it_r != range.end(); ++it_r) {

        std::vector<double> mGhost = privGetMZScan(*it_r);
        //cout << "mGhost" << endl;
        //printVec(mGhost);

        std::vector<double> iGhost = privGetIScan(*it_r);
        //cout << "iGhost" << endl;
       // printVec(iGhost);

        double left = mzApex - 0.015;
        double right = mzApex + 0.015;
        std::vector<int> leftLogic = mGhost >= left;
        std::vector<int> rightLogic = mGhost<= right;

        //cout << "Left is " << left << endl;
        std::vector<int> twos = leftLogic + rightLogic;
        std::vector<int> neighborIdx = twos == 2;

     //cout << "My First Neighbor  is " << mGhost.at(neighborIdx.at(0)) << endl;

        std::vector<double> iGhostSub;
        iGhostSub = copySubIdx(iGhost, neighborIdx);
       //cout << "iSub" << endl;
        //printVec(iGhostSub);

        if (neighborIdx.size() > 0) {
            std::vector<double>::iterator it_imax;
            it_imax = max_element(iGhostSub.begin(), iGhostSub.end());
            int imaxIdx = int(it_imax -iGhostSub.begin());

            //cout << "imaxIdx is "  << imaxIdx << endl;
            //cout << "neignboridx is  "  <<  neighborIdx.at(0) << endl;
            //cout << "val of mz is " << mz.at(neighborIdx.at(imaxIdx)) << endl;
             //cout << "val of int is " << *it_imax << endl;
            iMaxFeat.push_back(*it_imax);
            mMaxFeat.push_back(mGhost.at(neighborIdx.at(imaxIdx)));
        }
    }
   //cout << "Range is " << endl;
   // printList(range);
   //cout << "Final Picks for MZ is " << endl;
   // printList(mMaxFeat);
   //cout << " Final Picks for Intensity is " << endl;
   // printList(iMaxFeat);
    initMZS2 = computeAnySampVar(mMaxFeat);
    initIS2  = computeAnySampVar(iMaxFeat);

  // printVec(intensity);

    intensity = transformIntensity(intensity); //take the sqrt

 // printVec(intensity);
}
/*
double DataKeeper::findMaxIdx() {

    double maxIdx = 0;
    double max = intensity.at(0);

    for (unsigned int i = 0; i < intensity.size(); i++) {
           if (intensity.at(i) >
    }


}
*/

/*std::vector<double> DataKeeper::getMzTrFootprints(int s) {

}*/


void DataKeeper::printVec(const std::vector<double> & myvec) {
    for (size_t i = 0; i < myvec.size(); i++) {
        //cout << myvec.at(i) << "\t";
    }
    //cout << "\n";
}

void DataKeeper::printList(const std::list<int> & mylist) {
    std::list<int>::const_iterator it;
    for (it = mylist.begin(); it != mylist.end(); ++it) {
        //cout << *it  << "\t";
    }
    //cout << "\n";
}

void DataKeeper::printList(const std::list<double> & mylist) {
    std::list<double>::const_iterator it;
    for (it = mylist.begin(); it != mylist.end(); ++it) {
        //cout << *it  << "\t";
    }
   //cout << "\n";
}

std::vector<double> DataKeeper::transformIntensity(std::vector<double> & A) {
    unsigned int i;
    for (i = 0; i < A.size(); ++i) {
      A[i] = sqrt(A.at(i));
    }
    return A;
}

void  DataKeeper::transformIntensityR() {
    int i;
    for (i = 0; i < nmz; ++i) {
      pinten[i] = sqrt(pinten[i]);
    }
}
