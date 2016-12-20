//SegProc.cpp

#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "OpOverload.h"
#include "TrMgr.h"
#include "SegProc.h"

#include <R.h>
#include <Rdefines.h>

using namespace std;
SegProc::SegProc(int otn) : origTrNum(otn) {

    segClusters = std::vector<int> (origTrNum, 0);
    unionIdx.push_back(0);
}

SegProc::~SegProc() {
}

void SegProc::groupSegments(TrMgr & busybody) {

    //cout << "PIC Counts Before SegProc: " << busybody.getPicCounts() << endl;

    std::list<int> candIdx;

    //sort the pic According to retention time starting values
    //(i.e. scan numbers)

    int ppm = busybody.getPpm();

    std::vector<int>::iterator it_i;
    std::vector<int> picIdx = busybody.getPicIdx();
    int i  = -1; //use for segCluster indexing
    for (it_i = picIdx.begin(); it_i != picIdx.end(); ++it_i) {
        ++i;
        candIdx.clear(); //candidates are different for each tracker
        //mz mean check
        double mzTol =  busybody.getTracker(*it_i)->getXbar() *  ppm / 1e6;
        std::vector<int>::iterator it_j;
        for (it_j = picIdx.begin(); it_j != picIdx.end(); ++it_j) {
            if (*it_i == *it_j) { continue; }
            double jmeandiff = fabs(busybody.getTracker(*it_i)->getXbar() -
                    busybody.getTracker(*it_j)->getXbar());
            //cout << " " <<  endl;
            if (jmeandiff < mzTol) {
                candIdx.push_back(*it_j);
            //    cout << "*****New Means Comparison *****" << endl;
            //    cout << "j xbar: " <<  busybody.getTracker(*it_j)->getXbar() << "\t";
             //   cout << "i xbar: " <<  busybody.getTracker(*it_i)->getXbar() << endl;
            }
        }

        //cout << "candIdx.size(): " << candIdx.size() << endl;
        //cout << "segClusters.at(i):  " <<  segClusters.at(i) << endl;
        //cout << "i: " << i << endl;
        if (candIdx.size() == 0 || segClusters.at(i) != 0) { continue; }
        segClusters[i] = 1;

        size_t tmpSegCounts = segIdx.size(); //current segments processed

        //retention time check
        //cout << "**Consider Edge Connection**" <<  endl;
        std::list<int> edges;
        std::list<int>::iterator it;
        for (it = candIdx.begin(); it != candIdx.end(); ++it) {
            //cout << "*it start " <<
          //      busybody.getTracker(*it)->getStartScanIdx() << endl;
            //cout << "*i start " <<
            //    busybody.getTracker(*it_i)->getStartScanIdx() << endl;
            //cout << "*i stop " <<
              //  busybody.getTracker(*it_i)->getStopScanIdx() << endl;

            if (busybody.getTracker(*it)->getStartScanIdx() >
                    busybody.getTracker(*it_i)->getStartScanIdx() &&
                busybody.getTracker(*it)->getStartScanIdx() -
                busybody.getTracker(*it_i)->getStopScanIdx() < MAXGAP) {
                    edges.push_back(*it);
                    //cout << "seed is: " << *it_i << endl;
                    //cout << "edge pushed on: " << *it << endl;
            }
        }
        //no edges between trackers were found
        if (edges.size() == 0) { continue; }

        compareMeans(busybody, *it_i, edges, segIdx.size());

        //update the union information and add the seed to the trk_idx
        if (segIdx.size() > tmpSegCounts) {
            segIdx.push_back(*it_i);
            unionIdx.push_back(segIdx.size());
        }
    }
}

void SegProc::collapseSubsets() {

    //cout << "Before Collapsing, segIdx.size(): " << segIdx.size() << endl;
    //cout << "Before Collapsing, unionIdx.size(): " << unionIdx.size() << endl;
    if (segIdx.size() == 0) {
        //cout << "No trackers with same means" << endl;
        return;
    }

    int segCounts = 0;
    int cleanCounts = 1;
    cleanSegIdx = std::vector<int> (segIdx.size(), 0);
    cleanUnionIdx = std::vector<int> (unionIdx.size(), 0);
    std::vector<int> noDups(segIdx.size(), 1);

    for (size_t i = 0; i < unionIdx.size() - 1; ++i) {

        if (noDups.at(i) == 0) { continue; }

        //jo int combinedUnions = 0;
        std::vector<int> iUIdx = createSequence(unionIdx.at(i), unionIdx.at(i + 1) - 1, 1);
        std::vector<int> iSegs2Union = copySubIdx(segIdx, iUIdx);

        //the union of all the subsets
        //std::vector<int> iSegs2Union = iSegs2Union;

       /*cout << "Create Final Union i: " << endl;
        cout << "iSegs2Union: " << endl;
        for (size_t test1 = 0; test1 < iSegs2Union.size(); test1++) {
            cout << iSegs2Union.at(test1) << " ";
        }
        cout << endl;*/

        for (size_t j = i + 1; j < unionIdx.size() - 1; ++j) {

            if (i == j) { continue; }

            std::vector<int> jUIdx = createSequence(unionIdx.at(j), unionIdx.at(j + 1) - 1, 1);
            std::vector<int> jSIdx = copySubIdx(segIdx, jUIdx);

            /*cout << "jSIdx: " << endl;
            for (size_t test1 = 0; test1 < jSIdx.size(); test1++) {
                cout << jSIdx.at(test1) << " ";
            }
            cout << endl;*/

            int jCounts = jSIdx.size();

            int maxSet = max(iSegs2Union.size(), jSIdx.size());
            std::vector<int> shareSegs(maxSet, -1);
            std::vector<int>::iterator it_s;

            sort(iSegs2Union.begin(), iSegs2Union.end());
            sort(jSIdx.begin(), jSIdx.end());
            it_s = set_intersection(iSegs2Union.begin(), iSegs2Union.end(),
                    jSIdx.begin(), jSIdx.end(),
                    shareSegs.begin());
            int interEleNum =  int(it_s - shareSegs.begin());
            if (interEleNum > 0 && // there was an intersection of elements
                    jCounts != interEleNum) { //there is not perfect containment

                //cout << "...unite them at J: " << j << endl;
                //jo combinedUnions = 1;
                noDups.at(j) = 0;
                std::vector<int> tmpVec(iSegs2Union.size() + jSIdx.size() - 1, -1);

                std::vector<int>::iterator it_u;
                it_u = set_union(iSegs2Union.begin(), iSegs2Union.end(),
                        jSIdx.begin(), jSIdx.end(),
                        tmpVec.begin());
                tmpVec.erase(it_u, tmpVec.end());
                iSegs2Union = tmpVec;

                /*cout << "iSegs2Union: " << endl;
                for (size_t test1 = 0; test1 < iSegs2Union.size(); test1++) {
                    cout << iSegs2Union.at(test1) << " ";
                }
                cout << endl;*/
            }
        }
        //add to clean data structures
        for (size_t k = 0; k < iSegs2Union.size(); ++k) {
            cleanSegIdx[segCounts + k] = iSegs2Union.at(k);
        }
        segCounts += iSegs2Union.size();
        cleanUnionIdx[cleanCounts] = segCounts;
        cleanCounts++;
    }
    if (cleanCounts == 1) { return; }
    cleanSegIdx.assign(cleanSegIdx.begin(), cleanSegIdx.begin() + segCounts);
    cleanUnionIdx.assign(cleanUnionIdx.begin(), cleanUnionIdx.begin() + cleanCounts);

    /*cout << "The final cleanSegIdx size: " << cleanSegIdx.size() << endl;
    cout << "The final cleanUnionIndex size: " << cleanUnionIdx.size() << endl;
    cout << "cleanUnionIndex: " << endl;
    for (size_t test1 = 0; test1 < cleanUnionIdx.size(); test1++) {
        cout << cleanUnionIdx.at(test1) << " ";
    }
    cout << endl;
    cout << "cleanSegIdx: " << endl;
    for (size_t test1 = 0; test1 < cleanSegIdx.size(); test1++) {
        cout << cleanSegIdx.at(test1) << " ";
    }
    cout << endl;*/

    //Defensive programming test
    /*int segS = cleanSegIdx.size();
    cout << "Seg S: " << segS << "\t";
    std::vector<int>::iterator testU = unique(cleanSegIdx.begin(), cleanSegIdx.end());
    int segSU = int(testU - cleanSegIdx.begin());
    cout << "segSU: " << segSU << endl;
    //assert(segS == segSU);*/

}

void SegProc::compareMeans(TrMgr & busybody, const int seed, const std::list<int> edges, const int & segCounts) {

    std::list<int>::const_iterator it;
    for (it = edges.begin(); it != edges.end(); ++it) {
        double varRatio = busybody.getTracker(seed)->getS2() / busybody.getTracker(*it)->getS2();
        if (varRatio < TROBUST1 ||
            varRatio > TROBUST2) {

            ttestEq(busybody.getTracker(seed)->getXbar(),
                    busybody.getTracker(*it)->getXbar(),
                    busybody.getTracker(seed)->getTrLen(),
                    busybody.getTracker(*it)->getTrLen(),
                    busybody.getTracker(seed)->getS2(),
                    busybody.getTracker(*it)->getS2());

        }
        else {
            ttestWelch(busybody.getTracker(seed)->getXbar(),
                    busybody.getTracker(*it)->getXbar(),
                    busybody.getTracker(seed)->getTrLen(),
                    busybody.getTracker(*it)->getTrLen(),
                    busybody.getTracker(seed)->getS2(),
                    busybody.getTracker(*it)->getS2());
        }

        p  = 2*pt(fabs(t), v, 0, 0);
        if (p < ALPHA) { continue; }

        //if it arrived here -> no significant difference between means

        segIdx.push_back(*it);
//        segClusters[*it] = *it; //mark as having been checked (might be
//        faulty);

    }
}

void SegProc::solderSegs(TrMgr & busybody) {


    if (cleanUnionIdx.size() == 0) { return; }

    std::vector<int> erasedPicIdx (cleanSegIdx.size() - cleanUnionIdx.size() + 1, 0);
    //cout << "erasedPicIdx.size(): " << erasedPicIdx.size() << endl;
    int j = 0;
    for (size_t i = 0; i < (cleanUnionIdx.size() - 1); ++i) {

        //cout << "i: " << i << endl;
        std::vector<int> uIdx = createSequence(cleanUnionIdx.at(i),
                cleanUnionIdx.at(i + 1) - 1, 1);
        std::vector<int> sIdx = copySubIdx(cleanSegIdx, uIdx);

        //combine all other segments
        //let the last element be the one to combine all others
        std::vector<int>::iterator it;
        int appIdx = sIdx.back();
        for(it = sIdx.begin(); it != (sIdx.end() - 1); ++it) {
            //cout << "Deleting this pic idx: " << *it << endl;
            //int testOSize = busybody.getTracker(appIdx)->getScanList().size();
            list<int> sl = busybody.getTracker(*it)->getScanList();
            list<int> cl = busybody.getTracker(*it)->getCentroidList();
            list<double> ml = busybody.getTracker(*it)->getMzList();
            list<double> il = busybody.getTracker(*it)->getIntensityList();
	    busybody.getTracker(appIdx)->appendToTracker(sl, cl, ml, il);

            //assert((testOSize + sl.size()) == busybody.getTracker(appIdx)->getScanList().size());
            erasedPicIdx[j] = *it;
            j++;
        }
    }
    //erase the PicElements after the iteration has happened
    //cout << "The memory corruption begins here" << endl;
    busybody.erasePicElements(erasedPicIdx);
    Rprintf("\n The number of ROI'S that collapsed into a larger ROI: %d\n",j);
    //cout << "erasedPicIdx.back: " << erasedPicIdx.back() << endl;
}

void SegProc::segsToFile(TrMgr & busybody) {

    if (cleanUnionIdx.size() == 0) { return; }

    ofstream unionfile_idx("unionfile_idx.txt");

    ofstream feat_idx("seg_feat_idx.txt");
    ofstream scan_idx("seg_scan_idx.txt");
    ofstream cent_idx("seg_cent_idx.txt");

    feat_idx << 0 << endl;;
    int totalCents = 0;
    //iterate over all unions
    for (size_t i = 0; i < (cleanUnionIdx.size() - 1); ++i) {
        unionfile_idx << cleanUnionIdx.at(i)  << endl;

        std::vector<int> uIdx = createSequence(cleanUnionIdx.at(i),
                cleanUnionIdx.at(i + 1) - 1, 1);
        std::vector<int> sIdx = copySubIdx(cleanSegIdx, uIdx);
        //iterate over particular union
        std::vector<int>::iterator it;
        for(it = sIdx.begin(); it != sIdx.end(); ++it) {
            totalCents +=  busybody.getTracker(*it)->getTrLen();
            feat_idx << totalCents << endl;

            std::list<int> iscanList = busybody.getTracker(*it)->getScanList();
            std::list<int> icentList = busybody.getTracker(*it)->getCentroidList();

            std::list<int>::iterator it_s;
            std::list<int>::iterator it_c = icentList.begin();
            //write data to file
            for(it_s = iscanList.begin(); it_s != iscanList.end(); ++it_s) {
                scan_idx << *it_s << endl;
                cent_idx << *it_c  + 1 << endl;
                ++it_c;
            }
        }
    }
    unionfile_idx << cleanUnionIdx.back()  << endl;
}

void SegProc::ttestEq(double xbar1, double xbar2, double n1, double n2, double s12, double s22) {
    v = n1 + n2 - 2;
    double sp2 = ((n1 - 1)*s12 + (n2 - 1)*s22)/v;
    t = (xbar1 - xbar2)/sqrt(sp2*(1/n1 + 1/n2));
}

void SegProc::ttestWelch(double xbar1, double xbar2, double n1, double n2, double s12, double s22) {

    t = (xbar1 - xbar2) / sqrt(s12/n1 + s22/n2);

    double vnumer = (s12 / n1 + s22 / n2 ) * (s12 / n1 + s22 / n2 ); //square the numerator
    double vdenom = (s12 * s12) / (n1*n1 * (n1 - 1)) +
        (s22 * s22) / (n2*n2 * (n2 - 1));

    v = vnumer/vdenom;
}

double SegProc::pt(double x, double n, int lower_tail, int log_p) {
    /* return  P[ T <= x ]  where
     *  * T ~ t_{n}  (t distrib. with n degrees of freedom).
     *
     *   *  --> ./pnt.c for NON-central
     *    */
    double val, nx;
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n))
        return x + n;
#endif
    if (n <= 0.0) ML_ERR_return_NAN;

    if(!R_FINITE(x))
        return (x < 0) ? R_DT_0 : R_DT_1;
    if(!R_FINITE(n))
        return pnorm(x, 0.0, 1.0, lower_tail, log_p);

#ifdef R_version_le_260
    if (n > 4e5) { /*-- Fixme(?): test should depend on `n' AND `x' ! */
        /* Approx. from  Abramowitz & Stegun
         * 26.7.8 (p.949) */
        val = 1./(4.*n);
        return pnorm(x*(1. - val)/sqrt(1. + x*x*2.*val), 0.0, 1.0,
                lower_tail, log_p);
    }
#endif

    nx = 1 + (x/n)*x;
    /* FIXME: This test is probably losing
     * rather than gaining precision,
     *      * now that pbeta(*, log_p = TRUE)
     *      is much better.
     *           * Note however that a version
     *           of this test *is* needed for
     *           x*x > D_MAX */
    if(nx > 1e100) { /* <==>  x*x > 1e100 * n  */
        /* Danger of underflow. So use
         * Abramowitz & Stegun 26.5.4
         *     pbeta(z, a, b) ~ z^a(1-z)^b
         *     / aB(a,b) ~ z^a / aB(a,b),
         *         with z = 1/nx,  a =
         *         n/2,  b= 1/2 :
         *          */
        double lval;
        lval = -0.5*n*(2*log(fabs(x)) - log(n))
            - lbeta(0.5*n, 0.5) - log(0.5*n);
        val = log_p ? lval : exp(lval);
    } else {
        val = (n > x * x)
            ? pbeta (x * x / (n + x * x), 0.5, n / 2., /*lower_tail*/0, log_p)
            : pbeta (1. / nx,             n / 2., 0.5, /*lower_tail*/1, log_p);
    }

    /* Use "1 - v"  if  lower_tail  and  x
     * > 0 (but not both):*/
    if(x <= 0.)
    lower_tail = !lower_tail;

    if(log_p) {
        if(lower_tail) return log1p(-0.5*exp(val));
        else return val - M_LN2; /* = log(.5* pbeta(....)) */
    }
    else {
        val /= 2.;
        return R_D_Cval(val);
    }
}
