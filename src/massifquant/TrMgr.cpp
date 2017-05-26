//TrMgr.cpp
#include <algorithm>
#include <fstream>
#include <iostream>
#include "Tracker.h"
#include "TrMgr.h"

using namespace std;

TrMgr::TrMgr(int sidx, const double mi,
                 const int ml, const double cmm,
                 const double mass_acc, const double ct,
                 const int sB)
  : currScanIdx(sidx),
    minIntensity(mi), minTrLen(ml),
    currMissedMax(cmm),
    ppm(mass_acc), criticalT(ct),
    scanBack(sB) {

        initCounts = 0;
}

TrMgr::~TrMgr() {
    int i;
    for(i = 0; i < initCounts; i++) {
        delete trks[i];
    }
}

void TrMgr::setDataScan(const std::vector<double> & mdat,
          const std::vector<double> & idat) {
    //cout << "Setting A New Scan" << endl;
    mData = mdat;
    iData = idat;
}

void TrMgr::setCurrScanIdx(const int sidx) {
    currScanIdx = sidx;
}

void TrMgr::setPredDatIdx(const std::list<int> & pdi) {
    predDatIdx = pdi;
}

void TrMgr::setFoundActIdx(const std::list<int> & fai) {
    foundActIdx = fai;
}

void TrMgr::setMissActIdx(const std::list<int> & mai) {
    missActIdx = mai;
}


void TrMgr::setPredDist(const std::vector<double> & dist) {
    predDist = dist;
}

void TrMgr::setActIdx(const std::vector<int> & ai) {
    actIdx = ai;
}

int TrMgr::getPicCounts() {
    return picIdx.size();
}

int TrMgr::getActiveCounts() {
    return actIdx.size();
}

Tracker* TrMgr::getTracker(int i) {
    return trks[i];
}

std::vector<int> TrMgr::getPicIdx() {
    return picIdx;
}


double TrMgr::getPpm() {
    return ppm;
}

std::vector<double> TrMgr::iterOverFeatures(int i, double * scanTime) {
        return trks[picIdx.at(i)]->getFeatureInfo(scanTime);
}


void TrMgr::predictScan(const vector<double> & mzScan, const vector<double> & intenScan) {

  iData = intenScan;
  mData = mzScan;
  unsigned int i;
  predDatIdx.clear();
  foundActIdx.clear();
  missActIdx.clear();
  predDist.clear();
  int centIdx = -1;
  for (i = 0; i < actIdx.size(); i++) {
    //cout << "ActIdx: " << actIdx.at(i) << endl;
    trks[actIdx.at(i)]->predictCentroid();
    centIdx = trks[actIdx.at(i)]->claimDataIdx(mData,iData,predDist,minTrLen, scanBack);
    //build list of indices corresponding to found or missed
    if (centIdx > -1) {
      foundActIdx.push_back(actIdx.at(i));
      predDatIdx.push_back(centIdx);
    }
    else {
      missActIdx.push_back(actIdx.at(i));
      predDatIdx.push_back(-1);
    }
  }
}

void TrMgr::competeAct() {
  list<int> dupPredDatIdx = predDatIdx; //get copy of original

//  cout << "**********Does New Scan Have Duplicates?***********" << endl;
  //cout << "Act Idx is " << endl;
  /*for (size_t test_act = 0; test_act < actIdx.size(); test_act++) {
      cout << actIdx.at(test_act) << "  ";
  }*/
 //cout << endl;

 //cout << "predDatIdx with misses  is " << endl;
 //printList(dupPredDatIdx);

  predDatIdx = excludeMisses(predDatIdx);
  std::list<int> cpNoMisses = predDatIdx;
 //cout << "Without Misses DatIdx " << endl;
  //printList(predDatIdx);

  int origLen  = predDatIdx.size();
  predDatIdx.sort();
  predDatIdx.unique();
  //cout << "predDatIdx after unique is " << endl;
  //printList(predDatIdx);
  int chgdLen =  predDatIdx.size();
  //check to see if there are duplicates


  if (origLen == chgdLen) {


 //cout << "**********It had no duplicates***********" << endl;
 //cout << "Found List " << endl;
 //printList(foundActIdx);
  //cout << "Buggy Corresponding Found PredDatIdx" << endl;
  //printList(predDatIdx);
  //cout << "Better Corresponding Found PredDatIdx" << endl;
  //printList(cpNoMisses);
  predDatIdx = cpNoMisses;
  //cout << "Missed List " << endl;
  //printList(missActIdx);


      return;
  }
      //cout << "**********It had duplicates***********" << endl;
  std::list<int> restructFoundIdx;
  std::list<int> moreMissIdx;

  std::list<int>::iterator it;
  for (it = predDatIdx.begin(); it != predDatIdx.end(); ++it) {


    std::vector<int> collideIdx = dupPredDatIdx == *it;

    if (collideIdx.size() == 0) {  continue; }
    //cout << "A competing dat pt considered" << endl;
    //cout << "collide idx for datptidx " << *it << endl;
    /*for (size_t test = 0; test < collideIdx.size(); test++) {
        cout << collideIdx.at(test) << "  ";
    }
    cout << endl;     */
    if (collideIdx.size() == 1) {
      restructFoundIdx.push_back(actIdx.at(collideIdx.at(0)));
      continue;
    }
    /*not working*/ // findMinIdx(predDist, collideIdx);
    std::vector<double> subPredDist = copySubIdx(predDist, collideIdx);
    std::vector<double>::iterator bestIt;
    bestIt = min_element(subPredDist.begin(), subPredDist.end());
    int  bestIdx = int(bestIt - subPredDist.begin());
    //cout << "predDist" << endl;
   /* for(size_t test = 0; test < predDist.size(); test++) {
        cout << predDist.at(test) << " ";
    }*/
    //cout << endl;
    //cout << "best subPredIdx : " << bestIdx << endl;
    //cout << "Best Accepted Collide Idx: " << collideIdx.at(bestIdx) << endl;
    restructFoundIdx.push_back(actIdx.at(collideIdx.at(bestIdx)));
    moreMissIdx = collideIdx != collideIdx.at(bestIdx);
    //append moreMissIdx to missActIdx
    std::list<int>::iterator it_miss;
    for (it_miss = moreMissIdx.begin(); it_miss != moreMissIdx.end(); ++it_miss) {
//        //cout << "Missed Collide Idx: " << collideIdx.at(*it_miss) << endl;
        missActIdx.push_back(actIdx.at(collideIdx.at(*it_miss)));
    }
  }
  foundActIdx = restructFoundIdx;
  //cout << "Found List " << endl;
  //printList(foundActIdx);
  //cout << "Corresponding Found PredDatIdx" << endl;
  //printList(predDatIdx);
  //cout << "Missed List " << endl;
  //printList(missActIdx);
}

void TrMgr::manageMissed() {
    std::list<int>::iterator it;
  //  //cout << "Init Counts: " << initCounts << endl;
    for (it = missActIdx.begin(); it != missActIdx.end(); ++it) {
        trks[*it]->incrementMiss();
        //cout << "active tracker id: " << *it << " missed" << endl;
        //consider adding the extra checks like matlab code
        if (trks[*it]->getCurrMissed() > currMissedMax ||
            trks[*it]->getCurrMissed() > trks[*it]->getTrLen() ||
            (trks[*it]->getPredCounts()/2) > trks[*it]->getTrLen()) {
            judgeTracker(*it);
	}
    }
}

void TrMgr::judgeTracker(const int & i) {
    //Perform serial criteria ordered by computational complexity

    //get the index for final deletion
    std::vector<int> subActIdx = actIdx == i;
    //length check
    if (trks[i]->getTrLen() < minTrLen) {
        //cout << "Deleting on account of length: ActIdx is " << i << endl;
        actIdx.erase(actIdx.begin() + subActIdx.at(0));
        delete trks[i];
        trks[i] = NULL;
        return;
    }
    //intensity checks
    std::list<double> iList = trks[i]->getIntensityList();
    double maxHeight = *max_element(iList.begin(), iList.end());
    if (minIntensity > maxHeight) {
        //cout << "Deleting on account of intensity: ActIdx is " << i << endl;
        actIdx.erase(actIdx.begin() + subActIdx.at(0));
        delete trks[i];
        trks[i] = NULL;
        return;
    }
    if (hasMzDeviation(i)) {
        //cout << " Deleted on account of mean deviance" << endl;
        actIdx.erase(actIdx.begin() + subActIdx.at(0));
        delete trks[i];
        trks[i] = NULL;
        return;
    }

    /*
   if (isSeizmo(i)) {
       //cout << " Deleted on account of seizmo" << endl;
       actIdx.erase(actIdx.begin() + subActIdx.at(0));
       delete trks[i];
       trks[i] = NULL;
       return;
   }
    */

    if (scanBack == 1) {
        bool change = trks[i]->performScanBack();
        if (change == true) {
            trks[i]->computeMyXbar();
        }
    }

    //a real pic, have a good retirement.
    picIdx.push_back(i);
    //take off the active list
    //cout << "Its Miss sub index is : " << subActIdx.at(0) << endl;
    actIdx.erase(actIdx.begin() + subActIdx.at(0));
    return;
}

void TrMgr::manageTracked() {



    //assume the two iterators are same size
    std::list<int>::iterator it_f;
    std::list<int>::iterator it_d = predDatIdx.begin();
    for (it_f = foundActIdx.begin(); it_f != foundActIdx.end(); ++it_f) {

        trks[*it_f]->makeZeroCurrMissed();
        trks[*it_f]->incrementTrLen();
        trks[*it_f]->innovateCentroid(mData.at(*it_d),
                                    iData.at(*it_d),
                                    currScanIdx, *it_d);
        //identify for exclusion from new trackers initialized
        mData[*it_d] = CLAIMEDPT;
        iData[*it_d] = CLAIMEDPT;
        ++it_d;
    }
}

void TrMgr::initTrackers(const double & q_int, const double & q_mz,
        const double & r_int, const double & r_mz, const int & sidx) {
//cout << "*********New Scan**************" << endl;
//cout << "scan: " << sidx + 1 << endl;
//cout << " PIC Counts: " << picIdx.size() << endl;
//cout << " Act Counts: " << actIdx.size() << endl;

    currScanIdx = sidx;
    unsigned int i;
    for(i = 0; i < mData.size(); i++) {
        if (mData.at(i) == CLAIMEDPT) { continue; }
//         trks[initCounts] =
       trks.push_back(new Tracker(mData.at(i),iData.at(i),
                currScanIdx, i,
                q_int, q_mz, r_int, r_mz, criticalT));
        actIdx.push_back(initCounts);
        ++initCounts;

        //make sure I allocated enough space
        //assert(initCounts < MAXTRKS);

    }
}

void TrMgr::removeOvertimers() {

    std::vector<int>::iterator it;
    for(it = actIdx.begin(); it != actIdx.end(); ++it) {
        //length check
        if (trks[*it]->getTrLen() < 5) {
            //cout << "Deleting on account of length: ActIdx is " << *it << endl;
            continue;
        }
        //intensity checks
        std::list<double> iList = trks[*it]->getIntensityList();
        double maxHeight = *max_element(iList.begin(), iList.end());

        //used to be default at 60, but that is an arbitrary cut off
        if (minIntensity > maxHeight) {
            //cout << "Deleting on account of intensity: ActIdx is " << *it << endl;
continue; }


        if (hasMzDeviation(*it)) {
//cout << "Deleting on account of mz deviation: ActIdx is " << *it << endl;
            continue; }

/*
        if (isSeizmo(*it)) {
//cout << "Deleting on account of seizmo: ActIdx is " << *it << endl;
            continue; }
*/

        if (scanBack == 1) {
            bool change = trks[*it]->performScanBack();
            if (change == true) {
                trks[*it]->computeMyXbar();
            }
        }

        picIdx.push_back(*it);
    }
    actIdx.clear();
}


void TrMgr::displayTracked() {
    std::vector<int>::iterator it_pic;
    for (it_pic = picIdx.begin(); it_pic != picIdx.end(); ++it_pic) {
        trks[*it_pic]->displayContents();
    }
}

std::list<int> TrMgr::excludeMisses(const std::list<int> & A) {

  std::list<int> clean;
  std::list<int>::const_iterator it;
  for (it = A.begin(); it != A.end(); ++it) {
    if (*it != -1)
      clean.push_back(*it);
  }
  return clean;
}

int TrMgr::findMinIdx(const std::vector<double> & d,
		   const std::vector<int> & idx) {
  std::vector<int>::const_iterator it;
  double tmpmin = d.at(0);
  int min_idx = 1;
  for (it = idx.begin(); it != idx.end(); ++it) {
    //cout << "preddist is" << d.at(*it) << endl;
    if (d.at(*it) < tmpmin) {
      min_idx = *it;
      tmpmin = d.at(*it);
    }
  }
  //cout << "min was " << min_idx << endl;
  return min_idx;
}

void TrMgr::writePICsToFile() {

    ofstream feat_idx("mq_feat_idx.txt");
    ofstream scan_idx("mq_scan_idx.txt");
    ofstream cent_idx("mq_cent_idx.txt");
    ofstream mz_file("mq_mz.txt");
    ofstream intensity_file("mq_intensity.txt");


    /*iterate over all found features listed in picIdx*/
    feat_idx << 0 << endl;;
    int totalCents = 0;
    for (unsigned i = 0; i < picIdx.size(); ++i) {
        totalCents +=  trks[picIdx.at(i)]->getTrLen();
        feat_idx << totalCents << endl;

        std::list<int> iscanList = trks[picIdx.at(i)]->getScanList();
        std::list<int> icentList = trks[picIdx.at(i)]->getCentroidList();
        std::list<double> imzList = trks[picIdx.at(i)]->getMzList();
        std::list<double> iintensityList = trks[picIdx.at(i)]->getIntensityList();

        std::list<int>::iterator it_s;
        std::list<int>::iterator it_c = icentList.begin();
        std::list<double>::iterator it_m = imzList.begin();
        std::list<double>::iterator it_i = iintensityList.begin();

        /*iterate over all scans, indices in tracker (i)*/
        for(it_s = iscanList.begin(); it_s != iscanList.end(); ++it_s) {
            scan_idx << *it_s  << endl;
            cent_idx << *it_c  + 1 << endl;
            mz_file  << *it_m << endl;
            intensity_file << (*it_i)*(*it_i) << endl;
            ++it_c;
            ++it_m;
            ++it_i;
        }
    }
}

void TrMgr::sortPicIdx() {

    //buildl the map of (start values => picIdx);
    //the map will sort itself by keys automatically
    for (size_t i = 0; i < picIdx.size(); ++i) {
        int istart = trks[picIdx.at(i)]->getStartScanIdx();
        startMap[istart] = picIdx.at(i);
    }


    //cout << "picIdx: " << endl;
    for (size_t t = 0; t < picIdx.size(); ++t) {
        //cout << picIdx.at(t) << " ";
    }
    //cout << endl;

    //cout << "startMap: " << endl;
    std::map<int, int>::iterator itm;
    for (itm = startMap.begin(); itm != startMap.end();  ++itm) {
        //cout << itm->first << " " << itm->second << endl;
    }
    //cout << endl;
    //cout << "sorted picIdx: " << endl;

    //restructure the picIdx;
    std::map<int, int>::iterator it;
    int j = 0;
    for (it = startMap.begin(); it != startMap.end(); ++it) {
        picIdx[j] = it->second;
        //cout << picIdx.at(j) << " ";
        j++;
    }
    //cout << endl;





}

void TrMgr::erasePicElements(const std::vector<int> & eIdx) {

    for(size_t i = 0; i < eIdx.size(); ++i) {
        std::vector<int>::iterator it = find(picIdx.begin(), picIdx.end(), eIdx.at(i));
        //if it did find it
        if (it != picIdx.end()) {
            delete trks[eIdx.at(i)];
            trks[eIdx.at(i)] = NULL;
            picIdx.erase(it);
        }
    }
}

std::list<double> TrMgr::diff(const std::list<double> vec) {

    std::list<double> delta;
    std::list<double>::const_iterator it = vec.begin();
    size_t i = 0;
    while (i != vec.size() - 1) {
        delta.push_back(*it - *(++it));
        //cout << "delta.back(): " << delta.back() << endl;
        i++;
    }
    return delta;
}


bool TrMgr::hasMzDeviation(int i) {
    trks[i]->computeMyXbar();
    trks[i]->computeMyS2();
    double mzTol = ppm*trks[i]->getXbar()/10e5;
   //cout << "mzTol is: " << mzTol << endl;
    double peakDeviance = fabs(computeAnyXbar(diff(trks[i]->getMzList())));
    //cout << "Peak Deviance is " << peakDeviance << endl;
    if (peakDeviance > mzTol)
        return true;
    else
        return false;
}

// bool  TrMgr::isSeizmo(int i) {

//     std::list<double> mzList = trks[i]->getMzList();
//     std::vector<double> mzVec(mzList.begin(), mzList.end());
//     std::vector<double> rmz = mzVec; //make copy
//     int midIdx = int(mzList.size()/2);
//     int n = mzList.size() - midIdx;
//     for (int i = 0; i < 3; ++i) {
//         random_shuffle ( rmz.begin(), rmz.end() );
//         std::vector<double> seizmo(n);
//         int k = 0;
//         for (size_t j = midIdx; j < mzVec.size(); ++j) {
//             seizmo[k] = fabs(rmz.at(j) - mzVec.at(j));
//         }
//         for (size_t z = 0; z < seizmo.size(); ++z) {
//             if (seizmo.at(z) > 0.01) { return true; }
//         }
//     }
//     return false;
// }

void TrMgr::shiftUpIndices(const int i) {

    for (size_t j = 0; j < actIdx.size(); ++j) {
        if (i >= actIdx.at(j)) {
            actIdx[j] -= - 1;
        }
    }

    for (size_t j = 0; j < picIdx.size(); ++j) {
        if (i >= picIdx.at(j)) {
            picIdx[j] -= - 1;
        }
    }
}
