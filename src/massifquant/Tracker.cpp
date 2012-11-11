
//Tracker.cpp
#include <stdio.h>
#include <math.h>
#include <algorithm>

#include "OpOverload.h"
#include "Tracker.h"

#include <R.h>
#include <Rdefines.h>

using namespace std;

Tracker::Tracker(const double & init_cent_m,
        const double & init_cent_i,
        const int & scan_num,
        const int & cent_num,
        const double & q_int, const double & q_mz,
        const double & r_int, const double & r_mz,
        const double & ct) {

    mzList.push_back(init_cent_m);
    intensityList.push_back(init_cent_i);

    scanList.push_back(scan_num);
    centroidList.push_back(cent_num);

    q_val_i = q_int;
    q_val_m = q_mz;
    irR = r_int;
    mrR = r_mz;
    criticalT = ct;

    predCounts = 0;
    trLen = 1;
    currMissed = 0;
    mzXbar = 0;
    mzS2 = 0;

    //NO DEPENDENCY INIT
    irXhat = vector<double>(DIMV, 0);
    irF = vector<double>(DIMM, 0);
    irFt = vector<double>(DIMM, 0);
    irH = vector<double>(DIMV, 0);
    irQ = vector<double>(DIMM, 0);
    irP = vector<double>(DIMM, 0);

    mrXhat = vector<double>(DIMV, 0);
    mrF = vector<double>(DIMM, 0);
    mrFt = vector<double>(DIMM, 0);
    mrH = vector<double>(DIMV, 0);
    mrQ = vector<double>(DIMM, 0);
    mrP = vector<double>(DIMM, 0);


    mrXhat[0] = init_cent_m;
    irXhat[0] = init_cent_i;

    //intensity model
    irF[0] = irF[1] = irF[3] = 1;
    irF[2] = 0;

    irFt[0] = irFt[2] = irFt[3] = 1;
    irFt[1] = 0;

    irH[0] = 1; irH[1] = 0;

    irQ[0] = irQ[3] = q_int;
    irQ[1] = irQ[2] = 0;

    irP[0] = q_int;
    irP[3] = P_INIT_I;
    irP[1] = irP[2] = 0;

    //mz model
    mrF[0] = mrF[1] = mrF[3] = 1;
    mrF[2] = 0;

    mrFt[0] = mrFt[2] = mrFt[3] = 1;
    mrFt[1] = 0;

    mrH[0] = 1; mrH[1] = 0;

    mrQ[0] = mrQ[3] =  mrQ[1] = mrQ[2] = 0;

    mrP[0] = q_mz;
    mrP[3] = P_INIT_MZ;
    mrP[1] = mrP[2] = 0;

}

//Destructor
Tracker::~Tracker() {
}

//Debugging
void Tracker::displayContents() {
  //cout << "*******Tracker Contents********" << endl;
  //cout << "predCounts is : " << predCounts << endl;
  //cout << "trLen is : " << trLen << endl;
  //cout << "currMissed is : " << currMissed << endl;
  //cout << "critical is : " << criticalT << endl;
  //cout << "mXhat is : " << mXhat(0) << endl;
  //cout << "mP is : " << mP(0, 0) << endl;
  //cout << "iXhat is : " << iXhat(0) << endl;
  //cout << "iP is : " << iP(0,0) << endl;

  //cout << "Scan List "  << endl;
  //printList(scanList);
  //cout << "Cent List "  << endl;
  //printList(centroidList);
  //cout << "MZ List "  << endl;
  //printList(mzList);
  //cout << "Intensity List "  << endl;
  //printList(intensityList);
}

void Tracker::incrementMiss() {
  currMissed += 1;
}

int Tracker::getCurrMissed() {
  return currMissed;
}

void Tracker::makeZeroCurrMissed() {
   currMissed = 0;
}

void Tracker::setXhat(double m, double i) {
	mrXhat[0] = m;
	irXhat[0] = i;
}

void Tracker::incrementTrLen() {
    trLen += 1;
}

int Tracker::getTrLen() {
 return trLen;
}
int Tracker::getPredCounts() {
    return predCounts;
}

double Tracker::getXbar() {
    return mzXbar;
}

double Tracker::getS2() {
    return mzS2;
}

std::list<int> Tracker::getScanList() {
    return scanList;
}

std::list<int> Tracker::getCentroidList() {
    return centroidList;
}

void  Tracker::appendToTracker(const std::list<int> & sl,
		const std::list<int> & cl,
		const std::list<double> & ml,
		const std::list<double> & il) {
    scanList.insert(scanList.end(), sl.begin(), sl.end());
    centroidList.insert(centroidList.end(), cl.begin(), cl.end());
    mzList.insert(mzList.end(), ml.begin(), ml.end());
    intensityList.insert(intensityList.end(), il.begin(), il.end());
    trLen = scanList.size();
}

int Tracker::getStartScanIdx() {
    return scanList.back();
}

int Tracker::getStopScanIdx() {
    return scanList.front();
}

std::list<double> Tracker::getIntensityList() {
    return intensityList;
}

std::list<double> Tracker::getMzList() {
    return mzList;
}

void Tracker::predictCentroid() {

    /*NO DEPENDENCY VERSION*/
    //mass
    //  mP = mF*mP*mFt + mQ;
    mrP = mrF*mrP*mrFt;
    mrXhat = multiplyMatVec(mrF, mrXhat);
    //intensity
    irP = irF*irP*irFt + irQ;
    irXhat = multiplyMatVec(irF, irXhat);

    predCounts += 1;
}

void Tracker::innovateCentroid(const double & my,
        const double & iy,
        const int scanIdx,
        const int centIdx) {

    /*NO DEPENDENCY VERSION*/
    std::vector<double> mrK(2, 0);
    //mass
    //step1
    mrK[0] = mrP[0] * (1/(mrP[0] + mrR));
   mrK[1] = mrP[2] * (1/(mrP[0] + mrR));
   //step2
    mrXhat[1] =  mrXhat[1] + mrK[1]*(my - mrXhat[0]);
    mrXhat[0] = mrXhat[0] + mrK[0]*(my - mrXhat[0]);
   //step3
    //    /*intermediate vector*/
    std::vector<double> IminusK(4,0);
    IminusK[0] = 1 - mrK[0];
    IminusK[1] = 0;
    IminusK[2] = 0 - mrK[1];
    IminusK[3] = 1;
    mrP = IminusK*mrP;
       //intensity
    std::vector<double> irK(2);
    //step1
    irK[0] = irP[0] * (1/(irP[0] + irR));
    irK[1] = irP[2] * (1/(irP[0] + irR));
    //step2
    //do this first before irXhat[0] changes
    irXhat[1] =  irXhat[1] + irK[1]*(iy - irXhat[0]);
    irXhat[0] = irXhat[0] + irK[0]*(iy - irXhat[0]);
    //step3
    IminusK[0] = 1 - irK[0];
    IminusK[1] = 0;
    IminusK[2] = 0 - irK[1];
    IminusK[3] = 1;
    irP = IminusK*irP;

    //record keeping
    scanList.push_back(scanIdx);
    centroidList.push_back(centIdx);
    mzList.push_back(my);
    intensityList.push_back(iy);

}

int Tracker::claimDataIdx(const std::vector<double> & mData,
			  const std::vector<double> & iData,
			  std::vector<double> & predDist, int minTrLen, int scanBack) {

  //The returned data point index
  int centIdx;

  //Marginal Error
  double mzErrMg = sqrt(mrP[0])*criticalT;
  double left = mrXhat[0] - mzErrMg;
  double right = mrXhat[0] + mzErrMg;

  if ( (trLen >= minTrLen  - 1) && (scanBack == 1) ) {
     lowerList.push_back(left);
     upperList.push_back(right);
  }

  std::vector<double>::const_iterator low,up;
  low = lower_bound(mData.begin(), mData.end(), left);
  up = upper_bound(mData.begin(), mData.end(), right);
  std::vector<int> neighborIdxR;
  int lowint = int(low - mData.begin());
  int upint = int(up - mData.begin());
  if (lowint != upint) {
      neighborIdxR = createSequence(lowint, upint - 1, 1);
  }
  else {
      predDist.push_back(-1);
      centIdx = -1;
      return centIdx;
  }

  //Distance Metric R
  std::vector<double> mSubDataR = copySubIdx(mData, neighborIdxR);
  std::vector<double> iSubDataR = copySubIdx(iData, neighborIdxR);
  std::vector<double> distR = measureDist(mSubDataR, iSubDataR);
  //will be changed inside findMin b/c passed as reference.
  unsigned int centSubIdxR;
  double bestDistR = findMin(distR, centSubIdxR);
  predDist.push_back(bestDistR);
  int centIdxR = neighborIdxR.at(centSubIdxR);

  return centIdxR;
}


std::vector<double> Tracker::getFeatureInfo(double * scanTime) {

    std::vector<double> featInfo(INFOSIZE);
    //1-mz center coordinate is mean
    featInfo[0]=mzXbar;
    //2-mz min
    featInfo[1]= *min_element(mzList.begin(), mzList.end());
    //3-mz max
    featInfo[2]= *max_element(mzList.begin(), mzList.end());
    //4-length
    featInfo[3] = scanList.size();
    //5-scan min
    featInfo[4]  = double(*min_element(scanList.begin(), scanList.end()));
    //featInfo[4]  = scanTime[*min_element(scanList.begin(), scanList.end())];
    //6-scan max
    featInfo[5] = double(*max_element(scanList.begin(), scanList.end()));
    //featInfo[5] = scanTime[*max_element(scanList.begin(), scanList.end())];
    //7-integrated (not normalized intensity)
    std::list<double>::iterator it_i;
    double area = 0;
    double maxInten = 0;
    /*note that sqrt transformation is one to one for values >= 0*/
    //Rprintf("\nIntensity Values\n");
    for(it_i = intensityList.begin(); it_i != intensityList.end(); ++it_i) {
        if (maxInten < *it_i) {
            maxInten = *it_i;
        }
        //Rprintf("%f\t", (*it_i)*(*it_i));
        area += (*it_i)*(*it_i); //resquare it to original form
    }
    featInfo[6]=area;
    featInfo[7]=maxInten * maxInten;

   /*verify some of my findings*/
   /* Rprintf("\nThe Scan Times\n");
    for (int i = scanList.back(); i <= scanList.front(); ++i) {
        Rprintf("%f  ", scanTime[i - 1]);
    }
    Rprintf("\n");

    Rprintf("\nThe Results\n");
    std::vector<double>::iterator it_v;
    for(it_v = featInfo.begin(); it_v != featInfo.end(); ++it_v) {
        Rprintf("%f\t", *it_v);
    }
    */
    return featInfo;
}

std::vector<double> Tracker::measureDist(const vector<double> & mSubData,
					 const vector<double> & iSubData) {

    std::vector<double> d;

    vector<double> mNumerator = mSubData - mrXhat[0];
    vector<double> iNumerator = iSubData - irXhat[0];

    vector<double> mDist = dottimes(mNumerator, mNumerator)/sqrt(mrP[0]);
    vector<double> iDist = dottimes(iNumerator, iNumerator)/sqrt(irP[0]);
    d = dotadd(mDist, iDist);
    return d;
}

double Tracker::findMin(const std::vector<double> & d,
			unsigned int & centIdx) {

  double myMin = d.at(0);
  centIdx = 0;
  unsigned int i;
  for(i = 0; i < d.size(); i++) {
    if (d.at(i) < myMin) {
      myMin = d.at(i);
      centIdx = i;
    }
  }
  return myMin;
}

double Tracker::computeMyXbar() {
    //make it a weighted mean
    std::list<double>::iterator it_m;
    std::list<double>::iterator it_w = intensityList.begin();
    double intenSum = 0;
    for (it_m = mzList.begin(); it_m != mzList.end(); ++it_m) {
        double inten2 = (*it_w) * (*it_w);
        mzXbar += (*it_m) * inten2;
        intenSum += inten2;
        ++it_w;
    }
    mzXbar = mzXbar/intenSum; //mzList.size();
    return mzXbar;
}


double Tracker::getLowerXbar() {
    //defensive measure
    if (lowerList.size() == 0) {
        return mzXbar - 0.1;
    }
    else {
        return  computeAnyXbar(lowerList);
    }
}

double Tracker::getUpperXbar() {
    //defensive measure
    if (upperList.size() == 0) {
        return mzXbar + 0.1;
    }
    else {
        return  computeAnyXbar(upperList);
    }
}

bool Tracker::performScanBack() {

    double lower = getLowerXbar();
    double upper = getUpperXbar();

    std::list<double>::iterator it_m = mzList.begin();
    std::list<double>::iterator it_i = intensityList.begin();
    std::list<int>::iterator it_s = scanList.begin();
    std::list<int>::iterator it_c = centroidList.begin();

    int delCount = 0;
    while (it_m != mzList.end()) {
        //is it outside the converged bounds?
        if (*it_m < lower || upper < *it_m) {
            it_m = mzList.erase(it_m);
            it_i = intensityList.erase(it_i);
            it_s = scanList.erase(it_s);
            it_c = centroidList.erase(it_c);
            delCount++;
        }
        else {
            ++it_m;
            ++it_i;
            ++it_s;
            ++it_c;
        }
    }
    //Rprintf("%d\t", delCount);

    if (delCount > 0) {
        trLen = int(mzList.size());
        return true;
    }
    else {
        return false;
    }
}


double Tracker::computeMyS2() {

    std::list<double>::iterator it;
    for (it = mzList.begin(); it != mzList.end(); ++it) {
        mzS2 += (*it - mzXbar) * (*it - mzXbar);
    }
    mzS2 = mzS2/(mzList.size() - 1);//to make it unbiased
    return mzS2;
}



double Tracker::approxMassAccuracy() {

    std::list<double> diff;

    std::list<double>::iterator it_m;
    for (it_m = mzList.begin(); it_m != mzList.end(); ++it_m) {
        each_ppm.push_back((fabs(*it_m - mzXbar) * MILLION)/mzXbar);

        //cout << "each ppm: " << (fabs(*it_m - mzXbar) * MILLION)/mzXbar << endl;
        //cout << "Mz: " << *it_m << endl;
    }

    //double theMassAcc = computeAnyXbar(each_ppm);
    //cout << "mass accuracy: " << theMassAcc << endl;
    return massAcc;
}
