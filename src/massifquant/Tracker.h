#ifndef TR_h
#define TR_h

#include<iostream>
#include <vector>
#include <list>
//Tracker.h

//Global Constants
const int DIMV = 2;
const int DIMM = 4;
const int INFOSIZE = 8;
const double P_INIT_I = 10000;
const double P_INIT_MZ = 0.000001; //10^(-6);
const double MILLION = 1000000;
class Tracker {

    private:

        std::list<int> centroidList;
        std::list<int> scanList;

        std::list<double> intensityList;
        std::list<double> mzList;
        std::list<double> lowerList;
        std::list<double> upperList;
        std::list<double> each_ppm;

        int predCounts;
        int trLen;
        int currMissed;
        double criticalT;

        double mzXbar; //mean
        double mzS2;   //variance
        double massAcc;

        //Model Specs
        /*intensity*/
        std::vector<double> irXhat;
        std::vector<double> irF;
        std::vector<double> irFt;
        std::vector<double> irH;
        std::vector<double> irQ;
        double irR;
        std::vector<double> irP;

        double q_val_i; //process uncertainty scalar
        double r_val_i;
        double p_val_i;

        double q_val_m; //process uncertainty scalar
        double r_val_m;
        double p_val_m;

        /*mass*/
        std::vector<double> mrXhat;
        std::vector<double> mrF;
        std::vector<double> mrFt;
        std::vector<double> mrH;
        std::vector<double> mrQ;
        double mrR;
        std::vector<double> mrP;

        //methods
        double getLowerXbar();

        double getUpperXbar();

        std::vector<double> measureDist(const std::vector<double> & mSubData,
                const std::vector<double> & iSubData);

        double findMin(const std::vector<double> & d,
                unsigned int & centIdx);
    public:

        //Constructor
        Tracker(const double & init_cent_m,
                const double & init_cent_i,
                const int & scan_num,
                const int & cent_num,
                const double & q_int, const double & q_mz,
                const double & r_int, const double & r_mz,
                const double & ct);

        //Destructor
        ~Tracker();

        void incrementMiss();

        int getCurrMissed();

        void makeZeroCurrMissed();

        void incrementTrLen();

        int getTrLen();

        int getPredCounts();

        double getXbar();

        double getS2();

        std::list<int> getScanList();

        std::list<int> getCentroidList();

        void appendToTracker(const std::list<int> & sl,
                const std::list<int> & cl,
		const std::list<double> & ml,
		const std::list<double> & il);

        int getStartScanIdx();

        int getStopScanIdx();

        void setXhat(double m, double i);

        double computeMyXbar();

        double computeMyS2();

        double approxMassAccuracy();

        std::list<double> getIntensityList();

        std::list<double> getMzList();

        void displayContents();

        void predictCentroid();

        void innovateCentroid(const double & my,
                const double & iy,
                const int scanIdx,
                const int centIdx);

        int claimDataIdx(const std::vector<double> & mData,
                const std::vector<double>  & iData,
                std::vector<double> & predDist,
                int minTrLen, int scanBack);

        bool performScanBack();

        std::vector<double> getFeatureInfo(double * scanTime);
};


void printList(const std::list<int> & mylist);

void printList(const std::list<double> & mylist);

#endif
