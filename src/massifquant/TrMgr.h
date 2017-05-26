#ifndef TRMGR_h
#define TRMGR_h

#include <iostream>
#include <algorithm>
#include <list>
#include <vector>
#include <map>
#include "DataKeeper.h"
#include "OpOverload.h"
#include "Tracker.h"

//const int MAXTRKS = 1e6;
const double CLAIMEDPT = -1;

class TrMgr {

    private:
        int currScanIdx;
        double minIntensity;
        int minTrLen;
        int currMissedMax;
        double ppm;
        double criticalT;
        int scanBack;

        std::vector<double> iData; //intensity for a given scan
        std::vector<double> mData; //mz

        std::vector<Tracker*> trks; //old -> trks[MAXTRKS];
        int initCounts;
        std::vector<int> actIdx;
        std::vector<int> picIdx;
        std::map<int, int> startMap; //first is startScanIdx,
                                    //second is subidx of picIdx.
        int picCounts;
        int actCounts;
        //prediction info
        std::list<int> predDatIdx; //store data points corresponding to claimed tr
        std::list<int> initDatIdx; //mark points for new trackers to be init.
        std::vector<double> predDist; //store distance from claimed tr pred
        std::list<int> foundActIdx; //active index of trs that found
        std::list<int> missActIdx; //active index of trs that missed


        std::list<int> excludeMisses(const std::list<int> & A);

        int findMinIdx(const std::vector<double> & d,
                const std::vector<int> & idx);

        void judgeTracker(const int & i);

        std::list<double> diff(const std::list<double> vec);

        bool hasMzDeviation(int i);

        /* bool isSeizmo(int i); */

    public:


        bool customSort(int i, int j);

        TrMgr(int sidx, const double mi,
                const int ml, const double cmm,
                const double mass_acc, const double ct, const int sB);

        ~TrMgr();

        void setDataScan(const std::vector<double> & mdat,
                const std::vector<double> & idat);
        void setCurrScanIdx(const int sidx);
        void setPredDatIdx(const std::list<int> & pdi);

        void setFoundActIdx(const std::list<int> & fai);
        void setMissActIdx(const std::list<int> & mai);
        void setPredDist(const std::vector<double> & dist);
        void setActIdx(const std::vector<int> & ai);

        int getPicCounts();

        int getActiveCounts();

        Tracker* getTracker(int i);

        std::vector<int> getPicIdx();

        double getPpm();

        std::vector<double> iterOverFeatures(int i, double * scanTime);

        void predictScan(const std::vector<double> & mzScan, const std::vector<double> & intenScan);

        void competeAct();

        void manageMissed();

        void manageTracked();

        void initTrackers(const double & q_int, const double & q_mz,
                const double & r_int, const double & r_mz,
                const int & sidx);

        void removeOvertimers();

        void displayTracked();

        void writePICsToFile();

        void sortPicIdx();

        void erasePicElements(const std::vector<int> & eIdx);

        void shiftUpIndices(const int i);

};


#endif
