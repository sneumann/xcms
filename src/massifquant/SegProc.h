#ifndef SPROC_h
#define SPROC_h

#include "nmath.h"
#include "dpq.h"

#include <list>
#include <vector>

const int MAXGAP = 5;
const double TROBUST1 = 0.5;
const double TROBUST2 = 2;
const double ALPHA = 0.001;

class SegProc {

    private:

        double origTrNum; //num of trackers before soldering of trackers
        double t; //the test statistic from a student t distribution
        double v; //degrees of freedom for student t distribution
        double p; //the probability of obtaining a test statistic at least
                  //as extreme as the one that was actually observed, assuming
                  //that the null hypothesis is true (wiki)

        std::vector<int> segClusters;
        std::vector<int> segIdx;//the list of tracker indices to be soldered together
        std::vector<int> unionIdx; //the delimiation of unions of soldered trks

        //after cleaned subsets
        std::vector<int> cleanSegIdx;//the list of tracker indices to be soldered together
        std::vector<int> cleanUnionIdx; //the delimiation of unions of soldered trks



        void compareMeans(TrMgr & busybody, const int seed, const std::list<int> edges, const int & segCounts);

        //equal variances
        void ttestEq(double xbar1, double xbar2, double n1, double n2, double s12, double s22);

        //unequal variance - Welch-Satterwaite Approx.
        void ttestWelch(double xbar1, double xbar2, double n1, double n2, double s12, double s22);

        double pt(double x, double n, int lower_tail, int log_p);

    public:

        SegProc(int otn);

        ~SegProc();

        void groupSegments(TrMgr & busybody);

        void collapseSubsets();

        void solderSegs(TrMgr & busybody);

        void segsToFile(TrMgr & busybody);
};

#endif
