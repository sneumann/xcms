#ifndef _DYNPROG_H
#define _DYNPROG_H

#include "math.h"

#include "vec.h"
#include "mat.h"

using namespace VEC;


class DynProg {
    private:
        float DEFAULT_GAP_PENALTY_SLOPE;
    public:
        int stupid;
        MatF* _smat;
        MatF _asmat;
        MatI _tb;
        MatI _tbpath;
        MatF _tbscores;  // the tbpath with the score at each position
        MatI _gapmat;
        VecI _mCoords;
        VecI _nCoords;
        VecF _sCoords;
        float _bestScore; // the scores at each m,n coordinate!
        float _prob;

        DynProg() { DEFAULT_GAP_PENALTY_SLOPE = 2.f; }

        // If gap_penalty array len = 0, then a linear gap penalty based on the
        // average matrix score will be used
        // neither diag or gap factor can be 0.0 for minimization
        void find_path(MatF &smat, VecF &gap_penalty, int minimize=0, float diag_factor=2.f,
		       float gap_factor=1.f, int local=0, float init_penalty=0.0f);
        // If gap_penalty array len = 0, then a linear gap penalty based on the
        // average matrix score will be used
        // a gap is introduced without adding in the score of the matrix
        // at that index
        //void find_path_with_gaps(MatF &smat, VecF &gap_penalty, int minimize=0, int local=0, float init_penalty=0.0f);
        void default_gap_penalty(MatF &smat, VecF &out);

        ~DynProg() {}
        // x are the times along the n axis of the tbpath
        // fx are the equivalents along the y axis of the tbpath

        // RETURNS the actual number of internal anchors used
        void warp_map(VecI &mOut, VecI &nOut, float percent_anchors, int minimize=0);
        void best_anchors(VecI &mBijShort, VecI &nBijShort, VecF &sBijShort, VecI &mCoords, VecI &nCoords, VecI &mOut, VecI &nOut, int num_internal_anchors);
        void best_anchors(VecI &mCoordsBijShort, VecI &nCoordsBijShort, VecF &sCoordsBijShort, VecI &mOut, VecI &nOut, int num_internal_anchors);
        void bijective_anchors(VecI &mCoords, VecI &nCoords, VecF &scores, VecI &mBijShort, VecI &nBijShort, VecF &sBijShort);

        // Calculates the sum of the sq of the residuals of the warped nVals
        float sum_sq_res_yeqx(VecF &m_tm, VecF &n_tm, VecI &mWarpMap, VecI &nWarpMap, VecF &mVals, VecF &nVals);
        void path_accuracy_details(VecF &mWarpMapFt, VecF &nWarpMapFt, VecF &mVals, VecF &nVals, VecF &sq_res_yeqx, VecF &abs_diff, int linear_interp=0);
        void path_accuracy(VecF &mWarpMapFt, VecF &nWarpMapFt, VecF &mVals, VecF &nVals, float &sum_sq_res_yeqx, float &avg_sq_res_yeqx, float &sum_abs_diff, float &avg_abs_diff, int linear_interp=0);
        void path_accuracy(VecF &m_tm, VecF &n_tm, VecI &mWarpMap, VecI &nWarpMap, VecF &mVals, VecF &nVals, float &sum_sq_res_yeqx, float &avg_sq_res_yeqx, float &sum_abs_diff, float &avg_abs_diff, int linear_interp=0);
        void _max(float diag, float top, float left, float &val, int &pos);
        void _min(float diag, float top, float left, float &val, int &pos);
        float _global_max(MatF& asmat, int& m_index, int& n_index);
        float _max_right(MatF& asmat, int& m_index);
        float _max_bottom(MatF& asmat, int& n_index);
        float _global_min(MatF& asmat, int& m_index, int& n_index);
        float _min_right(MatF& asmat, int& m_index);
        float _min_bottom(MatF& asmat, int& n_index);
        static int exponential_less_before(float order, int len,  float *expon, float *lessbefore) {
            float val;
            for (int i = 0; i < len; i++) {
                val = pow(i,order);
                expon[i] = val;
                lessbefore[i] = val - expon[i - 1];
            }
			return 1;
        }

        // Score matrices arranged like this
        // nCoords -> scans along the x axis
        // mCoords run |
        //             V
        //             scans along the y axis
        void score_product(MatF &mCoords, MatF &nCoords, MatF &scores);
        void score_covariance(MatF &mCoords, MatF &nCoords, MatF &scores);
        void score_pearsons_r(MatF &mCoords, MatF &nCoords, MatF &scores);

        void score_pearsons_r_opt(MatF &mCoords, MatF &nCoords, MatF &scores);

        void score_mutual_info(MatF &mCoords, MatF &nCoords, MatF &scores, int num_bins=2);
        void score_euclidean(MatF &mCoords, MatF &nCoords, MatF &scores);
        // convenience method for scoring
        void score(MatF &mCoords, MatF &nCoords, MatF &scores, const char *type, int mi_num_bins=2);

//   DynProg::expandFlag(mat1, 2, 1)
//
//       before              after
//    2 0 0 0 0 0         2 2 2 0 0 0
//    0 2 0 0 0 0         2 2 2 2 0 0
//    0 0 2 0 0 0    =>   2 2 2 2 2 0
//    0 0 0 2 0 0         0 2 2 2 2 2
//    0 0 0 0 2 0         0 0 2 2 2 2
//    0 0 0 0 0 2         0 0 0 2 2 2
        static void expandFlag(MatI &flagged, int flag, int numSteps, MatI &expanded);
        // val[i] -= val[i-1]  (inplace)
        void less_before(VecF &arr);
        // linear function mx + b where m is slope and b is y intercept
        // each value is less the value of of the array before it
        void linear_less_before(float m, float b, int len, VecF &lessbefore);
        // linear function mx + b where m is slope and b is y intercept
        void linear(float m, float b, int len,  VecF &arr);
};

#endif
