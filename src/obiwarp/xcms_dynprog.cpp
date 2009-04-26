#include <cstdlib>
#include <cstring>
#include <iostream>
#include<math.h>

#include "xcms_dynprog.h"
#include "vec.h"
#include "mat.h"
#include "assert.h"

#ifndef min
	#define min(a,b) ( ( (a) < (b) ) ? (a) : (b) )
#endif
#ifndef max
	#define max(a,b) ( ( (a) > (b) ) ? (a) : (b) )
#endif


#include <fstream>
using namespace std;

using namespace VEC;

int pngCnt = 0;
float _LOG2 = logf(2);

/************************************************************
 * DECLARE HELPER FUNCTIONS
 ************************************************************/
float sum(MatF &mat, int rowNum);
float sumXSquared(MatF &mat, int rowNum);
float sumOfProducts(MatF &mat1, int rowNum1, MatF &mat2, int rowNum2);
void _subtract(MatF &mat, int rowNum, float val, MatF &minused);
float entropy(MatF &mat, int rowNum, int numBins, float minVal, float scaleFactor, MatI &indArray);
void entropyXY(MatI &binIndX, MatI &binIndY, VecF &entropyX, VecF &entropyY, MatF &scores, int numBins);

void _traceback(MatI &tb, MatF &smat, int m, int n, MatI &tbpath, VecI &equiv1, VecI &equiv2, VecF &scores);

// NEED to redo these guys and affirm correctness
/*
//void DynProg::warp_diag(VecI &mCoords, VecI &nCoords, MatF &mMat, MatF &warpedOut) {
void DynProg::warp(VecI &mCoords, VecI &nCoords, VecF &mVec, VecF &warpedOut, bool mCoord_row_nums) {
    int i;
    VecF mCoordsF(mCoords.length());
    VecF nCoordsF(nCoords.length());
    for (i=0;i<mCoords.length();++i) { mCoordsF[i] = (float)(mCoords[i]); }
    for (i=0;i<nCoords.length();++i) { nCoordsF[i] = (float)(nCoords[i]); }
    warp(mCoordsF, nCoordsF, mVec, warpedOut, mCoord_row_nums);
}

void DynProg::warp(VecF &mCoords, VecF &nCoords, VecF &mVec, VecF &warpedOut, bool mCoord_row_nums) {
    int i;

    VecF new_y;
    // get the indices array:
    VecF new_x(mVec.length());
    for (i = 0; i < mVec.length(); ++i) { new_x[i] = (float)i; }
    VecF::chfe(mCoords, nCoords, new_x, new_y);

    int nCoords_last_val = (int)(nCoords[nCoords.length()-1]);  
    //printf("ncoordslastval: %d\n", nCoords_last_val);
    VecF warp_x(nCoords_last_val+1); // nCoords_last_val is the last INDEX!
    
    if (mCoord_row_nums) {
        int mCoords_last_val = (int)mCoords[mCoords.length()-1];
        VecF tmp_warp_x(mCoords_last_val+1);
        double inc = (double)nCoords_last_val/(double)mCoords_last_val;
        for (i = 0; i < mCoords_last_val+1; ++i) {
            tmp_warp_x[i] = ((double)i)*inc;
        }
        warp_x.take(tmp_warp_x);
    }
    else {
        for (i = 0; i < nCoords_last_val+1; ++i) {
            warp_x[i] = (float)i;
        }
    }
     
    VecF::chfe(new_y, mVec, warp_x, warpedOut);
}



void DynProg::warp(VecI &mCoords, VecI &nCoords, MatF &mMat, MatF &warpedOut, bool mCoord_row_nums) {
    int i;
    VecF mCoordsF(mCoords.length());
    VecF nCoordsF(nCoords.length());
    for (i=0;i<mCoords.length();++i) { mCoordsF[i] = (float)(mCoords[i]); }
    for (i=0;i<nCoords.length();++i) { nCoordsF[i] = (float)(nCoords[i]); }

    VecF new_y;
    // get the indices array:
    VecF new_x(mMat.rows());
    for (i = 0; i < mMat.rows(); ++i) { new_x[i] = (float)i; }
    VecF::chfe(mCoordsF, nCoordsF, new_x, new_y);

    int nCoords_last_val = (int)(nCoords[nCoords.length()-1]);  
    //printf("ncoordslastval: %d\n", nCoords_last_val);
    VecF warp_x(nCoords_last_val+1); // nCoords_last_val is the last INDEX!
    
    if (mCoord_row_nums) {
        int mCoords_last_val = (int)mCoords[mCoords.length()-1];
        VecF tmp_warp_x(mCoords_last_val+1);
        double inc = (double)nCoords_last_val/(double)mCoords_last_val;
        for (i = 0; i < mCoords_last_val+1; ++i) {
            tmp_warp_x[i] = ((double)i)*inc;
        }
        warp_x.take(tmp_warp_x);
    }
    else {
        for (i = 0; i < nCoords_last_val+1; ++i) {
            warp_x[i] = (float)i;
        }
    }
     
    MatF trans;
    mMat.transpose(trans);
    
    int transcols = trans.cols();
    int transrows = trans.rows();

    int warp_x_length = warp_x.length(); 
    MatF newValsTrans(transrows, warp_x_length); 

    for (i = 0; i < transrows; ++i) {
        VecF warp_y(transcols, newValsTrans.pointer(i), 1);
        VecF chrm_y(transcols, trans.pointer(i), 1);
        VecF::chfe(new_y, chrm_y, warp_x, warp_y);
        //warp_y.print(); 
    } 
    newValsTrans.transpose(warpedOut);
}
*/

void DynProg::warp_map(VecI &mOut, VecI &nOut, float percent_anchors, int minimize) {

    // If minimize, then we maximize the negative!
    if (minimize) {
        _sCoords*= -1.f;
    }

    VecI mBijS, nBijS;
    VecF sBijS;
    bijective_anchors(_mCoords, _nCoords, _sCoords, mBijS, nBijS, sBijS);
    // percent_anchors:
    // 0 = 2 anchors on the ends (0 internal anchors)
    // 100 = all points (total_size - 2)
    
    // num_internal spots
    double D_num_internal_anchs = (percent_anchors/100)*(double)mBijS.size();  

    // Round to nearest int
    int num_internal_anchs = (int)D_num_internal_anchs;
    if ((int)D_num_internal_anchs != (int)(D_num_internal_anchs + 0.5)) {
        num_internal_anchs++; 
    }

    if (minimize) {
        _sCoords*= -1.f;
    }
    best_anchors(mBijS, nBijS, sBijS, _mCoords, _nCoords, mOut, nOut, num_internal_anchs);
}


void DynProg::bijective_anchors(VecI &mCoords, VecI &nCoords, VecF &scores, VecI &mBijShort, VecI &nBijShort, VecF &sBijShort) {
    
    int i;
    // number of true peak equivalencies (not counting gaps) 
    // guarantees that the first and last equivalencies will be the same
    //puts("equivs:"); mCoords.print(); nCoords.print(); scores.print();
    int lengthEquiv = mCoords.dim();

    int *nC_arr = new int[lengthEquiv];  // Max length...
    int *mC_arr = new int[lengthEquiv];  // Max length...
    float *sC_arr = new float[lengthEquiv];  // Max length...

    // monoMaps do NOT include the first and last equivalents!

    // find the length of the equivalents for MONOTONIC pairs (gaps collapse)

    // Eliminate the equivalents that have the same m or n as the two corners
    int *mCoordsS_arr = new int[lengthEquiv-2];
    int *nCoordsS_arr = new int[lengthEquiv-2];
    float *scoresS_arr = new float[lengthEquiv-2];
    int last_index = lengthEquiv - 1;
    int startm = mCoords[0];
    int startn = nCoords[0];
    int endm = mCoords[last_index];
    int endn = nCoords[last_index];
    int newLengthS = 0;
    for (i = 1; i < lengthEquiv; ++i) {
        if (mCoords[i] == startm || mCoords[i] == endm || nCoords[i] == startn || nCoords[i] == endn) {
            // Toss out those guys
        }
        else {
            mCoordsS_arr[newLengthS] = mCoords[i];
            nCoordsS_arr[newLengthS] = nCoords[i];
            scoresS_arr[newLengthS] = scores[i];
            ++newLengthS;
        }
    }

    VecI mCoordsS(newLengthS, mCoordsS_arr);
    VecI nCoordsS(newLengthS, nCoordsS_arr);
    VecF scoresS(newLengthS, scoresS_arr);

    //puts("WITHOUT COMPETTION"); mCoordsS.print(); nCoordsS.print(); scoresS.print();

    //puts("MARCHING");
    int I = -1;
    int prevN = -1;
    int prevM = -1;
    int prevN2 = -2;
    int prevM2 = -2;
    for (i = 0; i < newLengthS; ++i) {
        if (prevN != nCoordsS[i] && prevM != mCoordsS[i]) {  //diagonal
            //printf("%d %d DIAGONAL NEW COMPETITION\n", i, I);
            ++I;
            mC_arr[I] = mCoordsS[i];
            nC_arr[I] = nCoordsS[i];
            sC_arr[I] = scoresS[i];
        }
        else if (prevN != nCoordsS[i]) {  // going right
            if (prevN2 == prevN) {  // This is a CHANGE in direction to down
                //printf("mC_arr[I] = %d mCoordsS[i] = %d\n", mC_arr[I],mCoordsS[i]);
                if (mC_arr[I] == mCoordsS[i]) {  // is this m coord spoken for
                    //printf("%d %d M COORD SPOKEN FOR \n", i, I);
                }
                else {  // start a new competition
                    //printf("%d %d going right NEW COMPETITION \n", i, I);
                    ++I;
                    mC_arr[I] = mCoordsS[i];
                    nC_arr[I] = nCoordsS[i];
                    sC_arr[I] = scoresS[i];
                }
            }
            else {  // continuing right
                if (scoresS[i] >= sC_arr[I]) {  // is this a better score?
                    //printf("%d %d A BETTER SCORE! (CONTINUING RIGHT) \n", i, I);
                    mC_arr[I] = mCoordsS[i];
                    nC_arr[I] = nCoordsS[i];
                    sC_arr[I] = scoresS[i];
                }
            }
        }
        else {  // if (prevM != mCoordsS[I])  // going down
            if (prevM2 == prevM) {  // This is a CHANGE in direction to down
                if (nC_arr[I] == nCoordsS[i]) {  // is this n coord spoken for
                    //printf("%d %d N COORD SPOKEN FOR \n", i, I);
                }
                else {  // start a new competition
                    //printf("%d %d going DOWN NEW COMPETITION \n", i, I);
                    ++I;
                    mC_arr[I] = mCoordsS[i];
                    nC_arr[I] = nCoordsS[i];
                    sC_arr[I] = scoresS[i];
                }
            }
            else {  // continuing down
                if (scoresS[i] >= sC_arr[I]) {  // is this a better score?
                    //printf("%d %d A BETTER SCORE! (CONTINUING DOWN) \n", i, I);
                    mC_arr[I] = mCoordsS[i];
                    nC_arr[I] = nCoordsS[i];
                    sC_arr[I] = scoresS[i];
                }
            }
        }
        prevN2 = prevN;
        prevM2 = prevM;
        prevN = nCoordsS[i];
        prevM = mCoordsS[i];
    }

    //puts("FIRST ROUND"); VecI mC(I, mC_arr, 1); VecI nC(I, nC_arr, 1); VecF sC(I, sC_arr, 1);
    //mC.print(); nC.print(); sC.print();

    //printf("FINAL I: %d\n", I);

    // mC_arr = short!
    mBijShort.take(I, mC_arr);
    nBijShort.take(I, nC_arr);
    sBijShort.take(I, sC_arr);
    //puts("mBijshort: ");mBijShort.print();
    //puts("nBijshort: ");nBijShort.print();
    //puts("sBijshort: ");sBijShort.print();
}


void DynProg::best_anchors(VecI &mBijShort, VecI &nBijShort, VecF &sBijShort, VecI &mCoords, VecI &nCoords, VecI &mOut, VecI &nOut, int num_internal_anchors) {
    // Need mCoords and nCoords simply to set the ends of the anchors!

    int i;
    int shortLength = mBijShort.size();
    int trueLength = shortLength + 2;
    // Check that the num_internal_anchors+2
    if (num_internal_anchors+2 > trueLength) {
        std::cerr << "changing " << num_internal_anchors << " num_internal_anchors to " << trueLength -2 << " to be inbounds";
        num_internal_anchors = trueLength - 2;
    }

    // Setup the new equivalency arrays
    int newLen = num_internal_anchors+2;
    //printf("NEWLEN: %d\n", newLen);

    VecI newMcoords(newLen);
    VecI newNcoords(newLen);

    // Set the ends:
    newMcoords[0] = mCoords[0];
    newNcoords[0] = nCoords[0];
    newMcoords[newLen-1] = mCoords[mCoords.size()-1];
    newNcoords[newLen-1] = nCoords[nCoords.size()-1];

    int numSegments = num_internal_anchors;
    float segLength = ((float)(shortLength)/numSegments);  

    int fmark;
    int smark;
    float maxScore;

    //printf("numSegments: %d\n", numSegments);
    for (int outer = 0; outer < numSegments; outer++) {
        //printf("##############################################\n");
        //printf("OUT LOOP: %d\n", outer);

        fmark = (int)(((double)(outer+1))*(double)segLength); 
        smark = (int)(((double)(outer))*(double)segLength); 
        maxScore = sBijShort[smark];
        for (i = smark; i < fmark; ++i) {
            //printf("************************\n");
            //printf("smark: %d\n", smark);
            //printf("fmark: %d\n", fmark);
            //printf("I: %d\n", i);
            //printf("************************\n");
            if (sBijShort[i] >= maxScore) {
                maxScore = sBijShort[i];
                //printf("maxscore: %f\n", maxScore);
                newMcoords[outer+1] = mBijShort[i];
                newNcoords[outer+1] = nBijShort[i];
            }
        }
    }

    //warpMap = tmpWarpMap;
    mOut.take(newMcoords);
    nOut.take(newNcoords);
}

// n and m are the starting index for traceback
// Sets tbpath to be a matrix of 1's and 0's showing the path
// scores are the scores at each equivalent point

void _traceback(MatI &tb, MatF &smat, int m, int n, MatI &tbpath, VecI &mCoord, VecI &nCoord, VecF &scores) {
    int i;
    // Vec2I tbvec   to add in future for the equivalencies
    int rows = tb.rows();
    int cols = tb.cols();
    int *n_eqr = new int[rows+cols];
    int *m_eqr = new int[rows+cols];
    float *score_pathr = new float[rows+cols];

    int cnt = 0;
    // Go until we jump off the matrix:
    while (m != -1 && n != -1) {
        //if (cnt > rows+cols) { printf("rows: %d cols: %d\n", rows, cols); printf("cnt too big!"); exit(1); }
        n_eqr[cnt] = n;
        m_eqr[cnt] = m;
        tbpath(m,n) = 1;
        score_pathr[cnt] = smat(m,n);

        int val = tb(m,n);  
        if (val == 0) { // Diag
            m -= 1;
            n -= 1;
        }
        else if (val == 1) { // UP
            m -= 1;
        }
        else {  // val == 2  // Left
            n -= 1;
        }
        cnt++;
    }

    int *tmpEquiv_m = new int[cnt];
    int *tmpEquiv_n = new int[cnt];
    float *tmpScores = new float[cnt];

    // Reverse the arrays
    int rev = cnt - 1;
    // This could be made a little faster...
    for (i = 0; i < cnt; ++i) {
        tmpEquiv_m[i] = m_eqr[rev];
        tmpEquiv_n[i] = n_eqr[rev];
        tmpScores[i] = score_pathr[rev];
        rev--;
    }
    delete[] n_eqr;
    delete[] m_eqr;
    delete[] score_pathr;
    mCoord.take(cnt, tmpEquiv_m);
    nCoord.take(cnt, tmpEquiv_n);
    scores.take(cnt, tmpScores);
}

/*
// Replaces all locations marked by 1 in toReplace with values from a different
// location 
void DynProg::replaceAlignmentPathRandom(MatF& mat, MatI& toReplace) {
float *unrolledSqeezed = new float[mat.rows()*mat.cols()];
// Create a single vec of the OK vals:
unsigned int *_m = new unsigned int[mat.rows()*mat.cols()];
unsigned int *_n = new unsigned int[mat.rows()*mat.cols()];
unsigned int cnt = 0;
unsigned int otherCnt = 0;

    for (int m = 0; m < mat.rows(); ++m) {
        for (int n = 0; n < mat.cols(); ++n) {
            if (!(toReplace(m,n))) {
                unrolledSqeezed[cnt] = mat(m,n);
                ++cnt;
            }
            else {
                _m[otherCnt] = m;
                _n[otherCnt] = n;
                otherCnt++;
            }
        }
    }
    // DOING IT WITH GOOD RAND:
    // Turn this stuff off if we want scores to always be exactly the same...
    //time_t _seconds;
    //time(&_seconds);
    //goodRandSeed((unsigned int)_seconds);

    // DOING IT WITH rand:
    srand( time(NULL) );

    for (int i = 0; i < otherCnt; ++i) {
        //mat(_m[i],_n[i]) = unrolledSqeezed[goodRandI(0, cnt-1)];
		int _mind = _m[i];
		int _nind = _n[i];
        mat(_mind,_nind) = unrolledSqeezed[rand()%cnt];
    }
}
*/

// Fetch the highest value on the right and bottom sides and sets 
// m_index and n_index at that index
float DynProg::_global_max(MatF& asmat, int& m_index, int& n_index) {
    float max_score;
    float max_right = DynProg::_max_right(asmat, m_index);
    float max_bottom = DynProg::_max_bottom(asmat, n_index);
    if (max_right > max_bottom) {
        max_score = max_right;
        n_index = asmat.cols() - 1;
    }
    else {
        max_score = max_bottom;
        m_index = asmat.rows() - 1;
    }
    return max_score; 
}

// Fetch the highest value on the right and bottom sides and sets 
// m_index and n_index at that index
float DynProg::_global_min(MatF& asmat, int& m_index, int& n_index) {
    float min_score;
    float min_right = DynProg::_min_right(asmat, m_index);
    float min_bottom = DynProg::_min_bottom(asmat, n_index);
    if (min_right < min_bottom) {
        min_score = min_right;
        n_index = asmat.cols() - 1;
    }
    else {
        min_score = min_bottom;
        m_index = asmat.rows() - 1;
    }
    return min_score; 
}

// Returns the max value and also sets m_index at the max index (0 indexed)
float DynProg::_max_right(MatF& asmat, int& m_index) {
    int cols = asmat.cols();
    int rows = asmat.rows();
    int n = cols - 1;
    float max = asmat(0,n);
    for (int m = 0; m < rows; ++m) {
        if (asmat(m,n) >= max) {
            max = asmat(m,n);
            m_index = m;
        }
    }
    return max;
}

float DynProg::_max_bottom(MatF& asmat, int& n_index) {
    int cols = asmat.cols();
    int rows = asmat.rows();
    int m = rows - 1;
    float max = asmat(m,0);
    for (int n = 0; n < cols; ++n) {
        if (asmat(m,n) >= max) {
            max = asmat(m,n);
            n_index = n;
        }
    }
    return max;
}

// Returns the min value and also sets m_index at the min index (0 indexed)
float DynProg::_min_right(MatF& asmat, int& m_index) {
    int cols = asmat.cols();
    int rows = asmat.rows();
    int n = cols - 1;
    float min = asmat(0,n);
    for (int m = 0; m < rows; ++m) {
        if (asmat(m,n) <= min) {
            min = asmat(m,n);
            m_index = m;
        }
    }
    return min;
}

float DynProg::_min_bottom(MatF& asmat, int& n_index) {
    int cols = asmat.cols();
    int rows = asmat.rows();
    int m = rows - 1;
    float min = asmat(m,0);
    for (int n = 0; n < cols; ++n) {
        if (asmat(m,n) <= min) {
            min = asmat(m,n);
            n_index = n;
        }
    }
    return min;
}



// Precedence on tie goes to leftmost argument
// Finds the max value and positions (0,1,2)
// Will accept at least two NULL arguments in diag,top,left
void DynProg::_max(float diag, float top, float left, float &val, int &pos) {
    val = diag;
    if (diag >= top) {
        if (diag >= left) {
            // Already set val to diag
            pos = 0;
            return;
        }
        else {
            val = left;
            pos = 2;
        }
    }
    else if (top >= left) {
        val = top;
        pos = 1;
    }
    else {
        val = left;
        pos = 2;
    }
    //printf("max = %f pos = %d of %f,%f,%f\n", val, pos, diag, top, left);
}

// Precedence on tie goes to leftmost argument
// Finds the min value and positions (0,1,2)
// Will accept at least two NULL arguments in diag,top,left
void DynProg::_min(float diag, float top, float left, float &val, int &pos) {
    val = diag;
    if (diag <= top) {
        if (diag <= left) {
            // Already set val to diag
            pos = 0;
            return;
        }
        else {
            val = left;
            pos = 2;
        }
    }
    else if (top <= left) {
        val = top;
        pos = 1;
    }
    else {
        val = left;
        pos = 2;
    }
    //printf("min = %f pos = %d of %f,%f,%f\n", val, pos, diag, top, left);
}




/*
// Requires that the initial dynamic programming was previously done
// i.e. must have:
//      _smat (pointer)
//      _tbpath
// the gap_penalty, type, init_penalty, minimum should all be the same
// as was used to initialize the DynProg object
float DynProg::toProb(int halfWindow, short int numShuffles, char *type, float init_penalty, int minimum) {
    int _PATHFLAG = 1;
    MatI _expanded(_tbpath.rows(), _tbpath.cols());
    expandFlag(_tbpath, _PATHFLAG, halfWindow, _expanded);
    MatF _scorecpy = _smat->copy();
    replaceAlignmentPathRandom(_scorecpy, _expanded);
    VecF _unrolled(_scorecpy);
    VecF randomScores(numShuffles);
    for (int i = 0; i < numShuffles; ++i) {
        _unrolled.shuffle();
        //_unrolled.print();
        //_scorecpy.print();
        DynProg tmp;
        tmp.findPath(_scorecpy, type, init_penalty, minimum);
        randomScores[i] = tmp._bestScore;
    }
    //float mean;
    //float stddev;
    return randomScores.probOneSideRight(_bestScore);
    // Vec2 ---> 0 = n; 1 = m
}
*/

void DynProg::path_accuracy_details(VecF &mWarpMapFt, VecF &nWarpMapFt, VecF &mVals, VecF &nVals, VecF &sq_res_yeqx, VecF &abs_diff, int linear_interp) {
    VecF nValsNew;
    if (linear_interp) {
        VecF::linear_interp(mWarpMapFt, nWarpMapFt, mVals, nValsNew);    
    }
    else {
        VecF::chfe(mWarpMapFt, nWarpMapFt, mVals, nValsNew);    
    }

    // sq res
    VecF nVals_tmp1;
    nVals.copy(nVals_tmp1);
    nVals_tmp1 -= nValsNew;
    VecF nVals_tmp1_sq;
    nVals_tmp1.mul(nVals_tmp1, nVals_tmp1_sq);
    nVals_tmp1_sq /= 2;
    sq_res_yeqx.take(nVals_tmp1_sq);

    // abs diffs:
    VecF nVals_tmp2;
    nVals.copy(nVals_tmp2);
    nVals_tmp2 -= nValsNew;
    nVals_tmp2.abs_val();
    abs_diff.take(nVals_tmp2);
}

void DynProg::path_accuracy(VecF &mWarpMapFt, VecF &nWarpMapFt, VecF &mVals, VecF &nVals, float &sum_sq_res_yeqx, float &avg_sq_res_yeqx, float &sum_abs_diff, float &avg_abs_diff, int linear_interp) {

    VecF nValsNew;
    if (linear_interp) {
        VecF::linear_interp(mWarpMapFt, nWarpMapFt, mVals, nValsNew);    
    }
    else {
        VecF::chfe(mWarpMapFt, nWarpMapFt, mVals, nValsNew);    
    }

    //calculate the sum of the sq of the residuals:
    sum_sq_res_yeqx = VecF::sum_sq_res_yeqx(nVals, nValsNew);
    avg_sq_res_yeqx = sum_sq_res_yeqx/nVals.length(); 
    VecF diff;
    nVals.sub(nValsNew, diff);
    diff.abs_val();  
    sum_abs_diff = diff.sum();
    avg_abs_diff = sum_abs_diff/nVals.length();
}


void DynProg::path_accuracy(VecF &m_tm, VecF &n_tm, VecI &mWarpMap, VecI &nWarpMap, VecF &mVals, VecF &nVals, float &sum_sq_res_yeqx, float &avg_sq_res_yeqx, float &sum_abs_diff, float &avg_abs_diff, int linear_interp) {
    // mVals are the m time standard points
    // nVals are the n time standard points

    // warp the nVals:
    VecF _mWarpMapFt(mWarpMap.length());
    VecF _nWarpMapFt(nWarpMap.length());
    for (int i = 0; i < mWarpMap.length(); ++i) {
        if (mWarpMap[i] < 0 || mWarpMap[i] >= m_tm.length()) {
            std::cerr << "ASKING FOR VAL OUTSIDE RANGE, length: " << m_tm.length() << " requested: " << mWarpMap[i] << "\n";
        }
        _mWarpMapFt[i] = m_tm[mWarpMap[i]];
        _nWarpMapFt[i] = n_tm[nWarpMap[i]];
    }
    path_accuracy(_mWarpMapFt, _nWarpMapFt, mVals, nVals, sum_sq_res_yeqx, avg_sq_res_yeqx, sum_abs_diff, avg_abs_diff, linear_interp);
}


float DynProg::sum_sq_res_yeqx(VecF &m_tm, VecF &n_tm, VecI &mWarpMap, VecI &nWarpMap, VecF &mVals, VecF &nVals) {
    // mVals are the m time standard points
    // nVals are the n time standard points

    // warp the nVals:
    VecF mWarpMapFt(mWarpMap.length());
    VecF nWarpMapFt(nWarpMap.length());
    for (int i = 0; i < mWarpMap.length(); ++i) {
        if (mWarpMap[i] < 0 || mWarpMap[i] >= m_tm.length()) {
            std::cerr << "ASKING FOR VAL OUTSIDE RANGE, length: " << m_tm.length() << " requested: " << mWarpMap[i] << "\n";
        }
        mWarpMapFt[i] = m_tm[mWarpMap[i]];
        nWarpMapFt[i] = n_tm[nWarpMap[i]];
    }
    VecF nValsNew;
    VecF::chfe(mWarpMapFt, nWarpMapFt, mVals, nValsNew);    

    //calculate the sum of the sq of the residuals:
    return VecF::sum_sq_res_yeqx(nVals, nValsNew);
}

void DynProg::score_product(MatF &mCoords, MatF &nCoords, MatF &scores) {
    int i;
    int s_mlen = mCoords.rows();// s_cols = length_n
    int s_nlen = nCoords.rows();// s_rows = length_m  // Both rows and cols derived from # rows
    int cols = mCoords.cols();
    assert(cols == nCoords.cols());
    MatF tmp(s_mlen, s_nlen);
    for (int m = 0; m < s_mlen; ++m) {
        for (int n = 0; n < s_nlen; ++n) {
            float sum = 0.0;
            for (i = 0; i < cols; ++i) {
                sum += mCoords(m,i) * nCoords(n,i);
            }
            tmp(m,n) = sum;
            //printf("M: %d N: %d sum: %f\n", m, n, sum);
        }
    }
    scores.take(tmp);
}

void DynProg::score_covariance(MatF &mCoords, MatF &nCoords, MatF &scores) {
    // VERIFIED with Vec2D::pearsons_r subroutine!
    //Strategy is to cache the redundant values and only calculate
    //the MxN values that we must have
    int s_mlen = mCoords.rows();// s_cols = length_n
    int s_nlen = nCoords.rows();// s_rows = length_m  // Both rows and cols derived from # rows
    int cols = mCoords.cols();
    assert(cols == nCoords.cols());
    //printf("WORKING IN COVARIANCE\n");
    MatF tmp(s_mlen, s_nlen);

    double *sum_x = new double[s_nlen];
    double *sum_y = new double[s_mlen];
    int i;
    for (i = 0; i < s_nlen; ++i) {
        sum_x[i] = nCoords.sum(i);
    }
    for (i = 0; i < s_mlen; ++i) {
        sum_y[i] = mCoords.sum(i);
    }

    // CALCULATE ALL PAIR calculations
    for (int n = 0; n < s_nlen; ++n) {
        for (int m = 0; m < s_mlen; ++m) {
            tmp(m,n) = (sumOfProducts(mCoords, m, nCoords, n) -
                ((sum_x[n] * sum_y[m])/cols))/cols;
        }
    }
    delete[] sum_x;
    delete[] sum_y;
    scores.take(tmp);
}

void DynProg::score_pearsons_r(MatF &mCoords, MatF &nCoords, MatF &scores) {
    //Strategy is to cache the redundant values
    int s_nlen = nCoords.rows();// s_cols = length_n
    int s_mlen = mCoords.rows();// s_rows = length_m  // Both rows and cols derived from # rows
    int cols = mCoords.cols();
    assert(cols == nCoords.cols());
    MatF tmp(s_mlen, s_nlen);

    //printf("WORKING IN PEARSONS_R\n");
    float *bot_x = new float[s_nlen]; 
    float *bot_y = new float[s_mlen]; 
    // Sum(x^2)
    float *sum_x = new float[s_nlen]; 
    float *sum_y = new float[s_mlen]; 

    int i;
    for (i = 0; i < s_nlen; ++i) {
        sum_x[i] = nCoords.sum(i);
        //         sum(x^2)                  -    ((sum_x)^2/num_elements
        bot_x[i] = ( sumXSquared(nCoords,i) ) - ( ((sum_x[i])*sum_x[i])/cols);
    }
    for (i = 0; i < s_mlen; ++i) {
        sum_y[i] = mCoords.sum(i);
        //         sum(y^2)                  -    ((sum_y)^2/num_elements
        bot_y[i] = ( sumXSquared(mCoords,i) ) - ( ((sum_y[i])*sum_y[i])/cols );
    }

    // CALCULATE ALL PAIR calculations
    for (int n = 0; n < s_nlen; ++n) {
        for (int m = 0; m < s_mlen; ++m) {
            //        sum(X * Y) -    
            double top = sumOfProducts(mCoords, m, nCoords, n) -
                ((sum_x[n] * sum_y[m])/cols);
            //  (sum(x)      * sum(y))/num_elements
            double bot = sqrt(bot_x[n] * bot_y[m]);
            if (bot == 0) { tmp(m,n) = 0; }  // no undefined 
            else { tmp(m,n) = (float)(top/bot); }
        }
    }

    delete[] bot_x; delete[] bot_y; delete[] sum_x; delete[] sum_y;
    scores.take(tmp);
}
//***XCMS*******************************************************************
void DynProg::score_pearsons_r_opt(MatF &mCoords, MatF &nCoords, MatF &scores) {
    //Strategy is to cache the redundant values
    int s_nlen = nCoords.rows();// s_cols = length_n
    int s_mlen = mCoords.rows();// s_rows = length_m  // Both rows and cols derived from # rows
    int cols = mCoords.cols();
    assert(cols == nCoords.cols());
    MatF tmp(s_mlen, s_nlen);
    const int dfd = 10;
     
    //std::cout <<"s_nlen "<<s_nlen<<"s_mlen"<<s_mlen; 
    //printf("WORKING IN PEARSONS_R\n");
    float *bot_x = new float[s_nlen]; 
    float *bot_y = new float[s_mlen]; 
    // Sum(x^2)
    float *sum_x = new float[s_nlen]; 
    float *sum_y = new float[s_mlen]; 

    int i;
    for (i = 0; i < s_nlen; ++i) {
        sum_x[i] = nCoords.sum(i);
        //         sum(x^2)                  -    ((sum_x)^2/num_elements
        bot_x[i] = ( sumXSquared(nCoords,i) ) - ( ((sum_x[i])*sum_x[i])/cols);
    }
    for (i = 0; i < s_mlen; ++i) {
        sum_y[i] = mCoords.sum(i);
        //         sum(y^2)                  -    ((sum_y)^2/num_elements
        bot_y[i] = ( sumXSquared(mCoords,i) ) - ( ((sum_y[i])*sum_y[i])/cols );
    }
    //fill the matrix with infinity
    for (int m = 0; m < s_mlen; ++m) {
      for (int n = 0; n < s_nlen; ++n) {
       tmp(m,n) = INFINITY;
        }
    }

int diff = s_mlen-s_nlen;
    // CALCULATE REQUIRED PAIR calculations
    if (diff <= 0){
      for (int m = 0; m < s_mlen; ++m) {
     //  for (int n = s_nlen-s_nlen/dfd-m-1+diff; n < s_nlen; ++n) {
          for (int n = m-s_nlen/dfd; n < s_nlen/dfd+m-2*diff; ++n) {
          //  if(n>s_nlen-m+s_nlen/dfd||n<0||n>s_nlen) 
	    if(n<0||n>=s_nlen) 
            continue;
            
            //        sum(X * Y) -  
            double top = sumOfProducts(mCoords, m, nCoords, n) -
                ((sum_x[n] * sum_y[m])/cols);
            //  (sum(x)      * sum(y))/num_elements
            double bot = sqrt(bot_x[n] * bot_y[m]);
            if (bot == 0) { tmp(m,n) = 0; }  // no undefined 
            else { tmp(m,n) = (float)(top/bot); }
         //if((n==0)&&(m==0)){std::cout << bot << " " << top << " diff<=0" << std::endl;}
        }
      }
    }
    else
    for (int m = 0; m < s_mlen; ++m) {
       //for (int n = s_nlen-s_nlen/dfd-m-1; n < s_nlen; ++n) {
         for (int n = m-s_nlen/dfd; n < s_nlen/dfd+m+2*diff; ++n) {
          //if(n>s_nlen-m+s_nlen/dfd+diff||n<0||n>s_nlen) 
            if(n<0||n>=s_nlen)     
            continue;
            
            //        sum(X * Y) -  
            double top = sumOfProducts(mCoords, m, nCoords, n) -
                ((sum_x[n] * sum_y[m])/cols);
            //  (sum(x)      * sum(y))/num_elements
            double bot = sqrt(bot_x[n] * bot_y[m]);
            if (bot == 0) { tmp(m,n) = 0; }  // no undefined 
            else { tmp(m,n) = (float)(top/bot); }
          //if((n==0)&&(m==0)){std::cout << bot << " " << top << " diff>0" << std::endl;}
        }
      }

    delete[] bot_x; delete[] bot_y; delete[] sum_x; delete[] sum_y;
    scores.take(tmp);
    //tmp.print(0,0,1);
    //tmp.print(s_mlen-1,0,1);
    //tmp.print(0,s_nlen-1,1);
    //tmp.print(s_mlen-1,s_nlen-1,1);

}
//***End XCMS***********************************************************************

void DynProg::score_pearsons_r2(MatF &mCoords, MatF &nCoords, MatF &scores) {
    //Strategy is to cache the redundant values and only calculate
    int s_mlen = mCoords.rows();// s_cols = length_n
    int s_nlen = nCoords.rows();// s_rows = length_m  // Both rows and cols derived from # rows
    int cols = mCoords.cols();
    assert(cols == nCoords.cols());
    //printf("WORKING IN PEARSONS_R2\n");
    // CACHE all the values we can:
    float *bot_x = new float[s_nlen]; 
    float *bot_y = new float[s_mlen]; 
    // Sum(x^2)
    MatF tmp(s_mlen, s_nlen);

    MatF y_minus_mean(mCoords.rows(), mCoords.cols()); 
    MatF x_minus_mean(nCoords.rows(), nCoords.cols()); 

    int i;
    float colsF = (float)cols;
    for (i = 0; i < s_nlen; ++i) {
        float mean_x = nCoords.sum(i)/colsF;
        _subtract(nCoords, i, mean_x, x_minus_mean);
        bot_x[i] = sumXSquared(x_minus_mean, i); 
    }
    for (i = 0; i < s_mlen; ++i) {
        float mean_y = mCoords.sum(i)/colsF;
        _subtract(mCoords, i, mean_y, y_minus_mean);
        bot_y[i] = sumXSquared(y_minus_mean, i); 
    }

    // CALCULATE ALL PAIR calculations
    for (int n = 0; n < s_nlen; ++n) {
        for (int m = 0; m < s_mlen; ++m) {
            double sumOfProds = sumOfProducts(x_minus_mean, n, y_minus_mean, m);
            double top = sumOfProds * sumOfProds;
            double bot = bot_x[n] * bot_y[m];
            tmp(m,n) = top/bot;

            if (bot == 0) { tmp(m,n) = 0; }  // no undefined 
            else { tmp(m,n) = (float)(top/bot); }
        }
    }
    delete[] bot_x; delete[] bot_y;
    scores.take(tmp);
}

void DynProg::score_euclidean(MatF &mCoords, MatF &nCoords, MatF &scores) {
    int i;
    int s_mlen = mCoords.rows();// s_cols = length_n
    int s_nlen = nCoords.rows();// s_rows = length_m  // Both rows and cols derived from # rows
    int cols = mCoords.cols();
    assert(cols == nCoords.cols());
    MatF tmp(s_mlen, s_nlen);
    for (int m = 0; m < s_mlen; ++m) {
        for (int n = 0; n < s_nlen; ++n) {
            float sum = 0.0;
            for (i = 0; i < cols; ++i) {
                float diff = mCoords(m,i) - nCoords(n,i);
                sum += diff * diff;
            }
            tmp(m,n) = sqrt(sum);
            //printf("M: %d N: %d sum: %f\n", m, n, sum);
        }
    }
    scores.take(tmp);
}


void DynProg::score_mutual_info(MatF &mCoords, MatF &nCoords, MatF &scores, int num_bins) {
    //Strategy is to cache the redundant values
    int s_mlen = mCoords.rows();// s_rows = length_m  // Both rows and cols derived from # rows
    int s_nlen = nCoords.rows();// s_cols = length_n
    int cols = nCoords.cols();
    assert(cols == mCoords.cols());

    MatF tmpmat(s_mlen, s_mlen);

    // SETUP the default values:
    int MI_NUM_BINS = num_bins;
    float nCoordsMin, nCoordsMax, mCoordsMin, mCoordsMax;
    nCoords.min_max(nCoordsMin, nCoordsMax);
    mCoords.min_max(mCoordsMin, mCoordsMax);
    float MI_MAX_VAL = max(nCoordsMax, mCoordsMax); 
    float MI_MIN_VAL = min(nCoordsMin, mCoordsMin);
    float MI_SPAN = MI_MAX_VAL - MI_MIN_VAL;
    float MI_SCALE_FACTOR = MI_SPAN/MI_NUM_BINS;


    //printf("WORKING IN mutualInfo\n");
    //printf("MI_SCALE_FACTOR: %f MIMAX: %f MIN: %f\n", MI_SCALE_FACTOR, MI_MAX_VAL, MI_MIN_VAL);
    // CACHE all the values we can:
    VecF entropyX(s_nlen); 
    VecF entropyY(s_mlen); 
    MatI binIndNCoords(nCoords.rows(), nCoords.cols());
    MatI binIndMCoords(mCoords.rows(), mCoords.cols());

    assert(nCoords.cols() == mCoords.cols());

    int i;
    for (i = 0; i < nCoords.rows() ; ++i) {
        entropyX[i] = entropy(nCoords, i, MI_NUM_BINS, MI_MIN_VAL, MI_SCALE_FACTOR, binIndNCoords);
    }

    //printf("INDARR1:\n");
    //binIndNCoords.print();
    for (i = 0; i < mCoords.rows(); ++i) {
        entropyY[i] = entropy(mCoords, i, MI_NUM_BINS, MI_MIN_VAL, MI_SCALE_FACTOR, binIndMCoords);
    }
    //printf("INDARR2:\n");
    //binIndMCoords.print();

    // CALCULATE ALL PAIR calculations

    entropyXY(binIndNCoords, binIndMCoords, entropyX, entropyY, tmpmat, MI_NUM_BINS);
    scores.take(tmpmat);
}

void DynProg::score(MatF &mCoords, MatF &nCoords, MatF &scores, const char *type, int mi_num_bins) {
    if (!strcmp(type,"prd")) {
        score_product(mCoords,nCoords,scores);
    }
    else if (!strcmp(type,"cov")) {
        score_covariance(mCoords,nCoords,scores);
    }
    else if (!strcmp(type,"cor")) {
        score_pearsons_r(mCoords,nCoords,scores);
    }
//***XCMS***************************************
    else if (!strcmp(type,"cor_opt")) {
        score_pearsons_r_opt(mCoords,nCoords,scores);
    }
//***End XCMS***********************************
    else if (!strcmp(type,"euc")) {
        score_euclidean(mCoords,nCoords,scores);
    }
    // Not actually using these right now, but...
    else if (!strcmp(type,"pearsons_r2")) {
        score_pearsons_r2(mCoords,nCoords,scores);
    }
    else if (!strcmp(type,"mutual_info")) {
        score_mutual_info(mCoords,nCoords,scores,mi_num_bins);
    }
    else {
        printf("Unrecognized score type!: %s\n", type);
        exit(1);
    }
}


void DynProg::expandFlag(MatI &flagged, int flag, int numSteps, MatI &expanded) {
    int m_length = flagged.rows();
    int n_length = flagged.cols();
    MatI tmpExpanded(m_length, n_length);
    MatI definedExpanded(m_length, n_length, 0);
    int val = 0;
    int ml = 0;
    int mh = 0;
    int nl = 0;
    int nh = 0;
    for (int m = 0; m < m_length; ++m) {
        for (int n = 0; n < n_length; ++n) {
            val = flagged(m,n); 
            //printf("v: %d ", val);
            if (!definedExpanded(m,n)) {
                tmpExpanded(m,n) = val;
                definedExpanded(m,n) = 1;
            }
            if (val == flag) {
                ml = m - numSteps;
                mh = m + numSteps;
                nl = n - numSteps;
                nh = n + numSteps;
                if (ml < 0) { ml = 0; }
                if (mh >= m_length) { mh = m_length - 1; }
                if (nl < 0) { nl = 0; }
                if (nh >= n_length) { nh = n_length - 1; }
                for (int mm = ml; mm <= mh; mm++)  {
                    for (int nn = nl; nn <= nh; nn++) {
                        tmpExpanded(mm,nn) = flag;
                    }
                }
            }
        }
    }
    expanded = tmpExpanded;
} 

void entropyXY(MatI &binIndX, MatI &binIndY, VecF &entropyX, VecF &entropyY, MatF &scores, int numBins) {
    assert(binIndX.cols() == binIndY.cols());
    for (int m = 0; m < binIndY.rows(); ++m) {
        for (int n = 0; n < binIndX.rows(); ++n) {
            MatI counts(numBins, numBins,0);
            //printf("CoUNTs:\n");
            //counts.print();
			int i;
            for (i = 0; i < binIndX.cols(); ++i) {
                counts(binIndY(m,i),binIndX(n,i))++;
            }
            //printf("CoUNTs (after):\n");
            //counts.print();
            float totxy = (float)binIndY.cols();
            //printf("TOTXY: %f\n", totxy);
            float entropyXY = 0.f;
            float probxy;
            for (i = 0; i < numBins; ++i) {
                for (int j = 0; j < numBins; j++) {
                    probxy = ((float)(counts(i,j)))/totxy;
                    if (probxy != 0.f) {
                        entropyXY -= probxy * logf(probxy)/_LOG2;
                    }
                    //printf("PROBXY: %f\n", probxy);
                    //printf("logPROBXY: %f\n", log2(probxy));
                    //printf("entropy(ONGOING): %f\n",entropyXY); 

                }
            }
            //printf("M: %d N: %d\n", m,n);
            //printf("lYm: %d lxn: %d\n", entropyY.dim(),entropyX.dim());
            //printf("******************************\n");
            scores(m,n) = entropyY[m] + entropyX[n] - entropyXY;
            //printf("[ENTRY[%d]]%f + [ENTX[%d]]%f  - [entropyXY]%f\n = %f", m, entropyY[m], n, entropyX[n], entropyXY, scores(m,n));
        }
    }
}

// Calculate the entropy of a vector (scaleFactor is the span/numBins;
float entropy(MatF &mat, int rowNum, int numBins, float minVal, float scaleFactor, MatI &indArray) {
    VecI binArray(numBins, 0);
    //binArray = 0;  // zero it!
	int i;
    for (i = 0; i < mat.cols(); ++i) {
        int ind = (int)( ((float)mat(rowNum,i) - minVal)/scaleFactor );
        if (ind == numBins) {
            ind = numBins - 1;
        }
        binArray[ind]++;
        indArray(rowNum,i) = ind;
    }
    int tot = mat.cols();
    float probx;
    float entropy = 0;
    for (i = 0; i < numBins; ++i) {
        probx = ((float)binArray[i])/tot;
        //printf("PROBX: %f\n", probx);
        //printf("logPROBX: %f\n", log2(probx));
        if (probx != 0) {
            entropy -= probx * logf(probx)/_LOG2;
        }
    }
    //printf("ENTROPY: %f\n", entropy);
    return entropy;
}

//Subtract a value
void _subtract(MatF &mat, int rowNum, float val, MatF &minused) {
    float *matptr = mat.pointer(rowNum);
    float *minusedptr = minused.pointer(rowNum);
    for (int i = 0; i < mat.cols(); ++i) {
        minusedptr[i] = matptr[i] - val;
    }
}

//Sum of the products (i.e. the dot product at that row)
float sumOfProducts(MatF &mat1, int rowNum1, MatF &mat2, int rowNum2) {
    float *mat1ptr = mat1.pointer(rowNum1);
    float *mat2ptr = mat2.pointer(rowNum2);
    float sum = 0;
    for (int i = 0; i < mat1.cols(); ++i) {
        sum += mat1ptr[i] * mat2ptr[i]; 
    }
    return sum;
}

// Returns the sum of the square of the values in the row number
// could increase the speed here by getting the oneD version and doing pointer
// math(?)
float sumXSquared(MatF &mat, int rowNum) {
    float *matptr = mat.pointer(rowNum);
    float sum = 0;
    for (int i = 0; i < mat.cols(); ++i) {
        sum += matptr[i] * matptr[i]; 
    }
    return sum;
}

void DynProg::default_gap_penalty(MatF &smat, VecF &out) {
    int _length = smat.cols() + smat.rows(); 
    float _avg = smat.avg();
    linear_less_before(DEFAULT_GAP_PENALTY_SLOPE, _avg, _length, out);
}

void DynProg::less_before(VecF &arr) {
    for (int i = arr.len() - 1; i >= 1; --i) {
        arr[i] -= arr[i-1];
    }
}

// linear function given in terms of mx + b where m is slope and b is y
// intercept
void DynProg::linear(float m, float b, int len,  VecF &arr) {
    float *tmparr = new float[len];
    int ind = 0;
    while (ind < len) { 
        tmparr[ind] = m*ind + b;
        ++ind;
    }
    arr.take(len, tmparr);
}

void DynProg::linear_less_before(float m, float b, int veclen, VecF &lessbefore) {
    float *tmparr = new float[veclen];
    tmparr[0] = b;
    for (int i = 1; i < veclen; ++i) {
        tmparr[i] = m;
    }
    lessbefore.take(veclen, tmparr);
}

// gap penalty is ZERO indexed (i.e. the _first_ gap penalty is
// accessed at gap_penalty[0]
void DynProg::find_path(MatF& smat, VecF &gap_penalty, int minimize, float diag_factor, float gap_factor, int local, float init_penalty) {

    if (gap_penalty.len() == 0) {
        default_gap_penalty(smat, gap_penalty);
        if (minimize) { 
            gap_penalty *= -1.f;
        }
    }

    // Initialize matrices:
    int rows = smat.rows();
    int cols = smat.cols();
    MatF tmp_asmat(rows, cols);
    MatI tmp_tb(rows, cols);
    MatI tmp_tbpath(rows, cols,0);
    MatI tmp_gapmat(rows, cols);
    _smat = &smat; // save a pointer to smat

    // ********************************************************
    // * BEGIN CALC ADDITIVE SCORE MATRIX
    // ********************************************************
    int length_m = rows;
    int length_n = cols;

    //smat.print();
    // Initialize top left cell:
    tmp_asmat(0,0)  = smat(0,0);
    tmp_gapmat(0,0) = 0;
    tmp_tb(0,0)     = 0;

    // GLOBAL:
    if (!local) { 
            // FILL in the left and top sides with gaps!
            for (int m = 1; m < length_m; ++m) {
                tmp_asmat(m,0) =  (smat(m,0) * gap_factor) + tmp_asmat(m-1,0) - gap_penalty[tmp_gapmat(m-1,0)];
                tmp_gapmat(m,0) = m;
                tmp_tb(m,0) = 1;  // path is from above 
            }
            for (int n = 1; n < length_n; ++n) {
                tmp_asmat(0,n)= (smat(0,n) * gap_factor) + tmp_asmat(0,n - 1) - gap_penalty[tmp_gapmat(0,n - 1)];
                tmp_gapmat(0,n) = n;
                tmp_tb(0,n) = 2;  // path is from left 
            }
    }
    // LOCAL:
    else {   
        if (minimize) {
            for (int m = 1; m < length_m; ++m) {
                float diag = (smat(m,0) * diag_factor ) - init_penalty;  // drop in from left side
                float top  = (smat(m,0) * gap_factor) + tmp_asmat(m - 1,0) - gap_penalty[(tmp_gapmat(m - 1,0))];
                float best_val;
                int best_pos;

                // If minimize:
                if (diag <= top) { best_pos = 0; best_val = diag; }
                else { best_pos = 1; best_val = top; }

                // Set the gap matrix:
                if (best_pos == 1) { tmp_gapmat(m,0) = tmp_gapmat(m - 1,0) + 1; }
                else { tmp_gapmat(m,0) = 0; }  
                tmp_asmat(m,0) = best_val; tmp_tb(m,0)  = best_pos;
            }
            for (int n = 1; n < length_n; ++n) {
                float diag = (smat(0,n) * diag_factor) - init_penalty;  // drop in from top
                float left = (smat(0,n) * gap_factor) + tmp_asmat(0,n - 1) - gap_penalty[(tmp_gapmat(0,n -1))];
                float best_val;
                int best_pos;

                if (diag <= left) { best_pos = 0; best_val = diag; }
                else { best_pos = 2; best_val = left; }

                // Set the gap matrix:
                if (best_pos == 2) { tmp_gapmat(0,n) = tmp_gapmat(0,n - 1) + 1; }
                else { tmp_gapmat(0,n) = 0; }
                tmp_asmat(0,n) = best_val; tmp_tb(0,n)  = best_pos;
            }
        }
        else {  // maximize:

            for (int m = 1; m < length_m; ++m) {
                float diag = (smat(m,0) * diag_factor) - init_penalty;  // drop in from left side
                float top  = (smat(m,0) * gap_factor) + tmp_asmat(m - 1,0) - gap_penalty[(tmp_gapmat(m - 1,0))];
                float best_val;
                int best_pos;

                // If !minimize:
                if (diag >= top) { best_pos = 0; best_val = diag; }
                else { best_pos = 1; best_val = top; }

                // Set the gap matrix:
                if (best_pos == 1) { tmp_gapmat(m,0) = tmp_gapmat(m - 1,0) + 1; }
                else { tmp_gapmat(m,0) = 0; }  
                tmp_asmat(m,0) = best_val; tmp_tb(m,0)  = best_pos;
            }
            for (int n = 1; n < length_n; ++n) {
                float diag = (smat(0,n) * diag_factor) - init_penalty;  // drop in from top
                float left = (smat(0,n) * gap_factor) + tmp_asmat(0,n - 1) - gap_penalty[(tmp_gapmat(0,n -1))];
                float best_val;
                int best_pos;

                if (diag >= left) { best_pos = 0; best_val = diag; }
                else { best_pos = 2; best_val = left; }

                // Set the gap matrix:
                if (best_pos == 2) { tmp_gapmat(0,n) = tmp_gapmat(0,n - 1) + 1; }
                else { tmp_gapmat(0,n) = 0; }
                tmp_asmat(0,n) = best_val; tmp_tb(0,n)  = best_pos;
            }
        }
    }
    // COMPLETE the tmp_asmat:
    if (minimize) {  // complete the asmat with MINIMIZE!
        for (int m = 1; m < length_m; ++m) {
            for (int n = 1; n < length_n; ++n) {
                float best_val; int best_pos;
                float smat_at_ind = smat(m,n);
                // MINIMIZE: 
                float smat_at_ind_times_gap_factor = smat_at_ind * gap_factor;
                //DynProg::_min(diag, top, left, best_val, best_pos);
                DynProg::_min(
                        (smat_at_ind * diag_factor) + tmp_asmat(m-1,n-1),
                        smat_at_ind_times_gap_factor + tmp_asmat(m-1,n) - gap_penalty[tmp_gapmat(m-1,n)],
                        smat_at_ind_times_gap_factor + tmp_asmat(m,n-1) - gap_penalty[tmp_gapmat(m,n-1)],
                        best_val, best_pos
                        );
                // SET the gap_length_matrix
                if (best_pos == 1) { tmp_gapmat(m,n) = tmp_gapmat(m-1,n) + 1; }
                else if (best_pos == 2) { tmp_gapmat(m,n) = tmp_gapmat(m,n-1) + 1; }
                else { tmp_gapmat(m,n) = 0; }
                tmp_tb(m,n) = best_pos; tmp_asmat(m,n) = best_val;
            }
        }
    }
    else { // complete the asmat with MAXIMIZE!
        for (int m = 1; m < length_m; ++m) {
            for (int n = 1; n < length_n; ++n) {
                float best_val; int best_pos;
                float smat_at_ind = smat(m,n);
                // MAXIMIZE: 
                float smat_at_ind_times_gap_factor = smat_at_ind * gap_factor;
                //DynProg::_max(diag, top, left, best_val, best_pos);
                DynProg::_max(
                        ( smat_at_ind * diag_factor ) + tmp_asmat(m-1,n-1),
                        smat_at_ind_times_gap_factor + tmp_asmat(m-1,n) - gap_penalty[tmp_gapmat(m-1,n)],
                        smat_at_ind_times_gap_factor + tmp_asmat(m,n-1) - gap_penalty[tmp_gapmat(m,n-1)],
                        best_val, best_pos
                        );
                // SET the gap_length_matrix
                if (best_pos == 1) { tmp_gapmat(m,n) = tmp_gapmat(m-1,n) + 1; }
                else if (best_pos == 2) { tmp_gapmat(m,n) = tmp_gapmat(m,n-1) + 1; }
                else { tmp_gapmat(m,n) = 0; }
                tmp_tb(m,n) = best_pos; tmp_asmat(m,n) = best_val;
            }
        }
    }
    //  ************************************************************
    //  * END CALC ADD SCORE MATRIX
    //  ************************************************************ 

    float bestscore;
    int optimal_m;
    int optimal_n;
    if (local) {
        if (minimize) {
            bestscore = DynProg::_global_min(tmp_asmat, optimal_m, optimal_n);
        }
        else {
            bestscore = DynProg::_global_max(tmp_asmat, optimal_m, optimal_n);
        }
    }
    else {  //global
        bestscore = tmp_asmat(rows - 1, cols - 1);
        optimal_m = rows - 1;
        optimal_n = cols - 1;
    }

    // Reverse back the gap penalty
    if (minimize) { 
        gap_penalty *= -1.f;
    }

    _traceback(tmp_tb, smat, optimal_m, optimal_n, tmp_tbpath, _mCoords, _nCoords, _sCoords); 
    int _equivLastInd = _mCoords.dim()-1;
    _bestScore = tmp_asmat(_mCoords[_equivLastInd],_nCoords[_equivLastInd]);

    _asmat.take(tmp_asmat);
    _tb.take(tmp_tb);
    _tbpath.take(tmp_tbpath);
    _gapmat.take(tmp_gapmat);
}
