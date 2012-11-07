
#include <cstdlib>
#include "string.h"
#include "stdio.h"

#include "lmat.h"
#include "assert.h"
#include "math.h"
#include "vec.h"
#include "mat.h"

#include <R.h>

extern "C" {

  bool DEBUG = 0;
  using namespace VEC;

  LMat::LMat() : _mz_vals(0), _tm_vals(0) {
    _mz = new VecF();
    _tm = new VecF();
    _mat = new MatF();
  }

  LMat::~LMat() {
    delete _mz;
    delete _tm;
    delete _mat;
  }

  void LMat::set_from_xcms(int valuescantime, double *pscantime, int mzrange, double *mz, double *intensity) {
    delete _mz;
    delete _tm;
    delete _mat;

    // Get the time values:
    _tm_vals = valuescantime;

    float *tm_tmp = new float[_tm_vals];
    for(int i=0; i < _tm_vals; i++) {
      tm_tmp[i] = pscantime[i];
    }

    _tm = new VecF(_tm_vals, tm_tmp);

    // Get the mz values:
    _mz_vals = mzrange;
    float *mz_tmp = new float[_mz_vals];
    for(int i=0; i < _mz_vals; i++) {
      mz_tmp[i] = mz[i];
    }
    _mz = new VecF(_mz_vals, mz_tmp);

    // Read the matrix:
    int rows_by_cols = _tm_vals * _mz_vals;
    float *mat_tmp = new float[rows_by_cols];

    for(int i=0; i < rows_by_cols; i++) {
      mat_tmp[i] = intensity[i];
    }

    _mat = new MatF(_tm_vals, _mz_vals, mat_tmp);
  }

  void LMat::print_xcms() {
    float *mztmp = (float*)(*_mz);
    float *tmtmp = (float*)(*_tm);
    float *mattmp = (float*)(*_mat);
    int i;

    // The TIME vals:
    Rprintf("%d\n", _tm_vals);  // num of vals
    for (i = 0; i < _tm_vals - 1; ++i) {
      Rprintf("%f ",  tmtmp[i]);
    }
    Rprintf("%f\n",  tmtmp[i]);  // the last one

    // The M/Z vals:
    Rprintf("%d\n", _mz_vals);  // num of vals
    for (i = 0; i < _mz_vals - 1; ++i) {
      Rprintf("%f ",  mztmp[i]);
    }
    Rprintf("%f\n",  mztmp[i]);  // the last one
    for (int m = 0; m < _tm_vals; ++m) {
      int n;
      for (n = 0; n < _mz_vals - 1; ++n) {
	Rprintf("%f ", mattmp[m*_mz_vals+n]);
      }
      Rprintf("%f\n", mattmp[m*_mz_vals+n]);
    }

  }


  void LMat::mz_axis_vals(VecI &mzCoords, VecF &mzVals) {
    VecF tmp(mzCoords.length());
    for (int i = 0; i < mzCoords.length(); ++i) {
      if (mzCoords[i] < _mz_vals) {
	tmp[i] = (*_mz)[mzCoords[i]];
      }
      else {
	Rprintf("asking for mz value at index: %d (length: %d)\n", mzCoords[i], _mz_vals);
	R_ShowMessage("Serious error in obiwarp.");
      }
    }
    mzVals.take(tmp);
  }

  void LMat::tm_axis_vals(VecI &tmCoords, VecF &tmVals) {
    // Rprintf("tmCoords"); tmCoords.print();
    VecF tmp(tmCoords.length());
    //printf("tm_vals %d \n", _tm_vals);
    for (int i = 0; i < tmCoords.length(); ++i) {
      if (tmCoords[i] < _tm_vals) {
	tmp[i] = (*_tm)[tmCoords[i]];
	//printf("tmCoords[i] %d val out %f\n", tmCoords[i], tmp[i]);
      }
      else {
	Rprintf("asking for time value at index: %d (length: %d)\n", tmCoords[i], _tm_vals);
	R_ShowMessage("Serious error in obiwarp.");
      }
    }
    tmVals.take(tmp);
  }


  void LMat::chomp_plus_spaces( char *str ) {
    if( str ) {
      int len = strlen( str );
      if ( len <= 0 ) return;
      while ( --len ) {
	if ( str[len]=='\r' || str[len]=='\n' ) {
	  str[len] = 0;
	}
	else break;
      }
      // At this point len == strlen(str) - 1
      len = len+1;
      while ( --len ) {
	if ( str[len] != ' ' ) {
	  break;
	}
	else {
	  str[len] = 0;
	}
      }
    }
  }


  void LMat::warp_tm(VecF &selfTimes, VecF &equivTimes) {
    VecF out;
    VecF::chfe(selfTimes, equivTimes, *_tm, out, 1);  // run with sort option
    _tm->take(out);
  }

}
