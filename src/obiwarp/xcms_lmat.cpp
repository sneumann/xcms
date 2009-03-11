
#include <cstdlib>
#include "string.h"
#include "stdio.h"

#include "lmat.h"
#include "assert.h"
#include "math.h"
#include "vec.h"
#include "mat.h"


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



void LMat::set_from_binary_mat(const char *file) {
    delete _mz;
    delete _tm;
    delete _mat;
    _mat = new MatF;
    _mat->set_from_binary(file);
    _tm_vals = _mat->mlen();
    _mz_vals = _mat->nlen();
    _tm = new VecF(_tm_vals);
    _mz = new VecF(_mz_vals);
    int i;
    for (i = 0; i < _tm->size(); ++i) {
        _tm->at(i) = (float)i;
    }
    for (i = 0; i < _mz->size(); ++i) {
        _mz->at(i) = (float)i;
    }
}

void LMat::set_from_ascii_mat(const char *file) {
    delete _mz;
    delete _tm;
    delete _mat;
    _mat = new MatF;
    _mat->set_from_ascii(file);
    _tm_vals = _mat->mlen();
    _mz_vals = _mat->nlen();
    _tm = new VecF(_tm_vals);
    _mz = new VecF(_mz_vals);
    int i;
    for (i = 0; i < _tm->size(); ++i) {
        _tm->at(i) = (float)i;
    }
    for (i = 0; i < _mz->size(); ++i) {
        _mz->at(i) = (float)i;
    }
}

void LMat::set_from_binary(const char *file) {
    delete _mz;
    delete _tm;
    delete _mat;
    FILE *fh = fopen(file, "rb"); 
    if (fh == NULL) {
        printf("Could not open %s for reading\n", file);
        exit(1);
    }
    // Get the time values:
    fread(&_tm_vals, sizeof(int), 1, fh);
    float *tm_tmp = new float[_tm_vals];
    fread(tm_tmp, sizeof(float), _tm_vals, fh);
    _tm = new VecF(_tm_vals, tm_tmp);

    // Get the mz values:
    fread(&_mz_vals, sizeof(int), 1, fh);
    float *mz_tmp = new float[_mz_vals];
    fread(mz_tmp, sizeof(float), _mz_vals, fh);
    _mz = new VecF(_mz_vals, mz_tmp);

    // Read the matrix:
    int rows_by_cols = _tm_vals * _mz_vals;
    //printf("rbycools: %d\n", rows_by_cols);
    float *mat_tmp = new float[rows_by_cols];
    fread(mat_tmp, sizeof(float), rows_by_cols, fh);
    //puts("**********************************************************");
    //puts("THIS IS THE BINARY MAT READ in:");
    //printf("First: %d:%.0f\n", 0, mat_tmp[0]);
    //printf("Last: %d:%.0f\n", rows_by_cols, mat_tmp[rows_by_cols-1]);
    //puts("**********************************************************");
    _mat = new MatF(_tm_vals, _mz_vals, mat_tmp);
    fclose(fh);
}


//xcms----------------------------------------------------
//--------------------------------------------------------
void LMat::set_from_xcms(int valuescantime, double *pscantime, int mzrange, double *mz, double *intensity) {
    delete _mz;
    delete _tm;
    delete _mat;
    
    // Get the time values:
    _tm_vals = valuescantime;
//printf("hier:%d\n",_tm_vals);
    float *tm_tmp = new float[_tm_vals];
    for(int i=0; i < _tm_vals; i++)
       tm_tmp[i] = pscantime[i];
//printf("tm_tmp[i]:%f\n",tm_tmp[i]);}   
//   float tm_tmp = 42.42;
    _tm = new VecF(_tm_vals, tm_tmp);

    // Get the mz values:
    _mz_vals = mzrange;
    float *mz_tmp = new float[_mz_vals];
    for(int i=0; i < _mz_vals; i++)
       mz_tmp[i] = mz[i];  
    _mz = new VecF(_mz_vals, mz_tmp);

    // Read the matrix:
    int rows_by_cols = _tm_vals * _mz_vals;
    //printf("rbycools: %d\n", rows_by_cols);
    float *mat_tmp = new float[rows_by_cols];
    for(int i=0; i < rows_by_cols; i++)
       mat_tmp[i] = intensity[i];
    //puts("**********************************************************");
    //puts("THIS IS THE BINARY MAT READ in:");
    //printf("First: %d:%.0f\n", 0, mat_tmp[0]);
    //printf("Last: %d:%.0f\n", rows_by_cols, mat_tmp[rows_by_cols-1]);
    //puts("**********************************************************");
    _mat = new MatF(_tm_vals, _mz_vals, mat_tmp);
}

//--------------------------------------------------------
//end xcms------------------------------------------------
 
void LMat::write(const char *file) {
    FILE *fh;
    if (file != NULL) { fh = fopen(file, "wb"); }
    else { fh = stdout; }

    int sizeof_float = sizeof(float);
    int rows_by_cols = _mat->rows() * _mat->cols();
    fwrite(&_tm_vals, sizeof(int), 1, fh);
    fwrite((float*)(*_tm), sizeof_float, _tm_vals, fh);
    fwrite(&_mz_vals, sizeof(int), 1, fh);
    fwrite((float*)(*_mz), sizeof_float, _mz_vals, fh);
    fwrite((float*)(*_mat), sizeof_float, rows_by_cols, fh);

    if (file != NULL) { fclose(fh); }
}

void LMat::print(const char *file) {
    float *mztmp = (float*)(*_mz);
    float *tmtmp = (float*)(*_tm);
    float *mattmp = (float*)(*_mat);
    int i;

    FILE*fpt;
    if (file == NULL) { fpt = stdout; }
    else { fpt = fopen(file, "w"); }

    // The TIME vals:
    fprintf(fpt, "%d\n", _tm_vals);  // num of vals
    for (i = 0; i < _tm_vals - 1; ++i) {
        fprintf(fpt, "%f ",  tmtmp[i]);
    }
    fprintf(fpt, "%f\n",  tmtmp[i]);  // the last one

    // The M/Z vals:
    fprintf(fpt, "%d\n", _mz_vals);  // num of vals
    for (i = 0; i < _mz_vals - 1; ++i) {
        fprintf(fpt, "%f ",  mztmp[i]);
    }
    fprintf(fpt, "%f\n",  mztmp[i]);  // the last one
    for (int m = 0; m < _tm_vals; ++m) {
        int n;
        for (n = 0; n < _mz_vals - 1; ++n) {
            fprintf(fpt, "%f ", mattmp[m*_mz_vals+n]);
        }
        fprintf(fpt, "%f\n", mattmp[m*_mz_vals+n]);
    }

    if (file != NULL) { fclose(fpt); }
}
//XCMS*************************************************
void LMat::print_xcms() {
    float *mztmp = (float*)(*_mz);
    float *tmtmp = (float*)(*_tm);
    float *mattmp = (float*)(*_mat);
    int i;

     // The TIME vals:
    printf("%d\n", _tm_vals);  // num of vals
    for (i = 0; i < _tm_vals - 1; ++i) {
        printf("%f ",  tmtmp[i]);
    }
    printf("%f\n",  tmtmp[i]);  // the last one

    // The M/Z vals:
    printf("%d\n", _mz_vals);  // num of vals
    for (i = 0; i < _mz_vals - 1; ++i) {
        printf("%f ",  mztmp[i]);
    }
    printf("%f\n",  mztmp[i]);  // the last one
    for (int m = 0; m < _tm_vals; ++m) {
        int n;
        for (n = 0; n < _mz_vals - 1; ++n) {
            printf("%f ", mattmp[m*_mz_vals+n]);
        }
        printf("%f\n", mattmp[m*_mz_vals+n]);
    }

}
/**************
float LMat::back() {
    
    float *tmtmp = (float*)(*_tm);
    int i;

    // The TIME vals:
    //for (i = 0; i < _tm_vals - 1; ++i) {
    //    fprintf(fpt, "%f ",  tmtmp[i]);
    //}
    return tmtmp; 

    
}
//END XCMS *********************************************
*************/
void LMat::set_from_ascii(const char *file) {
    FILE *fpt;
    fpt = fopen(file,"rt");
    assert(fpt);
    set_tm_from_ascii(fpt);
    set_mz_from_ascii(fpt);
    set_mat_from_ascii(fpt,_tm->dim(), _mz->dim());
    fclose(fpt);
}


void LMat::set_mat_from_ascii(FILE *fpt, int rows, int cols) {
    MatF _mat_tmp(rows,cols);

    char *line = new char[cols*LEN_LARGEST_NUM];

    int row = 0;
    int col = 0;
    while( fgets(line,cols*LEN_LARGEST_NUM,fpt) ) {
        col = 0;
        chomp_plus_spaces(line);
        char *wordStart = line;
        for( char *c=line; *c; c++ ) {
            if( *c == ' ') {
                *c = 0;
                _mat_tmp(row,col) = (float)atof( wordStart );
                wordStart = c+1;
                col++;
            }
        }
        _mat_tmp(row,col) = (float)atof(wordStart);
        col++;
        row++;
    }
    assert(cols == col);
    assert(rows == row);

    _mat->take(_mat_tmp);

    delete [] line;
    line = 0;
}

void LMat::mz_axis_vals(VecI &mzCoords, VecF &mzVals) {
    VecF tmp(mzCoords.length());
    for (int i = 0; i < mzCoords.length(); ++i) {
         if (mzCoords[i] < _mz_vals) {
            tmp[i] = (*_mz)[mzCoords[i]];
        }
        else {
            printf("asking for mz value at index: %d (length: %d)\n", mzCoords[i], _mz_vals);
            exit(1);
        }
    }
    mzVals.take(tmp);
}

void LMat::tm_axis_vals(VecI &tmCoords, VecF &tmVals) {
   // puts("tmCoords"); tmCoords.print();
     VecF tmp(tmCoords.length());
    //printf("tm_vals %d \n", _tm_vals);
    for (int i = 0; i < tmCoords.length(); ++i) {
        if (tmCoords[i] < _tm_vals) {
            tmp[i] = (*_tm)[tmCoords[i]];
            //printf("tmCoords[i] %d val out %f\n", tmCoords[i], tmp[i]); 
        }
        else {
            printf("asking for time value at index: %d (length: %d)\n", tmCoords[i], _tm_vals);
            exit(1);
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

 
// File ptr should be on the line indicating the number of time points
// Next line holds the time values
void LMat::set_tm_from_ascii(FILE *fpt) {
    // Get number of m (time) vals:
    char *shortline = new char[LEN_LARGEST_NUM];
    fgets(shortline, LEN_LARGEST_NUM, fpt);
    int num_m_vals = atoi(shortline);
    delete[] shortline;
    
    float *_tmarr = new float[num_m_vals];

    char *line = new char[num_m_vals*LEN_LARGEST_NUM];
    fgets(line, num_m_vals*LEN_LARGEST_NUM, fpt);

    int ind = 0;
    char *wordStart = line;
    chomp_plus_spaces(line);
    for( char *c=line; *c; ++c ) {
        if( *c == ' ') {
            *c = 0;  // turn the space at the end of a into the end of string
            _tmarr[ind] = (float)atof( wordStart );
            wordStart = c+1;  // a new word begins!
            ind++;
        }
    }
    _tmarr[ind] = (float)atof( wordStart );  // grab the last one
    _tm_vals = num_m_vals;
    delete[] line;
    delete _tm;
    _tm = new VecF(num_m_vals, _tmarr);
}


// File ptr should be on the line indicating the number of mz points
// Next line holds the mz values
void LMat::set_mz_from_ascii(FILE *fpt) {
    // Get number of n (mz) vals:
    char *shortline = new char[LEN_LARGEST_NUM];
    fgets(shortline, LEN_LARGEST_NUM, fpt);
    int num_n_vals = atoi(shortline);
    delete[] shortline;
    
    float *_mzarr = new float[num_n_vals];

    char *line = new char[num_n_vals*LEN_LARGEST_NUM];
    fgets(line, num_n_vals*LEN_LARGEST_NUM, fpt);

    int ind = 0;
    char *wordStart = line;
    chomp_plus_spaces(line);
    for( char *c=line; *c; ++c ) {
        if( *c == ' ') {
            *c = 0;  // turn the space at the end of a into the end of string
            _mzarr[ind] = (float)atof( wordStart );
            wordStart = c+1;  // a new word begins!
            ind++;
        }
    }
    _mzarr[ind] = (float)atof( wordStart );  // grab the last one
    _mz_vals = num_n_vals;
    delete[] line;
    delete _mz;
    _mz = new VecF(num_n_vals, _mzarr);
}


void LMat::warp_tm(VecF &selfTimes, VecF &equivTimes) {
    VecF out;
    VecF::chfe(selfTimes, equivTimes, *_tm, out, 1);  // run with sort option
    _tm->take(out);
}



//int LMat::set_tm_from_ascii(FILE *fpt) {
//    int max_chars_line = LARGEST_NUM_MZ_VALS*LEN_LARGEST_NUM;
//    char *line = new char[max_chars_line];
//    fgets(line, max_chars_line, fpt);
//    chomp_plus_spaces(line);
//    int cnt = 0;
//	char *c;
//    for (c=line; *c; c++) {
//        if ( *c == ' ') {
//            cnt++; 
//        }
//    }
//    cnt++;
//    float *_tm_arr = new float[cnt];
//    char *wordStart = line;
//    int ind = 0;
//    for (c=line; *c; c++) {
//        if( *c == ' ' || *c == 0 ) {
//            *c = 0;
//            _tm_arr[ind] = (float)atof( wordStart );
//            wordStart = c+1;
//            ind++;
//        }
//    }
//    _tm_arr[ind] = (float)atof( wordStart );
//
//    // Set our object
//    delete _tm;  // delete the current one (no mem leaks!)
//    _tm = new VecF(cnt, _tm_arr);
//    _tm_vals = cnt;
//
//    delete[] line;
//    line = 0;
//	return 1;
//}


}
