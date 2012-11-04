
#ifndef _LMAT_H
#define _LMAT_H

#include "vec.h"
#include "mat.h"

extern "C"{
using namespace VEC;

class LMat {
    private:
#define LEN_LARGEST_NUM (30)
#define LARGEST_NUM_MZ_VALS (40000)
#define LARGEST_NUM_TIME_VALS (40000)
    public:
        int _mz_vals;
        int _tm_vals;

        // All constructors call new!
        // All swaps of these MUST delete their memory before swapping!
        MatF *_mat;
        VecF *_mz;
        VecF *_tm;

        LMat();
        // Takes a binary lmat file as input
        LMat(const char *file);

        ~LMat();
        int mzlen() { return _mz_vals; }
        int tmlen() { return _tm_vals; }
        int num_mz() { return _mz_vals; }
        int num_tm() { return _tm_vals; }
        MatF * mat() { return _mat; }
        VecF * mz() { return _mz; }
        VecF * tm() { return _tm; }

        float hi_mz() { return (*_mz)[_mz_vals-1]; }
        float lo_mz() { return (*_mz)[0]; }
        float hi_tm() { return (*_tm)[_tm_vals-1]; }
        float lo_tm() { return (*_tm)[0]; }
        void mz_axis_vals(VecI &mzCoords, VecF &mzVals);
        void tm_axis_vals(VecI &tmCoords, VecF &tmVals);

        void set_from_xcms(int valuescantime, double *pscantime, int mzrange,
			   double *mz, double *intensity);
        void print_xcms();

        // selfTimes and equivTimes are the anchor points for the warping
        // function..  warps the time values (not the actual data values)
        void warp_tm(VecF &selfTimes, VecF &equivTimes);

        // expects one line with the # mz vals and next with the vals
        void set_mz_from_ascii(FILE *fpt);
        // expects one line with the # tm vals and next with the vals
        void set_tm_from_ascii(FILE *fpt);
        // expects the matrix in ascii format
        void set_mat_from_ascii(FILE *ptr, int rows, int cols);
        // writes the lmat in binary to a file (or STDOUT if NULL)
        void write(const char *file=NULL);
        // writes the lmat in ascii to a file (or STDOUT if NULL)
        void print(const char *file=NULL);

        // obviously not the final resting place
        void chomp_plus_spaces( char *str);
};

#endif

}
