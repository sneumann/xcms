#ifndef _MAT_H
#define _MAT_H

#include "vec.h"

/*************************************************************
 * Creation from existing object/array is always shallow!.
 * Will delete any memory allocated.
 * Will NOT delete any memory not allocated.
 * If you want deep then use copy function!
 ************************************************************/

namespace VEC {

class MatI;
class MatF;
class MatD;

// BEGIN TEMPLATE

class MatI {

    public:
        // length
        int _m;
        int _n;
        VecI _dat;
        // Constructors:
        MatI();
        MatI(int m, int n);
        MatI(int m, int n, const int &val);


        // (copied from vec.h)
        // if (shallow == 1 (true)) then no memory is deleted upon destruction
        // if (shallow == 0 (false)) then delete[] is called
        // FOR THIS CONSTRUCTOR ONLY, there is no DEEP copying, EVER!
        MatI(int m, int n, int *arr, bool shallow=0);

        // (copied from vec.h)
        // if (shallow == 0 (false)) a DEEP copy is made of the data
        // if (shallow == 1 (true)) a copy of the pointer is made
        // if (shallow) then no memory is released upon destruction
        // shallow is used for a quick copy with which to work
        MatI(const MatI &A, bool shallow=0);

        operator int*() { return (int*)_dat; }
        operator const int*() { return (int*)_dat; }
        int* pointer() { return (int*)_dat; }
        int* pointer(int m) { return &_dat[m*_n]; }
        // creates vec objects
        // caller must have allocated the array for the vec objects
        // the data is a shallow copy!
        // transpose and call row_vecs for col_vecs!
        void row_vecs(int &cnt, VecI *vecs);

        MatI & operator=(const int &val);
        // DEEP
        MatI & operator=(MatI &A);
        ~MatI();
        // Deep copy unless shallow == true
        void copy(MatI &receiver, bool shallow=0) const;

        void set_from_ascii(std::ifstream &stream, int m, int n, MatI &out);
        void set_from_ascii(std::ifstream &stream, MatI &out);
        void set_from_ascii(const char *file, bool without_axes=0);
        void set_from_binary(const char *file);
        void file_rows_cols(std::ifstream &stream, int &rows, int &cols);
        // tnt_array2d_utils.h has a good example (use ifstream)

        // shallow copy and no ownership of memory
        void set(int m, int n, int *arr);
        // shallow copy and no ownership of memory
        void set(MatI &A);

        bool all_equal() {
            return _dat.all_equal();
        }

        // Deletes the object's memory (if not shallow) and takes ownership
        // of the array memory (we will call delete[])
        void take(int m, int n, int *arr);
        // Deletes previous memory (if not shallow) and takes ownership
        // of the other's memory.
        void take(MatI &A);

        // flattens the matrix and returns a vector
        void to_vec(VecI &outvec, bool shallow=0);

        bool operator==(const MatI &A);

        bool shallow() { return _dat.shallow(); }
        int dim1() const { return _m; }
        int dim2() const { return _n; }
        int mlen() const { return _m; }
        int nlen() const { return _n; }
        int rows() const { return _m; }
        int cols() const { return _n; }

        int& operator()(int m, int n) {
#ifdef JTP_BOUNDS_CHECK
            if (n < 0) { Rprintf("n < 0"); R_ShowMessage("Serious error in obiwarp.");
}
            if (n >= _n) { Rprintf("n >= _n"); R_ShowMessage("Serious error in obiwarp.");
}
            if (m < 0) { Rprintf("m < 0"); R_ShowMessage("Serious error in obiwarp.");
}
            if (m >= _m) { Rprintf("m >= _m"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[m*_n + n];
        }
        const int& operator()(int m, int n) const {
#ifdef JTP_BOUNDS_CHECK
            if (n < 0) { Rprintf("n < 0"); R_ShowMessage("Serious error in obiwarp.");
}
            if (n >= _n) { Rprintf("n >= _n"); R_ShowMessage("Serious error in obiwarp.");
}
            if (m < 0) { Rprintf("m < 0"); R_ShowMessage("Serious error in obiwarp.");
}
            if (m >= _m) { Rprintf("m >= _m"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[m*_n + n];
        }

        // NOTE: All assignment operators act on the caller!
        void operator+=(const MatI &A);
        void operator-=(const MatI &A);
        void operator*=(const MatI &A);
        void operator/=(const MatI &A);
        void operator+=(const int val) { _dat += val; }
        void operator-=(const int val) { _dat -= val; }
        void operator*=(const int val) { _dat *= val; }
        void operator/=(const int val) { _dat /= val; }


        void add(const MatI &toadd, MatI &out);
        void sub(const MatI &tosub, MatI &out);
        void mul(const MatI &tomul, MatI &out);
        void div(const MatI &todiv, MatI &out);

        // returns the transpose in out
        void transpose(MatI &out);

        void std_normal() { _dat.std_normal(); }
        void logarithm(double base) { _dat.logarithm(base); }
        void expand(MatI &result, int match, int expand_x_lt, int expand_x_rt, int expand_y_up, int expand_y_dn, int expand_diag_lt_up, int expand_diag_rt_up, int expand_diag_lt_dn, int expand_diag_rt_dn );

        void min_max(int &_min, int &_max) { _dat.min_max(_min,_max); }
        double avg() { return _dat.avg(); }
        //void operator++();
        //void operator--();

        int sum() { return _dat.sum(); } // return the sum of the entire matrix
        int sum(int m);  // return the sum of a given row
        // Returns in a vector all the values matching mask value
        void mask_as_vec(int return_val, MatI &mask, VecI &out);

        // prints the bare matrix as ascii
        void print(bool without_axes=0);
        void print(const char *file, bool without_axes=0);
        void print(std::ostream &fout, bool without_axes=0);

        // writes the matrix as binary (includes # rows first and # cols as
        // ints)
        void write(const char *file=NULL);


        // @TODO need to write these guys:
        // prints the matrix in binary format:
        // (int) num cols (int) num rows (int) data
//        void write(const char *file);
//        void write(std::ofstream &fout);

    private:
        void _copy(int *p1, const int *p2, int len) const;

}; // End class MatI



class MatD {

    public:
        // length
        int _m;
        int _n;
        VecD _dat;
        // Constructors:
        MatD();
        MatD(int m, int n);
        MatD(int m, int n, const double &val);


        // (copied from vec.h)
        // if (shallow == 1 (true)) then no memory is deleted upon destruction
        // if (shallow == 0 (false)) then delete[] is called
        // FOR THIS CONSTRUCTOR ONLY, there is no DEEP copying, EVER!
        MatD(int m, int n, double *arr, bool shallow=0);

        // (copied from vec.h)
        // if (shallow == 0 (false)) a DEEP copy is made of the data
        // if (shallow == 1 (true)) a copy of the pointer is made
        // if (shallow) then no memory is released upon destruction
        // shallow is used for a quick copy with which to work
        MatD(const MatD &A, bool shallow=0);

        operator double*() { return (double*)_dat; }
        operator const double*() { return (double*)_dat; }
        double* pointer() { return (double*)_dat; }
        double* pointer(int m) { return &_dat[m*_n]; }
        // creates vec objects
        // caller must have allocated the array for the vec objects
        // the data is a shallow copy!
        // transpose and call row_vecs for col_vecs!
        void row_vecs(int &cnt, VecD *vecs);

        MatD & operator=(const double &val);
        // DEEP
        MatD & operator=(MatD &A);
        ~MatD();
        // Deep copy unless shallow == true
        void copy(MatD &receiver, bool shallow=0) const;

        void set_from_ascii(std::ifstream &stream, int m, int n, MatD &out);
        void set_from_ascii(std::ifstream &stream, MatD &out);
        void set_from_ascii(const char *file, bool without_axes=0);
        void set_from_binary(const char *file);
        void file_rows_cols(std::ifstream &stream, int &rows, int &cols);
        // tnt_array2d_utils.h has a good example (use ifstream)

        // shallow copy and no ownership of memory
        void set(int m, int n, double *arr);
        // shallow copy and no ownership of memory
        void set(MatD &A);

        bool all_equal() {
            return _dat.all_equal();
        }

        // Deletes the object's memory (if not shallow) and takes ownership
        // of the array memory (we will call delete[])
        void take(int m, int n, double *arr);
        // Deletes previous memory (if not shallow) and takes ownership
        // of the other's memory.
        void take(MatD &A);

        // flattens the matrix and returns a vector
        void to_vec(VecD &outvec, bool shallow=0);

        bool operator==(const MatD &A);

        bool shallow() { return _dat.shallow(); }
        int dim1() const { return _m; }
        int dim2() const { return _n; }
        int mlen() const { return _m; }
        int nlen() const { return _n; }
        int rows() const { return _m; }
        int cols() const { return _n; }

        double& operator()(int m, int n) {
#ifdef JTP_BOUNDS_CHECK
            if (n < 0) { Rprintf("n < 0"); R_ShowMessage("Serious error in obiwarp.");
}
            if (n >= _n) { Rprintf("n >= _n"); R_ShowMessage("Serious error in obiwarp.");
}
            if (m < 0) { Rprintf("m < 0"); R_ShowMessage("Serious error in obiwarp.");
}
            if (m >= _m) { Rprintf("m >= _m"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[m*_n + n];
        }
        const double& operator()(int m, int n) const {
#ifdef JTP_BOUNDS_CHECK
            if (n < 0) { Rprintf("n < 0"); R_ShowMessage("Serious error in obiwarp.");
}
            if (n >= _n) { Rprintf("n >= _n"); R_ShowMessage("Serious error in obiwarp.");
}
            if (m < 0) { Rprintf("m < 0"); R_ShowMessage("Serious error in obiwarp.");
}
            if (m >= _m) { Rprintf("m >= _m"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[m*_n + n];
        }

        // NOTE: All assignment operators act on the caller!
        void operator+=(const MatD &A);
        void operator-=(const MatD &A);
        void operator*=(const MatD &A);
        void operator/=(const MatD &A);
        void operator+=(const double val) { _dat += val; }
        void operator-=(const double val) { _dat -= val; }
        void operator*=(const double val) { _dat *= val; }
        void operator/=(const double val) { _dat /= val; }


        void add(const MatD &toadd, MatD &out);
        void sub(const MatD &tosub, MatD &out);
        void mul(const MatD &tomul, MatD &out);
        void div(const MatD &todiv, MatD &out);

        // returns the transpose in out
        void transpose(MatD &out);

        void std_normal() { _dat.std_normal(); }
        void logarithm(double base) { _dat.logarithm(base); }
        void expand(MatD &result, double match, int expand_x_lt, int expand_x_rt, int expand_y_up, int expand_y_dn, int expand_diag_lt_up, int expand_diag_rt_up, int expand_diag_lt_dn, int expand_diag_rt_dn );

        void min_max(double &_min, double &_max) { _dat.min_max(_min,_max); }
        double avg() { return _dat.avg(); }
        //void operator++();
        //void operator--();

        double sum() { return _dat.sum(); } // return the sum of the entire matrix
        double sum(int m);  // return the sum of a given row
        // Returns in a vector all the values matching mask value
        void mask_as_vec(double return_val, MatI &mask, VecD &out);

        // prints the bare matrix as ascii
        void print(bool without_axes=0);
        void print(const char *file, bool without_axes=0);
        void print(std::ostream &fout, bool without_axes=0);

        // writes the matrix as binary (includes # rows first and # cols as
        // ints)
        void write(const char *file=NULL);


        // @TODO need to write these guys:
        // prints the matrix in binary format:
        // (int) num cols (int) num rows (double) data
//        void write(const char *file);
//        void write(std::ofstream &fout);

    private:
        void _copy(double *p1, const double *p2, int len) const;

}; // End class MatD



class MatF {

    public:
        // length
        int _m;
        int _n;
        VecF _dat;
        // Constructors:
        MatF();
        MatF(int m, int n);
        MatF(int m, int n, const float &val);


        // (copied from vec.h)
        // if (shallow == 1 (true)) then no memory is deleted upon destruction
        // if (shallow == 0 (false)) then delete[] is called
        // FOR THIS CONSTRUCTOR ONLY, there is no DEEP copying, EVER!
        MatF(int m, int n, float *arr, bool shallow=0);

        // (copied from vec.h)
        // if (shallow == 0 (false)) a DEEP copy is made of the data
        // if (shallow == 1 (true)) a copy of the pointer is made
        // if (shallow) then no memory is released upon destruction
        // shallow is used for a quick copy with which to work
        MatF(const MatF &A, bool shallow=0);

        operator float*() { return (float*)_dat; }
        operator const float*() { return (float*)_dat; }
        float* pointer() { return (float*)_dat; }
        float* pointer(int m) { return &_dat[m*_n]; }
        // creates vec objects
        // caller must have allocated the array for the vec objects
        // the data is a shallow copy!
        // transpose and call row_vecs for col_vecs!
        void row_vecs(int &cnt, VecF *vecs);

        MatF & operator=(const float &val);
        // DEEP
        MatF & operator=(MatF &A);
        ~MatF();
        // Deep copy unless shallow == true
        void copy(MatF &receiver, bool shallow=0) const;

        void set_from_ascii(std::ifstream &stream, int m, int n, MatF &out);
        void set_from_ascii(std::ifstream &stream, MatF &out);
        void set_from_ascii(const char *file, bool without_axes=0);
        void set_from_binary(const char *file);
        void file_rows_cols(std::ifstream &stream, int &rows, int &cols);
        // tnt_array2d_utils.h has a good example (use ifstream)

        // shallow copy and no ownership of memory
        void set(int m, int n, float *arr);
        // shallow copy and no ownership of memory
        void set(MatF &A);

        bool all_equal() {
            return _dat.all_equal();
        }

        // Deletes the object's memory (if not shallow) and takes ownership
        // of the array memory (we will call delete[])
        void take(int m, int n, float *arr);
        // Deletes previous memory (if not shallow) and takes ownership
        // of the other's memory.
        void take(MatF &A);

        // flattens the matrix and returns a vector
        void to_vec(VecF &outvec, bool shallow=0);

        bool operator==(const MatF &A);

        bool shallow() { return _dat.shallow(); }
        int dim1() const { return _m; }
        int dim2() const { return _n; }
        int mlen() const { return _m; }
        int nlen() const { return _n; }
        int rows() const { return _m; }
        int cols() const { return _n; }

        float& operator()(int m, int n) {
#ifdef JTP_BOUNDS_CHECK
            if (n < 0) { Rprintf("n < 0"); R_ShowMessage("Serious error in obiwarp.");
}
            if (n >= _n) { Rprintf("n >= _n"); R_ShowMessage("Serious error in obiwarp.");
}
            if (m < 0) { Rprintf("m < 0"); R_ShowMessage("Serious error in obiwarp.");
}
            if (m >= _m) { Rprintf("m >= _m"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[m*_n + n];
        }
        const float& operator()(int m, int n) const {
#ifdef JTP_BOUNDS_CHECK
            if (n < 0) { Rprintf("n < 0"); R_ShowMessage("Serious error in obiwarp.");
}
            if (n >= _n) { Rprintf("n >= _n"); R_ShowMessage("Serious error in obiwarp.");
}
            if (m < 0) { Rprintf("m < 0"); R_ShowMessage("Serious error in obiwarp.");
}
            if (m >= _m) { Rprintf("m >= _m"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[m*_n + n];
        }

        // NOTE: All assignment operators act on the caller!
        void operator+=(const MatF &A);
        void operator-=(const MatF &A);
        void operator*=(const MatF &A);
        void operator/=(const MatF &A);
        void operator+=(const float val) { _dat += val; }
        void operator-=(const float val) { _dat -= val; }
        void operator*=(const float val) { _dat *= val; }
        void operator/=(const float val) { _dat /= val; }


        void add(const MatF &toadd, MatF &out);
        void sub(const MatF &tosub, MatF &out);
        void mul(const MatF &tomul, MatF &out);
        void div(const MatF &todiv, MatF &out);

        // returns the transpose in out
        void transpose(MatF &out);

        void std_normal() { _dat.std_normal(); }
        void logarithm(double base) { _dat.logarithm(base); }
        void expand(MatF &result, float match, int expand_x_lt, int expand_x_rt, int expand_y_up, int expand_y_dn, int expand_diag_lt_up, int expand_diag_rt_up, int expand_diag_lt_dn, int expand_diag_rt_dn );

        void min_max(float &_min, float &_max) { _dat.min_max(_min,_max); }
        double avg() { return _dat.avg(); }
        //void operator++();
        //void operator--();

        float sum() { return _dat.sum(); } // return the sum of the entire matrix
        float sum(int m);  // return the sum of a given row
        // Returns in a vector all the values matching mask value
        void mask_as_vec(float return_val, MatI &mask, VecF &out);

        // prints the bare matrix as ascii
        void print(bool without_axes=0);

        void print(int __m, int __n,bool without_axes);

        void print(const char *file, bool without_axes=0);
        void print(std::ostream &fout, bool without_axes=0);

        // writes the matrix as binary (includes # rows first and # cols as
        // ints)
        void write(const char *file=NULL);


        // @TODO need to write these guys:
        // prints the matrix in binary format:
        // (int) num cols (int) num rows (float) data
//        void write(const char *file);
//        void write(std::ofstream &fout);

    private:
        void _copy(float *p1, const float *p2, int len) const;

}; // End class MatF

// END TEMPLATE

} // End namespace

#endif
