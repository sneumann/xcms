#ifndef _VEC_H
#define _VEC_H

#include <fstream>

/*************************************************************
 * Creation from existing object/array is always deep!.
 * If you want shallow, use a pointer!
 ************************************************************/


namespace VEC {

class VecF;
class VecD;
class VecI;

// BEGIN TEMPLATE

class VecI {

    protected:
        // length
        int _n;
        int *_dat;
        bool _shallow;

    public:
        // Constructors:
        VecI();
        // Data values are NOT set by default
        explicit VecI(int n);
        VecI(int n, const int &val);

        // if (shallow == 1 (true)) then no memory is deleted upon destruction
        // if (shallow == 0 (false)) then delete[] is called
        // FOR THIS CONSTRUCTOR ONLY, there is no DEEP copying, EVER!
        VecI(int n, int *arr, bool shallow=0);

        // if (shallow == 0 (false)) a DEEP copy is made of the data
        // if (shallow == 1 (true)) a copy of the pointer is made
        // if (shallow) then no memory is released upon destruction
        // shallow is used for a quick copy with which to work
        VecI(const VecI &A, bool shallow=0);

        operator int*() {
            if (_n > 0) {
                return &(_dat[0]);
            }
            else {
                return 0;
            }
        }
        operator const int*() {
            if (_n > 0) {
                return &(_dat[0]);
            }
            else {
                return 0;
            }
        }

        int first() { return _dat[0]; }
        int last() { return _dat[_n-1]; }

        int* pointer() { return &(_dat[0]); }

        // Returns the name of the class
        // del needs to be called
        char * class_name();
        // shallow ownership
        void set(int n, int *arr);
        // Deletes the object's previous memory (if not shallow) and takes
        // ownership of the array (destructor call delete[])
        // shallow in this context only refers to calling delete
        // no data is copied
        // deep ownership (no copy is performed)
        void take(int n, int *arr);

        void to_f(VecF &out);
        void to_i(VecI &out);

        // returns the first index at the value, else -1
        int index(int val);

        // shallow ownership
        void set(VecI &A);
        // Deletes previous memory (if not shallow) and takes ownership
        // of the other's memory.
        void take(VecI &A);
        VecI & operator=(const int &val);
        VecI & operator=(VecI &A);
        ~VecI();
        // A deep copy unless shallow is set
        void copy(VecI &receiver, bool shallow=0) const;

        bool operator==(const VecI &A);

        int length() const { return _n; }
        int len() const { return _n; }
        int size() const { return _n; }
        int dim() const { return _n; }
        int dim1() const { return _n; }
        // Returns in a vector all the values matching mask value
        void mask_as_vec(int return_val, VecI &mask, VecI &vec);

        // Returns true if all values are the same, false otherwise
        bool all_equal() {
            int _min, _max; min_max(_min, _max);
            if (_min == _max) { return 1; }
            else { return 0; }
        }

        int& operator[](int i) {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { Rprintf("index < 0 !"); R_ShowMessage("Serious error in obiwarp.");
}
            if (i >= _n) { Rprintf("i >= _n !"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[i];
        }
        const int& operator[](int i) const {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { Rprintf("index < 0 !"); R_ShowMessage("Serious error in obiwarp.");
}
            if (i >= _n) { Rprintf("i >= _n !"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[i];
        }

        int& at(int i) {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { Rprintf("index < 0 !"); R_ShowMessage("Serious error in obiwarp.");
}
            if (i >= _n) { Rprintf("i >= _n !"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[i];
        }
        const int& at(int i) const {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { Rprintf("index < 0 !"); R_ShowMessage("Serious error in obiwarp.");
}
            if (i >= _n) { Rprintf("i >= _n !"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[i];
        }

        bool shallow() {
            return _shallow;
        }

        // NOTE: All operators act on the caller!
        // Operators
        void operator+=(const VecI &A);
        void operator-=(const VecI &A);
        void operator*=(const VecI &A);
        void operator/=(const VecI &A);
        void operator+=(int val);
        void operator-=(int val);
        void operator*=(int val);
        void operator/=(int val);

        void add(const VecI &toadd, VecI &out);
        void sub(const VecI &tosub, VecI &out);
        void mul(const VecI &tomul, VecI &out);
        void div(const VecI &todiv, VecI &out);
        // This may be slow because we cast every value to double regardless
        void square_root();

        void logarithm(double base);
        void min_max(int &mn, int &mx);
        // alias for min_max
        void mn_mx(int &mn, int &mx) {min_max(mn,mx);}
        double avg() const;
        void hist(int num_bins, VecD &bins, VecI &freqs);
        void sample_stats(double &mean, double &std_dev);
        double prob_one_side_right(double x);
        int sum();
        int sum_of_sq();

        void abs_val();
        // converts the distribution of values into standard normal
        // Only for floating points right now!
        void std_normal();

        // uses quicksort to sort the values
        void sort();
        static int intCompare( const void *a, const void *b );

        // Removes the value at index and shortens the array by one
        // not shallow anymore regardless of previous state
        void remove(int index);

        // prints the vector (space delimited on one line without any length
        // header)
        void print(bool without_length=0);
        // prints the vector to the file with the length written on the line
        // before, unless without_length == true
        void print(const char *, bool without_length=0);
        // prints the vector to the filehandle with the length written on the
        // line before, unless without_length == true
        void print(std::ostream &fout, bool without_length=0);

        // CLASS FUNCTIONS:
        static int pchst(int arg1, int arg2) {
            if      (arg1*arg2 > 0) { return  1; }
            else if (arg1*arg2 < 0) { return -1; }
            else                    { return  0; }
        }

        static double pearsons_r(VecI &x, VecI &y);
        static double covariance(VecI &x, VecI &y);
        static double euclidean(VecI &x, VecI &y);
        static int dot_product(VecI &x, VecI &y);

        static void xy_to_x(VecI &x, VecI &y);
        static void x_to_xy(VecI &x, VecI &y);
        static void chim(VecI &x, VecI &y, VecI &out_derivs);

        static void calc_cubic_coeff(VecI &x, VecI &y, VecI &derivs, VecI &c2, VecI &c3);
        static void chfe(VecI &xin, VecI &yin, VecI &xe, VecI &out_ye, int sorted=0);
        //static void pchfe(VecI &xin, VecI &yin, VecI &XE, VecI &out_newy);
        // interpolates so that linearity is encouraged along x axis
        // if out_new_y.length() == 0 then new memory is allocated
        // otherwise, uses whatever memory is allocated in out_new_y
        static inline void chfev(int X1, int F1, int D1, int C2, int C3, int XE, int &FE);
        static inline void chfev_all(int X1, int X2, int F1, int F2, int D1, int D2, int XE, int &FE);
       // static void chfev(int X1, int X2, int F1, int F2, int D1, int D2, int NE, int *XE, int *FE, int *nlr, int &ierr);

        // interpolates so that linearity is encouraged along the xy line
        // if out_new_y.length() == 0 then new memory is allocated
        // otherwise, uses whatever memory is allocated in out_new_y
        static void chfe_xy(VecI &x, VecI &y, VecI &new_x, VecI &out_new_y, int sorted=0);

        static void linear_derivs(VecI &x, VecI &y, VecI &out_derivs);
        static void linear_interp(VecI &xin, VecI &yin, VecI &xe, VecI &out_ye, int sorted=0);
        //##### FOR ANY FUNCTION:
        //                B
        //               /|
        //              / |
        //             /  |
        //          c /   |a
        //           /    |
        //          /     |
        //       A --------C
        //             b
        //
        //  cPerp = Perpendicular from C to c
        //  (sin A)/a = (sin C)/c
        //  sin C = 1
        //  c = sqrt(a^2 + b^2)
        //  sin A = a/c
        //  cPerp = (a/c)*b      note:   = sin A * b
        //  cPerp^2 = (ab)^2/(a^2 + b^2)

        // g(x):  y = x
        // h(y):  x = y
        //  a = y - g(x)

        // b = actual x - (x = actual y)
        // a = actual y - (y = actual x)
        // avg residual = SUM(cPerp^2)/N
        //

        //###### FOR X=Y
        // residual^2 = 1/2((Y-X)^2)        note: [or X-Y]
        static double sum_sq_res_yeqx(VecI &x, VecI &y);
        // divides the sum of square of residuals by the length of the vector
        static double avg_sq_res_yeqx(VecI &x, VecI &y);

        // returns the average of the absolute values of the differences at
        // each index
        static double avg_abs_diff(VecI &x, VecI &y);
        static void rsq_slope_intercept(VecI &x, VecI &y, double &rsq, double &slope, double &y_intercept);

    private:
        void _copy(int *p1, const int *p2, int len) const;
        static void outliers_from_regression_line(VecI &x, VecI &y, VecI &indices_out);
        double _zScore(double mean, double sigma, double x);

}; // End class VecI



class VecD {

    protected:
        // length
        int _n;
        double *_dat;
        bool _shallow;

    public:
        // Constructors:
        VecD();
        // Data values are NOT set by default
        explicit VecD(int n);
        VecD(int n, const double &val);

        // if (shallow == 1 (true)) then no memory is deleted upon destruction
        // if (shallow == 0 (false)) then delete[] is called
        // FOR THIS CONSTRUCTOR ONLY, there is no DEEP copying, EVER!
        VecD(int n, double *arr, bool shallow=0);

        // if (shallow == 0 (false)) a DEEP copy is made of the data
        // if (shallow == 1 (true)) a copy of the pointer is made
        // if (shallow) then no memory is released upon destruction
        // shallow is used for a quick copy with which to work
        VecD(const VecD &A, bool shallow=0);

        operator double*() {
            if (_n > 0) {
                return &(_dat[0]);
            }
            else {
                return 0;
            }
        }
        operator const double*() {
            if (_n > 0) {
                return &(_dat[0]);
            }
            else {
                return 0;
            }
        }

        double first() { return _dat[0]; }
        double last() { return _dat[_n-1]; }

        double* pointer() { return &(_dat[0]); }

        // Returns the name of the class
        // del needs to be called
        char * class_name();
        // shallow ownership
        void set(int n, double *arr);
        // Deletes the object's previous memory (if not shallow) and takes
        // ownership of the array (destructor call delete[])
        // shallow in this context only refers to calling delete
        // no data is copied
        // deep ownership (no copy is performed)
        void take(int n, double *arr);

        void to_f(VecF &out);
        void to_i(VecI &out);

        // returns the first index at the value, else -1
        int index(double val);

        // shallow ownership
        void set(VecD &A);
        // Deletes previous memory (if not shallow) and takes ownership
        // of the other's memory.
        void take(VecD &A);
        VecD & operator=(const double &val);
        VecD & operator=(VecD &A);
        ~VecD();
        // A deep copy unless shallow is set
        void copy(VecD &receiver, bool shallow=0) const;

        bool operator==(const VecD &A);

        int length() const { return _n; }
        int len() const { return _n; }
        int size() const { return _n; }
        int dim() const { return _n; }
        int dim1() const { return _n; }
        // Returns in a vector all the values matching mask value
        void mask_as_vec(double return_val, VecI &mask, VecD &vec);

        // Returns true if all values are the same, false otherwise
        bool all_equal() {
            double _min, _max; min_max(_min, _max);
            if (_min == _max) { return 1; }
            else { return 0; }
        }

        double& operator[](int i) {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { Rprintf("index < 0 !"); R_ShowMessage("Serious error in obiwarp.");
}
            if (i >= _n) { Rprintf("i >= _n !"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[i];
        }
        const double& operator[](int i) const {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { Rprintf("index < 0 !"); R_ShowMessage("Serious error in obiwarp.");
}
            if (i >= _n) { Rprintf("i >= _n !"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[i];
        }

        double& at(int i) {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { Rprintf("index < 0 !"); R_ShowMessage("Serious error in obiwarp.");
}
            if (i >= _n) { Rprintf("i >= _n !"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[i];
        }
        const double& at(int i) const {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { Rprintf("index < 0 !"); R_ShowMessage("Serious error in obiwarp.");
}
            if (i >= _n) { Rprintf("i >= _n !"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[i];
        }

        bool shallow() {
            return _shallow;
        }

        // NOTE: All operators act on the caller!
        // Operators
        void operator+=(const VecD &A);
        void operator-=(const VecD &A);
        void operator*=(const VecD &A);
        void operator/=(const VecD &A);
        void operator+=(double val);
        void operator-=(double val);
        void operator*=(double val);
        void operator/=(double val);

        void add(const VecD &toadd, VecD &out);
        void sub(const VecD &tosub, VecD &out);
        void mul(const VecD &tomul, VecD &out);
        void div(const VecD &todiv, VecD &out);
        // This may be slow because we cast every value to double regardless
        void square_root();

        void logarithm(double base);
        void min_max(double &mn, double &mx);
        // alias for min_max
        void mn_mx(double &mn, double &mx) {min_max(mn,mx);}
        double avg() const;
        void hist(int num_bins, VecD &bins, VecI &freqs);
        void sample_stats(double &mean, double &std_dev);
        double prob_one_side_right(double x);
        double sum();
        double sum_of_sq();

        void abs_val();
        // converts the distribution of values into standard normal
        // Only for floating points right now!
        void std_normal();

        // uses quicksort to sort the values
        void sort();
        static int doubleCompare( const void *a, const void *b );

        // Removes the value at index and shortens the array by one
        // not shallow anymore regardless of previous state
        void remove(int index);


        // prints the vector (space delimited on one line without any length
        // header)
        void print(bool without_length=0);
        // prints the vector to the file with the length written on the line
        // before, unless without_length == true
        void print(const char *, bool without_length=0);
        // prints the vector to the filehandle with the length written on the
        // line before, unless without_length == true
        void print(std::ostream &fout, bool without_length=0);

        // CLASS FUNCTIONS:
        static int pchst(double arg1, double arg2) {
            if      (arg1*arg2 > 0) { return  1; }
            else if (arg1*arg2 < 0) { return -1; }
            else                    { return  0; }
        }

        static double pearsons_r(VecD &x, VecD &y);
        static double covariance(VecD &x, VecD &y);
        static double euclidean(VecD &x, VecD &y);
        static double dot_product(VecD &x, VecD &y);

        static void xy_to_x(VecD &x, VecD &y);
        static void x_to_xy(VecD &x, VecD &y);
        static void chim(VecD &x, VecD &y, VecD &out_derivs);

        static void calc_cubic_coeff(VecD &x, VecD &y, VecD &derivs, VecD &c2, VecD &c3);
        static void chfe(VecD &xin, VecD &yin, VecD &xe, VecD &out_ye, int sorted=0);
        //static void pchfe(VecD &xin, VecD &yin, VecD &XE, VecD &out_newy);
        // interpolates so that linearity is encouraged along x axis
        // if out_new_y.length() == 0 then new memory is allocated
        // otherwise, uses whatever memory is allocated in out_new_y
        static inline void chfev(double X1, double F1, double D1, double C2, double C3, double XE, double &FE);
        static inline void chfev_all(double X1, double X2, double F1, double F2, double D1, double D2, double XE, double &FE);
       // static void chfev(double X1, double X2, double F1, double F2, double D1, double D2, int NE, double *XE, double *FE, int *nlr, int &ierr);

        // interpolates so that linearity is encouraged along the xy line
        // if out_new_y.length() == 0 then new memory is allocated
        // otherwise, uses whatever memory is allocated in out_new_y
        static void chfe_xy(VecD &x, VecD &y, VecD &new_x, VecD &out_new_y, int sorted=0);

        static void linear_derivs(VecD &x, VecD &y, VecD &out_derivs);
        static void linear_interp(VecD &xin, VecD &yin, VecD &xe, VecD &out_ye, int sorted=0);
        //##### FOR ANY FUNCTION:
        //                B
        //               /|
        //              / |
        //             /  |
        //          c /   |a
        //           /    |
        //          /     |
        //       A --------C
        //             b
        //
        //  cPerp = Perpendicular from C to c
        //  (sin A)/a = (sin C)/c
        //  sin C = 1
        //  c = sqrt(a^2 + b^2)
        //  sin A = a/c
        //  cPerp = (a/c)*b      note:   = sin A * b
        //  cPerp^2 = (ab)^2/(a^2 + b^2)

        // g(x):  y = x
        // h(y):  x = y
        //  a = y - g(x)

        // b = actual x - (x = actual y)
        // a = actual y - (y = actual x)
        // avg residual = SUM(cPerp^2)/N
        //

        //###### FOR X=Y
        // residual^2 = 1/2((Y-X)^2)        note: [or X-Y]
        static double sum_sq_res_yeqx(VecD &x, VecD &y);
        // divides the sum of square of residuals by the length of the vector
        static double avg_sq_res_yeqx(VecD &x, VecD &y);

        // returns the average of the absolute values of the differences at
        // each index
        static double avg_abs_diff(VecD &x, VecD &y);
        static void rsq_slope_intercept(VecD &x, VecD &y, double &rsq, double &slope, double &y_intercept);

    private:
        void _copy(double *p1, const double *p2, int len) const;
        static void outliers_from_regression_line(VecD &x, VecD &y, VecI &indices_out);
        double _zScore(double mean, double sigma, double x);

}; // End class VecD



class VecF {

    protected:
        // length
        int _n;
        float *_dat;
        bool _shallow;

    public:
        // Constructors:
        VecF();
        // Data values are NOT set by default
        explicit VecF(int n);
        VecF(int n, const float &val);

        // if (shallow == 1 (true)) then no memory is deleted upon destruction
        // if (shallow == 0 (false)) then delete[] is called
        // FOR THIS CONSTRUCTOR ONLY, there is no DEEP copying, EVER!
        VecF(int n, float *arr, bool shallow=0);

        // if (shallow == 0 (false)) a DEEP copy is made of the data
        // if (shallow == 1 (true)) a copy of the pointer is made
        // if (shallow) then no memory is released upon destruction
        // shallow is used for a quick copy with which to work
        VecF(const VecF &A, bool shallow=0);

        operator float*() {
            if (_n > 0) {
                return &(_dat[0]);
            }
            else {
                return 0;
            }
        }
        operator const float*() {
            if (_n > 0) {
                return &(_dat[0]);
            }
            else {
                return 0;
            }
        }

        float first() { return _dat[0]; }
        float last() { return _dat[_n-1]; }

        float* pointer() { return &(_dat[0]); }

        // Returns the name of the class
        // del needs to be called
        char * class_name();
        // shallow ownership
        void set(int n, float *arr);
        // Deletes the object's previous memory (if not shallow) and takes
        // ownership of the array (destructor call delete[])
        // shallow in this context only refers to calling delete
        // no data is copied
        // deep ownership (no copy is performed)
        void take(int n, float *arr);

        void to_f(VecF &out);
        void to_i(VecI &out);

        // returns the first index at the value, else -1
        int index(float val);

        // shallow ownership
        void set(VecF &A);
        // Deletes previous memory (if not shallow) and takes ownership
        // of the other's memory.
        void take(VecF &A);
        VecF & operator=(const float &val);
        VecF & operator=(VecF &A);
        ~VecF();
        // A deep copy unless shallow is set
        void copy(VecF &receiver, bool shallow=0) const;

        bool operator==(const VecF &A);

        int length() const { return _n; }
        int len() const { return _n; }
        int size() const { return _n; }
        int dim() const { return _n; }
        int dim1() const { return _n; }
        // Returns in a vector all the values matching mask value
        void mask_as_vec(float return_val, VecI &mask, VecF &vec);

        // Returns true if all values are the same, false otherwise
        bool all_equal() {
            float _min, _max; min_max(_min, _max);
            if (_min == _max) { return 1; }
            else { return 0; }
        }

        float& operator[](int i) {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { Rprintf("index < 0 !"); R_ShowMessage("Serious error in obiwarp.");
}
            if (i >= _n) { Rprintf("i >= _n !"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[i];
        }
        const float& operator[](int i) const {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { Rprintf("index < 0 !"); R_ShowMessage("Serious error in obiwarp.");
}
            if (i >= _n) { Rprintf("i >= _n !"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[i];
        }

        float& at(int i) {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { Rprintf("index < 0 !"); R_ShowMessage("Serious error in obiwarp.");
}
            if (i >= _n) { Rprintf("i >= _n !"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[i];
        }
        const float& at(int i) const {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { Rprintf("index < 0 !"); R_ShowMessage("Serious error in obiwarp.");
}
            if (i >= _n) { Rprintf("i >= _n !"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
            return _dat[i];
        }

        bool shallow() {
            return _shallow;
        }

        // NOTE: All operators act on the caller!
        // Operators
        void operator+=(const VecF &A);
        void operator-=(const VecF &A);
        void operator*=(const VecF &A);
        void operator/=(const VecF &A);
        void operator+=(float val);
        void operator-=(float val);
        void operator*=(float val);
        void operator/=(float val);

        void add(const VecF &toadd, VecF &out);
        void sub(const VecF &tosub, VecF &out);
        void mul(const VecF &tomul, VecF &out);
        void div(const VecF &todiv, VecF &out);
        // This may be slow because we cast every value to double regardless
        void square_root();

        void logarithm(double base);
        void min_max(float &mn, float &mx);
        // alias for min_max
        void mn_mx(float &mn, float &mx) {min_max(mn,mx);}
        double avg() const;
        void hist(int num_bins, VecD &bins, VecI &freqs);
        void sample_stats(double &mean, double &std_dev);
        double prob_one_side_right(double x);
        float sum();
        float sum_of_sq();

        void abs_val();
        // converts the distribution of values into standard normal
        // Only for floating points right now!
        void std_normal();

        // uses quicksort to sort the values
        void sort();
        static int floatCompare( const void *a, const void *b );

        // Removes the value at index and shortens the array by one
        // not shallow anymore regardless of previous state
        void remove(int index);

        // prints the vector (space delimited on one line without any length
        // header)
        void print(bool without_length=0);
        // prints the vector to the file with the length written on the line
        // before, unless without_length == true
        void print(const char *, bool without_length=0);
        // prints the vector to the filehandle with the length written on the
        // line before, unless without_length == true
        void print(std::ostream &fout, bool without_length=0);

        void print_tm();
        float* back(){ return _dat; }

        // CLASS FUNCTIONS:
        static int pchst(float arg1, float arg2) {
            if      (arg1*arg2 > 0) { return  1; }
            else if (arg1*arg2 < 0) { return -1; }
            else                    { return  0; }
        }

        static double pearsons_r(VecF &x, VecF &y);
        static double covariance(VecF &x, VecF &y);
        static double euclidean(VecF &x, VecF &y);
        static float dot_product(VecF &x, VecF &y);

        static void xy_to_x(VecF &x, VecF &y);
        static void x_to_xy(VecF &x, VecF &y);
        static void chim(VecF &x, VecF &y, VecF &out_derivs);

        static void calc_cubic_coeff(VecF &x, VecF &y, VecF &derivs, VecF &c2, VecF &c3);
        static void chfe(VecF &xin, VecF &yin, VecF &xe, VecF &out_ye, int sorted=0);
        //static void pchfe(VecF &xin, VecF &yin, VecF &XE, VecF &out_newy);
        // interpolates so that linearity is encouraged along x axis
        // if out_new_y.length() == 0 then new memory is allocated
        // otherwise, uses whatever memory is allocated in out_new_y
        static inline void chfev(float X1, float F1, float D1, float C2, float C3, float XE, float &FE);
        static inline void chfev_all(float X1, float X2, float F1, float F2, float D1, float D2, float XE, float &FE);
       // static void chfev(float X1, float X2, float F1, float F2, float D1, float D2, int NE, float *XE, float *FE, int *nlr, int &ierr);

        // interpolates so that linearity is encouraged along the xy line
        // if out_new_y.length() == 0 then new memory is allocated
        // otherwise, uses whatever memory is allocated in out_new_y
        static void chfe_xy(VecF &x, VecF &y, VecF &new_x, VecF &out_new_y, int sorted=0);

        static void linear_derivs(VecF &x, VecF &y, VecF &out_derivs);
        static void linear_interp(VecF &xin, VecF &yin, VecF &xe, VecF &out_ye, int sorted=0);
        //##### FOR ANY FUNCTION:
        //                B
        //               /|
        //              / |
        //             /  |
        //          c /   |a
        //           /    |
        //          /     |
        //       A --------C
        //             b
        //
        //  cPerp = Perpendicular from C to c
        //  (sin A)/a = (sin C)/c
        //  sin C = 1
        //  c = sqrt(a^2 + b^2)
        //  sin A = a/c
        //  cPerp = (a/c)*b      note:   = sin A * b
        //  cPerp^2 = (ab)^2/(a^2 + b^2)

        // g(x):  y = x
        // h(y):  x = y
        //  a = y - g(x)

        // b = actual x - (x = actual y)
        // a = actual y - (y = actual x)
        // avg residual = SUM(cPerp^2)/N
        //

        //###### FOR X=Y
        // residual^2 = 1/2((Y-X)^2)        note: [or X-Y]
        static double sum_sq_res_yeqx(VecF &x, VecF &y);
        // divides the sum of square of residuals by the length of the vector
        static double avg_sq_res_yeqx(VecF &x, VecF &y);

        // returns the average of the absolute values of the differences at
        // each index
        static double avg_abs_diff(VecF &x, VecF &y);
        static void rsq_slope_intercept(VecF &x, VecF &y, double &rsq, double &slope, double &y_intercept);

    private:
        void _copy(float *p1, const float *p2, int len) const;
        static void outliers_from_regression_line(VecF &x, VecF &y, VecI &indices_out);
        double _zScore(double mean, double sigma, double x);

}; // End class VecF

// END TEMPLATE

} // End namespace

#endif
