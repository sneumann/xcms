
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include "vec.h"

#include <R.h> // for Rprintf

#ifndef min
	#define min(a,b) ( ( (a) < (b) ) ? (a) : (b) )
#endif
#ifndef max
	#define max(a,b) ( ( (a) > (b) ) ? (a) : (b) )
#endif


namespace VEC {

// BEGIN TEMPLATE


/****************************************************************
 * VecI
 ***************************************************************/

// Constructors:
VecI::VecI() : _n(0), _shallow(true) {
#ifdef JTP_DEBUG
    Rprintf("Creating DATA (NO ARGS)");
#endif
}

VecI::VecI(int n) : _n(n), _shallow(false) {
#ifdef JTP_BOUNDS_CHECK
    if (n < 0) { Rprintf("n < 0, exiting"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
    _dat = new int[_n];
#ifdef JTP_DEBUG
    Rprintf("Creating DATA(N)");
#endif
}

VecI::VecI(int n, const int &val) : _n(n), _shallow(false) {
    _dat = new int[_n];
    for (int i = 0; i < _n; ++i) {
        _dat[i] = val;
    }
#ifdef JTP_DEBUG
    Rprintf("Creating DATA(N,int)");
#endif
}

VecI::VecI(int n, int *arr, bool shallow) : _n(n), _dat(arr), _shallow(shallow) {
#ifdef JTP_DEBUG
    Rprintf("SHALLOW, (N,*ARR)");
#endif
}

VecI::VecI(const VecI &A, bool shallow) : _n(A._n), _shallow(shallow) {
    if (!shallow) {
        _dat = new int[_n];
        for (int i = 0; i < _n; ++i) {
            _dat[i] = A._dat[i];
        }
    }
    else {
        _dat = A._dat;
    }
#ifdef JTP_DEBUG
    Rprintf("created with VecI(const VecI &A)");
#endif
}

void VecI::to_f(VecF &out) {
    VecF _tmp(_n);
    for (int i = 0; i < _n; ++i) {
        _tmp[i] = (float)(_dat[i]);
    }
    out.take(_tmp);
}

void VecI::to_i(VecI &out) {
    VecI _tmp(_n);
    for (int i = 0; i < _n; ++i) {
        _tmp[i] = (int)(_dat[i]);
    }
    out.take(_tmp);
}


void VecI::set(int n, int *arr) {
    if (!_shallow) {
        delete[] _dat;
    }
    _dat = arr;
    _shallow = true;
    _n = n;
}

void VecI::take(int n, int *arr) {
    if (!_shallow) {
        delete[] _dat;
    }
    _dat = arr;
    _shallow = false;
    _n = n;
}

void VecI::set(VecI &A) {
    if (!_shallow) {
        delete[] _dat;
    }
    _dat = A._dat;
    _shallow = true;
    _n = A._n;
}


void VecI::take(VecI &A) {
    if (!_shallow) {
        delete[] _dat;
    }
    if (A._shallow) {
        Rprintf("Can't take ownership of memory of a shallow Vec!");
        R_ShowMessage("Serious error in obiwarp.");
    }
    _dat = A._dat;
    A._shallow = true;
    _shallow = false;
    _n = A._n;
}


bool VecI::operator==(const VecI &A) {
    // We don't care if one is shallow and the A is not!
    if (A._n == _n) { // Same size
        if (A._dat == _dat) { return true; }  // Same data
        else {
            for (int i = 0; i < _n; ++i) {
                if (A._dat[i] != _dat[i]) { return false; }
            }
            return true;
        }
    }
    else {
        return false;
    }
}

void VecI::copy(VecI &receiver, bool shallow) const {
    if (!receiver._shallow) {
        delete[] receiver._dat;
    }
    if (shallow) {
        receiver._dat = _dat;
        receiver._n = _n;
        receiver._shallow = true;
    }
    else {
        receiver._n = _n;
        receiver._dat = new int[_n];
        _copy(receiver._dat, _dat, _n);
        receiver._shallow = false;
    }
}

void VecI::_copy(int *p1, const int *p2, int len) const {
    // This is slightly faster on gcc Linux Mandrake
    for (int i = 0; i < len; ++i) {
        p1[i] = p2[i];
    }
}

VecI & VecI::operator=(const int &val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] = val;
    }
    return *this;
}

VecI & VecI::operator=(VecI &A) {
#ifdef JTP_DEBUG
    Rprintf("IN ASSIGNMENT OP tOP");
#endif
    if (this != &A) {
#ifdef JTP_DEBUG
        Rprintf("IN ASSIGNMENT OP MID");
#endif
        if (!_shallow) {
            delete[] _dat;
        }
        _n = A._n;
        _dat = new int[_n];
        _copy(_dat, A._dat, _n);
        _shallow = false;
    }
    return *this;
}

VecI::~VecI( ) {
#ifdef JTP_DEBUG
    Rprintf("DESTRUCTOR");
#endif
    if (!_shallow) {
#ifdef JTP_DEBUG
        Rprintf("DELETING DATA");
#endif
        delete[] _dat;
    }
}

/*************************
 * MATH OPERATORS
 ************************/

void VecI::operator+=(int val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] += val;
    }
}
void VecI::operator-=(int val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] -= val;
    }
}
void VecI::operator*=(int val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] *= val;
    }
}
void VecI::operator/=(int val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] /= val;
    }
}
void VecI::operator+=(const VecI &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] += A[i];
    }
}

void VecI::operator-=(const VecI &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] -= A[i];
    }
}

void VecI::operator*=(const VecI &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] *= A[i];
    }
}

void VecI::operator/=(const VecI &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] /= A[i];
    }
}

void VecI::add(const VecI &toadd, VecI &out) {
    if (toadd._n == _n) {
        int *_tmparr = new int[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] + toadd[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}

void VecI::sub(const VecI &tosub, VecI &out) {
    if (tosub._n == _n) {
        int *_tmparr = new int[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] - tosub[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}

void VecI::mul(const VecI &tomul, VecI &out) {
    if (tomul._n == _n) {
        int *_tmparr = new int[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] * tomul[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}


void VecI::div(const VecI &todiv, VecI &out) {
    if (todiv._n == _n) {
        int *_tmparr = new int[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] / todiv[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}

void VecI::square_root() {
    int *me = (int*)(*this);
    for (int i = 0; i < _n; ++i) {
        me[i] = (int)sqrt((double)me[i]);
    }
}

int VecI::sum() {
    int *me = (int*)(*this);
    int sum = 0;
    for( int n = 0; n < _n; n++) {
        sum += me[n];
    }
    return sum;
}

char * VecI::class_name() {
    char *name = new char[7];
    strcpy(name, "VecI");
    return name;
}


void VecI::abs_val() {
    for (int n = 0; n < _n; ++n) {
        if (_dat[n] < 0) { _dat[n] *= -1; }
    }
}


void VecI::std_normal() {
    // @TODO: would like avg and stdev to stay double, even for floats!
    (*this) -= (int)this->avg();
    double mean, stdev;
    this->sample_stats(mean, stdev);
    (*this) /= (int)stdev;
}


void VecI::remove(int index) {
    int *_tmp_arr = new int[_n - 1];
    int _cnt = 0;
    for (int i = 0; i < _n; ++i) {
        if (i != index) {
            _tmp_arr[_cnt] = _dat[i];
            ++_cnt;
        }
    }
    if (!_shallow) {
        delete []_dat;
    }
    _n = _n - 1;
    _dat = _tmp_arr;
    _shallow = false;
}

// Could be faster for ints
int VecI::intCompare( const void *a, const void *b ) {
    int c = *(int *)a - *(int *)b;
    if ( c < 0 ) return -1;
    if ( c > 0 ) return 1;
    return 0;
}

void VecI::sort() {
    qsort(_dat, _n, sizeof(int), intCompare);
}

int VecI::index(int val) {
    for (int i = 0; i < _n; ++i) {
        if (val == _dat[i]) {
            return i;
        }
    }
    return -1;
}

double VecI::avg() const {
    double total = 0;
    for( int n = 0; n < _n; ++n) {
        total += _dat[n];
    }
    return total/_n;
}

void VecI::sample_stats(double &mean, double &std_dev) {
    // Raw score method (below) is faster (~1.5X) than deviation score method
    // commonly used
    int *me = this->pointer();
    double _sum = 0.0;
    double _val;
    double _sumSq = 0.0;
    int _len = this->dim();
    for( int i=0; i<_len; ++i ) {
        _val = (double)me[i];
        _sum += _val;
        _sumSq += _val *_val;
    }
    double tmp = _sumSq - ((_sum * _sum)/_len);
    tmp /= _len>1 ? _len-1 : 1;
#ifdef WIN32
    std_dev = sqrt( tmp );
#else
    std_dev = std::sqrt( tmp );
#endif
    mean = _sum/_len;
}

double VecI::_zScore(double mean, double sigma, double x) {
    return (x - mean)/(sigma == 0.0 ? 1E-20: sigma);
}

void VecI::mask_as_vec(int return_val, VecI &mask, VecI &out) {
    if (mask.size() != _n) { Rprintf("mask.size() != this->length()"); R_ShowMessage("Serious error in obiwarp."); }
    int *me = (int*)(*this);
    int *maskptr = (int*)(mask);
    int *tmparr = new int[_n];
    int newcnt = 0;
    for (int i = 0; i < _n; ++i) {
        if (maskptr[i] == return_val) {
            tmparr[newcnt] = me[i];
            ++newcnt;
        }
    }
    out.take(newcnt, tmparr);
}

void VecI::hist(int num_bins, VecD &bins, VecI &freqs) {
    int i;

    // Create the scaling factor
    int _min; int _max;
    min_max(_min, _max);
    double dmin = (double)_min;
    double conv = ((double)num_bins)/(double)(_max - _min);

    // initialize arrays
    VecD _bins(num_bins);
    VecI _freqs(num_bins, 0);
    int _len = this->dim();
    int *me = this->pointer();

    // Create the histogram:
    for (i = 0; i < _len; ++i) {
        int index = (int)((me[i]-_min)*conv);
        if (index == num_bins) {
            --index;
        }
        _freqs[index]++;
    }

    // Create the bins:
    double iconv = 1.0/conv;
    for (i = 0; i < num_bins; ++i) {
        _bins[i] = ((i+0.5) * iconv) + dmin;  // avg
        // _bins[i] = ((i+0.5) * iconv) + dmin; //min
    }

    bins.take(_bins);
    freqs.take(_freqs);
}

void VecI::logarithm(double base) {
    int *me = (int*)(*this);
    for (int i = 0; i < _n; ++i) {
        //printf("ME: %f\n", me[i]);
        me[i] = (int)(log((double)(me[i]))/log(base));
        //printf("MELOGGED: %f\n", me[i]);
    }
}

void VecI::min_max(int &mn, int &mx) {
    int *me = (int*)(*this);
    mn = me[0];
    mx = me[0];
    for (int n = 0; n < _n; ++n) {
        mn = min(mn, me[n]);
        mx = max(mx, me[n]);
    }
}



void VecI::print(bool without_length ) {
    if (!without_length) {
        //std::cout << _n << std::endl;
    }
    int i;
    for (i = 0; i < _n - 1; ++i) {
        //std::cout << _dat[i] << " ";
    }
    //std::cout << _dat[i]; // the last one
    //std::cout << std::endl;
}

void VecI::print(const char *filename, bool without_length) {
    std::ofstream fh(filename);
    if (!fh) {
        //std::cout << "Error opening file " << filename << std::endl;
    }
    this->print(fh, without_length);
    fh.close();
}

void VecI::print(std::ostream &fout, bool without_length) {
    int i;
    if (!without_length) {
        fout << _n << std::endl;
    }
    for (i = 0; i < _n - 1; ++i) {
        fout << _dat[i] << " ";
    }
    fout << _dat[i];
    fout << std::endl;
}


// Class functions:
// THIS MUST BE FOR int AND DOUBLE ONLY!!!
// This is a fairly precise Fortran->C translation of the SLATEC chim code
// Evaluate the deriv at each x point
// return 1 if less than 2 data points
// return 0 if no errors
// ASSUMES monotonicity of the X data points !!!!!
// ASSUMES that this->length() >= 2
// If length == 1 then derivs[0] is set to 0
// If length == 0 then prints message to STDERR and returns;
void VecI::chim(VecI &x, VecI &y, VecI &out_derivs) {
    int *tmp_derivs = new int[x.length()];

#ifdef JTP_BOUNDS_CHECK
    if (x.length() != y.length()) { Rprintf("x.length() != y.length()"); R_ShowMessage("Serious error in obiwarp.");}
#endif
    int length = x.length();
    int del1;
    int del2;
    int h1;
    int h2;
    int hsum;
    int w1;
    int w2;
    int dmax;
    int dmin;
    int three = (int)3.0;
    int dsave;
    int drat1;
    int drat2;
    int hsumt3;

    int ierr = 0;
    int lengthLess1 = length - 1;

    if (length < 2) {
        if (length == 1) {
            tmp_derivs[0] = 0;
            return;
        }
        else {
	  Rprintf("trying to chim with 0 data points!\n");
        }
    }

    h1 = x[1] - x[0];
    del1 = (y[1] - y[0]) / h1;
    dsave = del1;

    // special case length=2 --use linear interpolation
    if (lengthLess1 < 2) {
        tmp_derivs[0] = del1;
        tmp_derivs[1] = del1;
        out_derivs.take(3, tmp_derivs);
        return;
    }

    // Normal case (length >= 3)
// 10

    h2 = x[2] - x[1];
    del2 = (y[2] - y[1]) / h2;

// SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
//     SHAPE-PRESERVING.

    hsum = h1 + h2;
    w1 = (h1 + hsum)/hsum;
    w2 = -h1/hsum;
    tmp_derivs[0] = (w1*del1) + (w2*del2);
    if ( pchst(tmp_derivs[0],del1) <= 0 ) {
        tmp_derivs[0] = (int)0;
    }
    else if ( pchst(del1,del2) < 0 ) {
        // need to do this check only if monotonicity switches
        dmax = three * del1;
        if (abs(tmp_derivs[0]) > abs(dmax)) {
            tmp_derivs[0] = dmax;
        }
    }

    int pchstval;
    int ind;

    for (ind = 1; ind < lengthLess1; ind++) {
        if (ind != 1) {
            h1 = h2;
            h2 = x[ind+1] - x[ind];
            hsum = h1 + h2;
            del1 = del2;
            del2 = (y[ind+1] - y[ind])/h2;
        }
// 40
        tmp_derivs[ind] = (int)0;

        pchstval = pchst(del1,del2);

// 45
        if (pchstval > 0) {
            hsumt3 = hsum+hsum+hsum;
            w1 = (hsum + h1)/hsumt3;
            w2 = (hsum + h2)/hsumt3;
            dmax = (int)max( abs(del1), abs(del2) );
            dmin = (int)min( abs(del1), abs(del2) );
            drat1 = del1/dmax;
            drat2 = del2/dmax;
            tmp_derivs[ind] = dmin/(w1*drat1 + w2*drat2);
        }
// 42
        else if (pchstval < 0 ) {
            ierr = ierr + 1;
            dsave = del2;
            continue;
        }
// 41
        else {  // equal to zero
            if (del2 == (int)0) { continue; }
            if (VecI::pchst(dsave,del2) < 0) { ierr = ierr + 1; }
            dsave = del2;
            continue;
        }
    }

// 50
    w1 = -h2/hsum;
    w2 = (h2 + hsum)/hsum;
    tmp_derivs[ind] = w1*del1 + w2*del2;
    if ( VecI::pchst(tmp_derivs[ind],del2) <= 0 ) {
        tmp_derivs[ind] = (int)0;
    }
    else if ( VecI::pchst(del1, del2) < 0) {
        // NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
        dmax = three*del2;
        if (abs(tmp_derivs[ind]) > abs(dmax)) {
            tmp_derivs[ind] = dmax;
        }
    }
    out_derivs.take(length, tmp_derivs);
    return;
}


void VecI::xy_to_x(VecI &x, VecI &y) {
    int *_x = (int*)x;
    int *_y = (int*)y;
    for (int i = 0; i < x.length(); i++) {
        _y[i] = _y[i] - _x[i];
    }
}

void VecI::x_to_xy(VecI &x, VecI &y) {
    int *_x = (int*)x;
    int *_y = (int*)y;
    for (int i = 0; i < x.length(); i++) {
        _y[i] = _y[i] + _x[i];
    }
}


void VecI::linear_derivs(VecI &x, VecI &y, VecI &out_derivs) {
    VecI tmp_d(x.size());
    for (int i = 0; i < x.size(); ++i) {
        tmp_d[i] = (y[i+1] - y[i]) / (x[i+1]-x[i]);
    }
    out_derivs.take(tmp_d);
}


void VecI::linear_interp(VecI &xin, VecI &yin, VecI &xe, VecI &out_ye, int sorted) {

    if (out_ye.size() == 0) {
        int *to_take = new int[xe.size()];
        out_ye.take(xe.size(), to_take);
    }

    // Calc the derivs:
    VecI derivs;
    VecI::linear_derivs(xin,yin,derivs);
    int i,j,ir;  // i indexes xin, j indexes xnew
    int ifirst = 0;

    // find the bounding points in xin
    int istart;
    int dt;

    if (sorted) {
        istart = 0;
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            istart = i;
            dt = xe[j] - xin[ifirst];  // diff in x, eval to input
            out_ye[j] = yin[ifirst] + (dt*derivs[ifirst]);
        }
    }
    else {

        // find the bounding points in xin
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            istart = 0;
            // @TODO: This should be a binary search:
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            dt = xe[j] - xin[ifirst];  // diff in x, eval to input
            out_ye[j] = yin[ifirst] + (dt * ((yin[ir] - yin[ifirst]) / (xin[ir]-xin[ifirst])) );
        }
    }
}

int VecI::sum_of_sq() {
    int *me = this->pointer();
    int total = 0;
    for( int n = 0; n < this->size(); n++) {
        total += me[n]*me[n];
    }
    return total;
}


double VecI::pearsons_r(VecI &x, VecI &y) {

    // Preparation:
    double sum_xTy = VecI::dot_product(x,y);
    double sum_x = x.sum();
    double sum_y = y.sum();
    // Could this step be sped up?
    double sum_x2 = x.sum_of_sq();
    double sum_y2 = y.sum_of_sq();
    int N = x.dim();

    // Here it is:
    // 'E' is Capital Sigma
    // r = EXY - (EXEY/N)
    //    -----------------
    //    sqrt( (EX^2 - (EX)^2/N) * (EY^2 - (EY)^2/N) )

    double top = sum_xTy - ((sum_x * sum_y)/N);
    double fbot = sum_x2 - ((sum_x*sum_x)/N);  //first part of bottom
    double sbot = sum_y2 - ((sum_y*sum_y)/N);  //second part of bottom
    return top / sqrt(fbot * sbot);

}

double VecI::covariance(VecI &x, VecI &y) {
    int i;
    int len = x.size();
    double mean_x = 0;
    double mean_y = 0;
    // get the means and x * y
    for (i = 0; i < len; ++i) {
        mean_x += x[i];
        mean_y += y[i];
    }
    mean_x /= len;
    mean_y /= len;
    double cov = 0;
    for (i = 0; i < len; ++i) {
        cov += (x[i] - mean_x) * (y[i] - mean_y);
    }
    return cov/len;
}

double VecI::euclidean(VecI &x, VecI &y) {
    VecF diff(x.size());
    double sum_of_diffs = 0;
    for (int i = 0; i < x.size(); ++i) {
        sum_of_diffs += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sqrt(sum_of_diffs);
}

int VecI::dot_product(VecI &x, VecI &y) {
    //assert(x.dim() == y.dim());
    int sum = 0;
    for (int i = 0; i < x.dim(); i++) {
        sum += (x[i] * y[i]);
    }
    return sum;
}

void VecI::chfe(VecI &xin, VecI &yin, VecI &xe, VecI &out_ye, int sorted) {
    //xin.print(); yin.print();

    if (out_ye.size() == 0) {
        int *to_take = new int[xe.size()];
        out_ye.take(xe.size(), to_take);
    }

    // Calc the derivs:
    VecI derivs;
    VecI::chim(xin,yin,derivs);
    int i,j,ir;  // i indexes xin, j indexes xnew
    int ifirst = 0;

    // find the bounding points in xin
    int istart;


    if (sorted) {
        VecI c2(xin.size());
        VecI c3(xin.size());
        calc_cubic_coeff(xin, yin, derivs, c2, c3);
        istart = 0;
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            istart = i;

            chfev(xin[ifirst], yin[ifirst], derivs[ifirst], c2[ifirst], c3[ifirst], xe[j], out_ye[j]);
        }
    }
    else {

        // find the bounding points in xin
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            istart = 0;
            // @TODO: This should be a binary search:
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }

            chfev_all(xin[ifirst], xin[ir], yin[ifirst], yin[ir], derivs[ifirst], derivs[ir], xe[j], out_ye[j]);
        }
    }
}




void VecI::calc_cubic_coeff(VecI &x, VecI &y, VecI &derivs, VecI &c2, VecI &c3) {

    //  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
    int DEL1, DEL2, DELTA, H;
    for (int i = 0; i < x.size() - 1; ++i) {
        H = x[i+1] - x[i];
        DELTA = (y[i+1] - y[i])/H;
        DEL1 = (derivs[i] - DELTA)/H;
        DEL2 = (derivs[i+1] - DELTA)/H;
        c2[i] = -(DEL1+DEL1 + DEL2);
        c3[i] = (DEL1 + DEL2)/H;
    }

}

void VecI::chfev_all(int X1, int X2, int F1, int F2, int D1, int D2, int XE, int &FE) {
    int C2, C3, DEL1, DEL2, DELTA, H, X;

    H = X2 - X1;

    //  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
    DELTA = (F2 - F1)/H;
    DEL1 = (D1 - DELTA)/H;
    DEL2 = (D2 - DELTA)/H;
    C2 = -(DEL1+DEL1 + DEL2);
    C3 = (DEL1 + DEL2)/H;

    X = XE - X1;

    FE = F1 + X*(D1 + X*(C2 + X*C3));
}


void VecI::chfev(int X1, int F1, int D1, int C2, int C3, int XE, int &FE) {
    int X;
    X = XE - X1;

    FE = F1 + X*(D1 + X*(C2 + X*C3));
}


void VecI::chfe_xy(VecI &x, VecI &y, VecI &new_x, VecI &out_new_y, int sorted) {
    VecI::xy_to_x(x,y);
    chfe(x,y,new_x, out_new_y, sorted);
    x_to_xy(new_x, out_new_y);
    VecI::x_to_xy(x,y);
}

double VecI::sum_sq_res_yeqx(VecI &x, VecI &y) {
    double __sum = 0.0;
    for (int i = 0; i < x.length(); ++i) {
        int diff = x[i] - y[i];
        __sum += 0.5*(diff*diff);
    }
    return __sum;
}

double VecI::avg_sq_res_yeqx(VecI &x, VecI &y) {
    return (sum_sq_res_yeqx(x,y))/x.length();
}

double VecI::avg_abs_diff(VecI &x, VecI &y) {
    double sum = 0.0;
    for (int i = 0; i < x.length(); ++i) {
      sum += (double)abs(x[i] - y[i]);
    }
    return sum/x.length();
}


void VecI::rsq_slope_intercept(VecI &x, VecI &y, double &rsq, double &slope, double &y_intercept) {
    int i;
    double mean_x = x.avg();
    double mean_y = y.avg();
    double sum_sq_res_xx = 0.0;
    double sum_sq_res_yy = 0.0;
    double sum_sq_res_xy = 0.0;

    for (i = 0; i < x.length(); ++i) {
        double x_minus_mean_i, y_minus_mean_i;
        x_minus_mean_i = ( (double)(x[i]) ) - mean_x;
        y_minus_mean_i = ( (double)(y[i]) ) - mean_y;
        sum_sq_res_xx += x_minus_mean_i*x_minus_mean_i;
        sum_sq_res_yy += y_minus_mean_i*y_minus_mean_i;
        sum_sq_res_xy += x_minus_mean_i*y_minus_mean_i;
    }
    slope = sum_sq_res_xy/sum_sq_res_xx;
    y_intercept = mean_y - (slope * mean_x);
    rsq = (sum_sq_res_xy*sum_sq_res_xy)/(sum_sq_res_xx*sum_sq_res_yy);
}





/****************************************************************
 * VecD
 ***************************************************************/

// Constructors:
VecD::VecD() : _n(0), _shallow(true) {
#ifdef JTP_DEBUG
    Rprintf("Creating DATA (NO ARGS)");
#endif
}

VecD::VecD(int n) : _n(n), _shallow(false) {
#ifdef JTP_BOUNDS_CHECK
    if (n < 0) { Rprintf("n < 0, exiting"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
    _dat = new double[_n];
#ifdef JTP_DEBUG
    Rprintf("Creating DATA(N)");
#endif
}

VecD::VecD(int n, const double &val) : _n(n), _shallow(false) {
    _dat = new double[_n];
    for (int i = 0; i < _n; ++i) {
        _dat[i] = val;
    }
#ifdef JTP_DEBUG
    Rprintf("Creating DATA(N,double)");
#endif
}

VecD::VecD(int n, double *arr, bool shallow) : _n(n), _dat(arr), _shallow(shallow) {
#ifdef JTP_DEBUG
    Rprintf("SHALLOW, (N,*ARR)");
#endif
}

VecD::VecD(const VecD &A, bool shallow) : _n(A._n), _shallow(shallow) {
    if (!shallow) {
        _dat = new double[_n];
        for (int i = 0; i < _n; ++i) {
            _dat[i] = A._dat[i];
        }
    }
    else {
        _dat = A._dat;
    }
#ifdef JTP_DEBUG
    Rprintf("created with VecD(const VecD &A)");
#endif
}

void VecD::to_f(VecF &out) {
    VecF _tmp(_n);
    for (int i = 0; i < _n; ++i) {
        _tmp[i] = (float)(_dat[i]);
    }
    out.take(_tmp);
}

void VecD::to_i(VecI &out) {
    VecI _tmp(_n);
    for (int i = 0; i < _n; ++i) {
        _tmp[i] = (int)(_dat[i]);
    }
    out.take(_tmp);
}


void VecD::set(int n, double *arr) {
    if (!_shallow) {
        delete[] _dat;
    }
    _dat = arr;
    _shallow = true;
    _n = n;
}

void VecD::take(int n, double *arr) {
    if (!_shallow) {
        delete[] _dat;
    }
    _dat = arr;
    _shallow = false;
    _n = n;
}

void VecD::set(VecD &A) {
    if (!_shallow) {
        delete[] _dat;
    }
    _dat = A._dat;
    _shallow = true;
    _n = A._n;
}


void VecD::take(VecD &A) {
    if (!_shallow) {
        delete[] _dat;
    }
    if (A._shallow) {
        Rprintf("Can't take ownership of memory of a shallow Vec!");
        R_ShowMessage("Serious error in obiwarp.");
    }
    _dat = A._dat;
    A._shallow = true;
    _shallow = false;
    _n = A._n;
}


bool VecD::operator==(const VecD &A) {
    // We don't care if one is shallow and the A is not!
    if (A._n == _n) { // Same size
        if (A._dat == _dat) { return true; }  // Same data
        else {
            for (int i = 0; i < _n; ++i) {
                if (A._dat[i] != _dat[i]) { return false; }
            }
            return true;
        }
    }
    else {
        return false;
    }
}

void VecD::copy(VecD &receiver, bool shallow) const {
    if (!receiver._shallow) {
        delete[] receiver._dat;
    }
    if (shallow) {
        receiver._dat = _dat;
        receiver._n = _n;
        receiver._shallow = true;
    }
    else {
        receiver._n = _n;
        receiver._dat = new double[_n];
        _copy(receiver._dat, _dat, _n);
        receiver._shallow = false;
    }
}

void VecD::_copy(double *p1, const double *p2, int len) const {
    // This is slightly faster on gcc Linux Mandrake
    for (int i = 0; i < len; ++i) {
        p1[i] = p2[i];
    }
}

VecD & VecD::operator=(const double &val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] = val;
    }
    return *this;
}

VecD & VecD::operator=(VecD &A) {
#ifdef JTP_DEBUG
    Rprintf("IN ASSIGNMENT OP tOP");
#endif
    if (this != &A) {
#ifdef JTP_DEBUG
        Rprintf("IN ASSIGNMENT OP MID");
#endif
        if (!_shallow) {
            delete[] _dat;
        }
        _n = A._n;
        _dat = new double[_n];
        _copy(_dat, A._dat, _n);
        _shallow = false;
    }
    return *this;
}

VecD::~VecD( ) {
#ifdef JTP_DEBUG
    Rprintf("DESTRUCTOR");
#endif
    if (!_shallow) {
#ifdef JTP_DEBUG
        Rprintf("DELETING DATA");
#endif
        delete[] _dat;
    }
}

/*************************
 * MATH OPERATORS
 ************************/

void VecD::operator+=(double val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] += val;
    }
}
void VecD::operator-=(double val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] -= val;
    }
}
void VecD::operator*=(double val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] *= val;
    }
}
void VecD::operator/=(double val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] /= val;
    }
}
void VecD::operator+=(const VecD &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] += A[i];
    }
}

void VecD::operator-=(const VecD &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] -= A[i];
    }
}

void VecD::operator*=(const VecD &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] *= A[i];
    }
}

void VecD::operator/=(const VecD &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] /= A[i];
    }
}

void VecD::add(const VecD &toadd, VecD &out) {
    if (toadd._n == _n) {
        double *_tmparr = new double[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] + toadd[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}

void VecD::sub(const VecD &tosub, VecD &out) {
    if (tosub._n == _n) {
        double *_tmparr = new double[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] - tosub[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}

void VecD::mul(const VecD &tomul, VecD &out) {
    if (tomul._n == _n) {
        double *_tmparr = new double[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] * tomul[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}


void VecD::div(const VecD &todiv, VecD &out) {
    if (todiv._n == _n) {
        double *_tmparr = new double[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] / todiv[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}

void VecD::square_root() {
    double *me = (double*)(*this);
    for (int i = 0; i < _n; ++i) {
        me[i] = (double)sqrt((double)me[i]);
    }
}


double VecD::sum() {
    double *me = (double*)(*this);
    double sum = 0;
    for( int n = 0; n < _n; n++) {
        sum += me[n];
    }
    return sum;
}

char * VecD::class_name() {
    char *name = new char[7];
    strcpy(name, "VecD");
    return name;
}


void VecD::abs_val() {
    for (int n = 0; n < _n; ++n) {
        if (_dat[n] < 0) { _dat[n] *= -1; }
    }
}


void VecD::std_normal() {
    // @TODO: would like avg and stdev to stay double, even for floats!
    (*this) -= (double)this->avg();
    double mean, stdev;
    this->sample_stats(mean, stdev);
    (*this) /= (double)stdev;
}


void VecD::remove(int index) {
    double *_tmp_arr = new double[_n - 1];
    int _cnt = 0;
    for (int i = 0; i < _n; ++i) {
        if (i != index) {
            _tmp_arr[_cnt] = _dat[i];
            ++_cnt;
        }
    }
    if (!_shallow) {
        delete []_dat;
    }
    _n = _n - 1;
    _dat = _tmp_arr;
    _shallow = false;
}

// Could be faster for ints
int VecD::doubleCompare( const void *a, const void *b ) {
    double c = *(double *)a - *(double *)b;
    if ( c < 0 ) return -1;
    if ( c > 0 ) return 1;
    return 0;
}

void VecD::sort() {
    qsort(_dat, _n, sizeof(double), doubleCompare);
}

int VecD::index(double val) {
    for (int i = 0; i < _n; ++i) {
        if (val == _dat[i]) {
            return i;
        }
    }
    return -1;
}

double VecD::avg() const {
    double total = 0;
    for( int n = 0; n < _n; ++n) {
        total += _dat[n];
    }
    return total/_n;
}

void VecD::sample_stats(double &mean, double &std_dev) {
    // Raw score method (below) is faster (~1.5X) than deviation score method
    // commonly used
    double *me = this->pointer();
    double _sum = 0.0;
    double _val;
    double _sumSq = 0.0;
    int _len = this->dim();
    for( int i=0; i<_len; ++i ) {
        _val = (double)me[i];
        _sum += _val;
        _sumSq += _val *_val;
    }
    double tmp = _sumSq - ((_sum * _sum)/_len);
    tmp /= _len>1 ? _len-1 : 1;
#ifdef WIN32
    std_dev = sqrt( tmp );
#else
    std_dev = std::sqrt( tmp );
#endif
    mean = _sum/_len;
}

double VecD::_zScore(double mean, double sigma, double x) {
    return (x - mean)/(sigma == 0.0 ? 1E-20: sigma);
}

void VecD::mask_as_vec(double return_val, VecI &mask, VecD &out) {
    if (mask.size() != _n) { Rprintf("mask.size() != this->length()"); R_ShowMessage("Serious error in obiwarp.");}
    double *me = (double*)(*this);
    int *maskptr = (int*)(mask);
    double *tmparr = new double[_n];
    int newcnt = 0;
    for (int i = 0; i < _n; ++i) {
        if (maskptr[i] == return_val) {
            tmparr[newcnt] = me[i];
            ++newcnt;
        }
    }
    out.take(newcnt, tmparr);
}

void VecD::hist(int num_bins, VecD &bins, VecI &freqs) {
    int i;

    // Create the scaling factor
    double _min; double _max;
    min_max(_min, _max);
    double dmin = (double)_min;
    double conv = ((double)num_bins)/(double)(_max - _min);

    // initialize arrays
    VecD _bins(num_bins);
    VecI _freqs(num_bins, 0);
    int _len = this->dim();
    double *me = this->pointer();

    // Create the histogram:
    for (i = 0; i < _len; ++i) {
        int index = (int)((me[i]-_min)*conv);
        if (index == num_bins) {
            --index;
        }
        _freqs[index]++;
    }

    // Create the bins:
    double iconv = 1.0/conv;
    for (i = 0; i < num_bins; ++i) {
        _bins[i] = ((i+0.5) * iconv) + dmin;  // avg
        // _bins[i] = ((i+0.5) * iconv) + dmin; //min
    }

    bins.take(_bins);
    freqs.take(_freqs);
}

void VecD::logarithm(double base) {
    double *me = (double*)(*this);
    for (int i = 0; i < _n; ++i) {
        //printf("ME: %f\n", me[i]);
        me[i] = (double)(log((double)(me[i]))/log(base));
        //printf("MELOGGED: %f\n", me[i]);
    }
}

void VecD::min_max(double &mn, double &mx) {
    double *me = (double*)(*this);
    mn = me[0];
    mx = me[0];
    for (int n = 0; n < _n; ++n) {
        mn = min(mn, me[n]);
        mx = max(mx, me[n]);
    }
}



void VecD::print(bool without_length ) {
    if (!without_length) {
        //std::cout << _n << std::endl;
    }
    int i;
    for (i = 0; i < _n - 1; ++i) {
        //std::cout << _dat[i] << " ";
    }
    //std::cout << _dat[i]; // the last one
    //std::cout << std::endl;
}

void VecD::print(const char *filename, bool without_length) {
    std::ofstream fh(filename);
    if (!fh) {
        //std::cout << "Error opening file " << filename << std::endl;
    }
    this->print(fh, without_length);
    fh.close();
}

void VecD::print(std::ostream &fout, bool without_length) {
    int i;
    if (!without_length) {
        fout << _n << std::endl;
    }
    for (i = 0; i < _n - 1; ++i) {
        fout << _dat[i] << " ";
    }
    fout << _dat[i];
    fout << std::endl;
}


// Class functions:
// THIS MUST BE FOR double AND DOUBLE ONLY!!!
// This is a fairly precise Fortran->C translation of the SLATEC chim code
// Evaluate the deriv at each x point
// return 1 if less than 2 data points
// return 0 if no errors
// ASSUMES monotonicity of the X data points !!!!!
// ASSUMES that this->length() >= 2
// If length == 1 then derivs[0] is set to 0
// If length == 0 then prints message to STDERR and returns;
void VecD::chim(VecD &x, VecD &y, VecD &out_derivs) {
    double *tmp_derivs = new double[x.length()];

#ifdef JTP_BOUNDS_CHECK
    if (x.length() != y.length()) { Rprintf("x.length() != y.length()"); R_ShowMessage("Serious error in obiwarp.");}
#endif
    int length = x.length();
    double del1;
    double del2;
    double h1;
    double h2;
    double hsum;
    double w1;
    double w2;
    double dmax;
    double dmin;
    double three = (double)3.0;
    double dsave;
    double drat1;
    double drat2;
    double hsumt3;

    int ierr = 0;
    int lengthLess1 = length - 1;

    if (length < 2) {
        if (length == 1) {
            tmp_derivs[0] = 0;
            return;
        }
        else {
	  Rprintf("trying to chim with 0 data points!\n");
        }
    }

    h1 = x[1] - x[0];
    del1 = (y[1] - y[0]) / h1;
    dsave = del1;

    // special case length=2 --use linear interpolation
    if (lengthLess1 < 2) {
        tmp_derivs[0] = del1;
        tmp_derivs[1] = del1;
        out_derivs.take(3, tmp_derivs);
        return;
    }

    // Normal case (length >= 3)
// 10

    h2 = x[2] - x[1];
    del2 = (y[2] - y[1]) / h2;

// SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
//     SHAPE-PRESERVING.

    hsum = h1 + h2;
    w1 = (h1 + hsum)/hsum;
    w2 = -h1/hsum;
    tmp_derivs[0] = (w1*del1) + (w2*del2);
    if ( pchst(tmp_derivs[0],del1) <= 0 ) {
        tmp_derivs[0] = (double)0;
    }
    else if ( pchst(del1,del2) < 0 ) {
        // need to do this check only if monotonicity switches
        dmax = three * del1;
        if (fabs(tmp_derivs[0]) > fabs(dmax)) {
            tmp_derivs[0] = dmax;
        }
    }

    int pchstval;
    int ind;

    for (ind = 1; ind < lengthLess1; ind++) {
        if (ind != 1) {
            h1 = h2;
            h2 = x[ind+1] - x[ind];
            hsum = h1 + h2;
            del1 = del2;
            del2 = (y[ind+1] - y[ind])/h2;
        }
// 40
        tmp_derivs[ind] = (double)0;

        pchstval = pchst(del1,del2);

// 45
        if (pchstval > 0) {
            hsumt3 = hsum+hsum+hsum;
            w1 = (hsum + h1)/hsumt3;
            w2 = (hsum + h2)/hsumt3;
            dmax = (double)max( fabs(del1), fabs(del2) );
            dmin = (double)min( fabs(del1), fabs(del2) );
            drat1 = del1/dmax;
            drat2 = del2/dmax;
            tmp_derivs[ind] = dmin/(w1*drat1 + w2*drat2);
        }
// 42
        else if (pchstval < 0 ) {
            ierr = ierr + 1;
            dsave = del2;
            continue;
        }
// 41
        else {  // equal to zero
            if (del2 == (double)0) { continue; }
            if (VecD::pchst(dsave,del2) < 0) { ierr = ierr + 1; }
            dsave = del2;
            continue;
        }
    }

// 50
    w1 = -h2/hsum;
    w2 = (h2 + hsum)/hsum;
    tmp_derivs[ind] = w1*del1 + w2*del2;
    if ( VecD::pchst(tmp_derivs[ind],del2) <= 0 ) {
        tmp_derivs[ind] = (double)0;
    }
    else if ( VecD::pchst(del1, del2) < 0) {
        // NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
        dmax = three*del2;
        if (fabs(tmp_derivs[ind]) > fabs(dmax)) {
            tmp_derivs[ind] = dmax;
        }
    }
    out_derivs.take(length, tmp_derivs);
    return;
}


void VecD::xy_to_x(VecD &x, VecD &y) {
    double *_x = (double*)x;
    double *_y = (double*)y;
    for (int i = 0; i < x.length(); i++) {
        _y[i] = _y[i] - _x[i];
    }
}

void VecD::x_to_xy(VecD &x, VecD &y) {
    double *_x = (double*)x;
    double *_y = (double*)y;
    for (int i = 0; i < x.length(); i++) {
        _y[i] = _y[i] + _x[i];
    }
}


void VecD::linear_derivs(VecD &x, VecD &y, VecD &out_derivs) {
    VecD tmp_d(x.size());
    for (int i = 0; i < x.size(); ++i) {
        tmp_d[i] = (y[i+1] - y[i]) / (x[i+1]-x[i]);
    }
    out_derivs.take(tmp_d);
}


void VecD::linear_interp(VecD &xin, VecD &yin, VecD &xe, VecD &out_ye, int sorted) {

    if (out_ye.size() == 0) {
        double *to_take = new double[xe.size()];
        out_ye.take(xe.size(), to_take);
    }

    // Calc the derivs:
    VecD derivs;
    VecD::linear_derivs(xin,yin,derivs);
    int i,j,ir;  // i indexes xin, j indexes xnew
    int ifirst = 0;

    // find the bounding points in xin
    int istart;
    double dt;

    if (sorted) {
        istart = 0;
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            istart = i;
            dt = xe[j] - xin[ifirst];  // diff in x, eval to input
            out_ye[j] = yin[ifirst] + (dt*derivs[ifirst]);
        }
    }
    else {

        // find the bounding points in xin
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            istart = 0;
            // @TODO: This should be a binary search:
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            dt = xe[j] - xin[ifirst];  // diff in x, eval to input
            out_ye[j] = yin[ifirst] + (dt * ((yin[ir] - yin[ifirst]) / (xin[ir]-xin[ifirst])) );
        }
    }
}

double VecD::sum_of_sq() {
    double *me = this->pointer();
    double total = 0;
    for( int n = 0; n < this->size(); n++) {
        total += me[n]*me[n];
    }
    return total;
}


double VecD::pearsons_r(VecD &x, VecD &y) {

    // Preparation:
    double sum_xTy = VecD::dot_product(x,y);
    double sum_x = x.sum();
    double sum_y = y.sum();
    // Could this step be sped up?
    double sum_x2 = x.sum_of_sq();
    double sum_y2 = y.sum_of_sq();
    int N = x.dim();

    // Here it is:
    // 'E' is Capital Sigma
    // r = EXY - (EXEY/N)
    //    -----------------
    //    sqrt( (EX^2 - (EX)^2/N) * (EY^2 - (EY)^2/N) )

    double top = sum_xTy - ((sum_x * sum_y)/N);
    double fbot = sum_x2 - ((sum_x*sum_x)/N);  //first part of bottom
    double sbot = sum_y2 - ((sum_y*sum_y)/N);  //second part of bottom
    return top / sqrt(fbot * sbot);

}

double VecD::covariance(VecD &x, VecD &y) {
    int i;
    int len = x.size();
    double mean_x = 0;
    double mean_y = 0;
    // get the means and x * y
    for (i = 0; i < len; ++i) {
        mean_x += x[i];
        mean_y += y[i];
    }
    mean_x /= len;
    mean_y /= len;
    double cov = 0;
    for (i = 0; i < len; ++i) {
        cov += (x[i] - mean_x) * (y[i] - mean_y);
    }
    return cov/len;
}

double VecD::euclidean(VecD &x, VecD &y) {
    VecF diff(x.size());
    double sum_of_diffs = 0;
    for (int i = 0; i < x.size(); ++i) {
        sum_of_diffs += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sqrt(sum_of_diffs);
}

double VecD::dot_product(VecD &x, VecD &y) {
    //assert(x.dim() == y.dim());
    double sum = 0;
    for (int i = 0; i < x.dim(); i++) {
        sum += (x[i] * y[i]);
    }
    return sum;
}

void VecD::chfe(VecD &xin, VecD &yin, VecD &xe, VecD &out_ye, int sorted) {
    //xin.print(); yin.print();

    if (out_ye.size() == 0) {
        double *to_take = new double[xe.size()];
        out_ye.take(xe.size(), to_take);
    }

    // Calc the derivs:
    VecD derivs;
    VecD::chim(xin,yin,derivs);
    int i,j,ir;  // i indexes xin, j indexes xnew
    int ifirst = 0;

    // find the bounding points in xin
    int istart;


    if (sorted) {
        VecD c2(xin.size());
        VecD c3(xin.size());
        calc_cubic_coeff(xin, yin, derivs, c2, c3);
        istart = 0;
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            istart = i;

            chfev(xin[ifirst], yin[ifirst], derivs[ifirst], c2[ifirst], c3[ifirst], xe[j], out_ye[j]);
        }
    }
    else {

        // find the bounding points in xin
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            istart = 0;
            // @TODO: This should be a binary search:
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }

            chfev_all(xin[ifirst], xin[ir], yin[ifirst], yin[ir], derivs[ifirst], derivs[ir], xe[j], out_ye[j]);
        }
    }
}




void VecD::calc_cubic_coeff(VecD &x, VecD &y, VecD &derivs, VecD &c2, VecD &c3) {

    //  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
    double DEL1, DEL2, DELTA, H;
    for (int i = 0; i < x.size() - 1; ++i) {
        H = x[i+1] - x[i];
        DELTA = (y[i+1] - y[i])/H;
        DEL1 = (derivs[i] - DELTA)/H;
        DEL2 = (derivs[i+1] - DELTA)/H;
        c2[i] = -(DEL1+DEL1 + DEL2);
        c3[i] = (DEL1 + DEL2)/H;
    }

}

void VecD::chfev_all(double X1, double X2, double F1, double F2, double D1, double D2, double XE, double &FE) {
    double C2, C3, DEL1, DEL2, DELTA, H, X;

    H = X2 - X1;

    //  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
    DELTA = (F2 - F1)/H;
    DEL1 = (D1 - DELTA)/H;
    DEL2 = (D2 - DELTA)/H;
    C2 = -(DEL1+DEL1 + DEL2);
    C3 = (DEL1 + DEL2)/H;

    X = XE - X1;

    FE = F1 + X*(D1 + X*(C2 + X*C3));
}


void VecD::chfev(double X1, double F1, double D1, double C2, double C3, double XE, double &FE) {
    double X;
    X = XE - X1;

    FE = F1 + X*(D1 + X*(C2 + X*C3));
}


void VecD::chfe_xy(VecD &x, VecD &y, VecD &new_x, VecD &out_new_y, int sorted) {
    VecD::xy_to_x(x,y);
    chfe(x,y,new_x, out_new_y, sorted);
    x_to_xy(new_x, out_new_y);
    VecD::x_to_xy(x,y);
}

double VecD::sum_sq_res_yeqx(VecD &x, VecD &y) {
    //// PDL way
    //return sum(0.5*(($y - $x)**2));
    double __sum = 0.0;
    for (int i = 0; i < x.length(); ++i) {
        double diff = x[i] - y[i];
        __sum += 0.5*(diff*diff);
    }
    return __sum;
}

double VecD::avg_sq_res_yeqx(VecD &x, VecD &y) {
    return (sum_sq_res_yeqx(x,y))/x.length();
}

double VecD::avg_abs_diff(VecD &x, VecD &y) {
    double sum = 0.0;
    for (int i = 0; i < x.length(); ++i) {
        sum += fabs(x[i] - y[i]);
    }
    return sum/x.length();
}


void VecD::rsq_slope_intercept(VecD &x, VecD &y, double &rsq, double &slope, double &y_intercept) {
    int i;
    double mean_x = x.avg();
    double mean_y = y.avg();
    double sum_sq_res_xx = 0.0;
    double sum_sq_res_yy = 0.0;
    double sum_sq_res_xy = 0.0;

    for (i = 0; i < x.length(); ++i) {
        double x_minus_mean_i, y_minus_mean_i;
        x_minus_mean_i = ( (double)(x[i]) ) - mean_x;
        y_minus_mean_i = ( (double)(y[i]) ) - mean_y;
        sum_sq_res_xx += x_minus_mean_i*x_minus_mean_i;
        sum_sq_res_yy += y_minus_mean_i*y_minus_mean_i;
        sum_sq_res_xy += x_minus_mean_i*y_minus_mean_i;
    }
    slope = sum_sq_res_xy/sum_sq_res_xx;
    y_intercept = mean_y - (slope * mean_x);
    rsq = (sum_sq_res_xy*sum_sq_res_xy)/(sum_sq_res_xx*sum_sq_res_yy);
}





/****************************************************************
 * VecF
 ***************************************************************/

// Constructors:
VecF::VecF() : _n(0), _shallow(true) {
#ifdef JTP_DEBUG
    Rprintf("Creating DATA (NO ARGS)");
#endif
}

VecF::VecF(int n) : _n(n), _shallow(false) {
#ifdef JTP_BOUNDS_CHECK
    if (n < 0) { Rprintf("n < 0, exiting"); R_ShowMessage("Serious error in obiwarp.");}
#endif
    _dat = new float[_n];
#ifdef JTP_DEBUG
    Rprintf("Creating DATA(N)");
#endif
}

VecF::VecF(int n, const float &val) : _n(n), _shallow(false) {
    _dat = new float[_n];
    for (int i = 0; i < _n; ++i) {
        _dat[i] = val;
    }
#ifdef JTP_DEBUG
    Rprintf("Creating DATA(N,float)");
#endif
}

VecF::VecF(int n, float *arr, bool shallow) : _n(n), _dat(arr), _shallow(shallow) {
#ifdef JTP_DEBUG
    Rprintf("SHALLOW, (N,*ARR)");
#endif
}

VecF::VecF(const VecF &A, bool shallow) : _n(A._n), _shallow(shallow) {
    if (!shallow) {
        _dat = new float[_n];
        for (int i = 0; i < _n; ++i) {
            _dat[i] = A._dat[i];
        }
    }
    else {
        _dat = A._dat;
    }
#ifdef JTP_DEBUG
    Rprintf("created with VecF(const VecF &A)");
#endif
}

void VecF::to_f(VecF &out) {
    VecF _tmp(_n);
    for (int i = 0; i < _n; ++i) {
        _tmp[i] = (float)(_dat[i]);
    }
    out.take(_tmp);
}

void VecF::to_i(VecI &out) {
    VecI _tmp(_n);
    for (int i = 0; i < _n; ++i) {
        _tmp[i] = (int)(_dat[i]);
    }
    out.take(_tmp);
}


void VecF::set(int n, float *arr) {
    if (!_shallow) {
        delete[] _dat;
    }
    _dat = arr;
    _shallow = true;
    _n = n;
}

void VecF::take(int n, float *arr) {
    if (!_shallow) {
        delete[] _dat;
    }
    _dat = arr;
    _shallow = false;
    _n = n;
}

void VecF::set(VecF &A) {
    if (!_shallow) {
        delete[] _dat;
    }
    _dat = A._dat;
    _shallow = true;
    _n = A._n;
}


void VecF::take(VecF &A) {
    if (!_shallow) {
        delete[] _dat;
    }
    if (A._shallow) {
        Rprintf("Can't take ownership of memory of a shallow Vec!");
        R_ShowMessage("Serious error in obiwarp.");
    }
    _dat = A._dat;
    A._shallow = true;
    _shallow = false;
    _n = A._n;
}


bool VecF::operator==(const VecF &A) {
    // We don't care if one is shallow and the A is not!
    if (A._n == _n) { // Same size
        if (A._dat == _dat) { return true; }  // Same data
        else {
            for (int i = 0; i < _n; ++i) {
                if (A._dat[i] != _dat[i]) { return false; }
            }
            return true;
        }
    }
    else {
        return false;
    }
}

void VecF::copy(VecF &receiver, bool shallow) const {
    if (!receiver._shallow) {
        delete[] receiver._dat;
    }
    if (shallow) {
        receiver._dat = _dat;
        receiver._n = _n;
        receiver._shallow = true;
    }
    else {
        receiver._n = _n;
        receiver._dat = new float[_n];
        _copy(receiver._dat, _dat, _n);
        receiver._shallow = false;
    }
}

void VecF::_copy(float *p1, const float *p2, int len) const {
    // This is slightly faster on gcc Linux Mandrake
    for (int i = 0; i < len; ++i) {
        p1[i] = p2[i];
    }
}

VecF & VecF::operator=(const float &val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] = val;
    }
    return *this;
}

VecF & VecF::operator=(VecF &A) {
#ifdef JTP_DEBUG
    Rprintf("IN ASSIGNMENT OP tOP");
#endif
    if (this != &A) {
#ifdef JTP_DEBUG
        Rprintf("IN ASSIGNMENT OP MID");
#endif
        if (!_shallow) {
            delete[] _dat;
        }
        _n = A._n;
        _dat = new float[_n];
        _copy(_dat, A._dat, _n);
        _shallow = false;
    }
    return *this;
}

VecF::~VecF( ) {
#ifdef JTP_DEBUG
    Rprintf("DESTRUCTOR");
#endif
    if (!_shallow) {
#ifdef JTP_DEBUG
        Rprintf("DELETING DATA");
#endif
        delete[] _dat;
    }
}

/*************************
 * MATH OPERATORS
 ************************/

void VecF::operator+=(float val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] += val;
    }
}
void VecF::operator-=(float val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] -= val;
    }
}
void VecF::operator*=(float val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] *= val;
    }
}
void VecF::operator/=(float val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] /= val;
    }
}
void VecF::operator+=(const VecF &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] += A[i];
    }
}

void VecF::operator-=(const VecF &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] -= A[i];
    }
}

void VecF::operator*=(const VecF &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] *= A[i];
    }
}

void VecF::operator/=(const VecF &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] /= A[i];
    }
}

void VecF::add(const VecF &toadd, VecF &out) {
    if (toadd._n == _n) {
        float *_tmparr = new float[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] + toadd[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}

void VecF::sub(const VecF &tosub, VecF &out) {
    if (tosub._n == _n) {
        float *_tmparr = new float[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] - tosub[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}

void VecF::mul(const VecF &tomul, VecF &out) {
    if (tomul._n == _n) {
        float *_tmparr = new float[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] * tomul[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}


void VecF::div(const VecF &todiv, VecF &out) {
    if (todiv._n == _n) {
        float *_tmparr = new float[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] / todiv[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}

void VecF::square_root() {
    float *me = (float*)(*this);
    for (int i = 0; i < _n; ++i) {
        me[i] = (float)sqrt((double)me[i]);
    }
}

float VecF::sum() {
    float *me = (float*)(*this);
    float sum = 0;
    for( int n = 0; n < _n; n++) {
        sum += me[n];
    }
    return sum;
}

char * VecF::class_name() {
    char *name = new char[7];
    strcpy(name, "VecF");
    return name;
}


void VecF::abs_val() {
    for (int n = 0; n < _n; ++n) {
        if (_dat[n] < 0) { _dat[n] *= -1; }
    }
}


void VecF::std_normal() {
    // @TODO: would like avg and stdev to stay double, even for floats!
    (*this) -= (float)this->avg();
    double mean, stdev;
    this->sample_stats(mean, stdev);
    (*this) /= (float)stdev;
}


void VecF::remove(int index) {
    float *_tmp_arr = new float[_n - 1];
    int _cnt = 0;
    for (int i = 0; i < _n; ++i) {
        if (i != index) {
            _tmp_arr[_cnt] = _dat[i];
            ++_cnt;
        }
    }
    if (!_shallow) {
        delete []_dat;
    }
    _n = _n - 1;
    _dat = _tmp_arr;
    _shallow = false;
}

// Could be faster for ints
int VecF::floatCompare( const void *a, const void *b ) {
    float c = *(float *)a - *(float *)b;
    if ( c < 0 ) return -1;
    if ( c > 0 ) return 1;
    return 0;
}

void VecF::sort() {
    qsort(_dat, _n, sizeof(float), floatCompare);
}

int VecF::index(float val) {
    for (int i = 0; i < _n; ++i) {
        if (val == _dat[i]) {
            return i;
        }
    }
    return -1;
}

double VecF::avg() const {
    double total = 0;
    for( int n = 0; n < _n; ++n) {
        total += _dat[n];
    }
    return total/_n;
}


void VecF::sample_stats(double &mean, double &std_dev) {
    // Raw score method (below) is faster (~1.5X) than deviation score method
    // commonly used
    float *me = this->pointer();
    double _sum = 0.0;
    double _val;
    double _sumSq = 0.0;
    int _len = this->dim();
    for( int i=0; i<_len; ++i ) {
        _val = (double)me[i];
        _sum += _val;
        _sumSq += _val *_val;
    }
    double tmp = _sumSq - ((_sum * _sum)/_len);
    tmp /= _len>1 ? _len-1 : 1;
#ifdef WIN32
    std_dev = sqrt( tmp );
#else
    std_dev = std::sqrt( tmp );
#endif
    mean = _sum/_len;
}

double VecF::_zScore(double mean, double sigma, double x) {
    return (x - mean)/(sigma == 0.0 ? 1E-20: sigma);
}

void VecF::mask_as_vec(float return_val, VecI &mask, VecF &out) {
    if (mask.size() != _n) { Rprintf("mask.size() != this->length()"); R_ShowMessage("Serious error in obiwarp.");}
    float *me = (float*)(*this);
    int *maskptr = (int*)(mask);
    float *tmparr = new float[_n];
    int newcnt = 0;
    for (int i = 0; i < _n; ++i) {
        if (maskptr[i] == return_val) {
            tmparr[newcnt] = me[i];
            ++newcnt;
        }
    }
    out.take(newcnt, tmparr);
}

void VecF::hist(int num_bins, VecD &bins, VecI &freqs) {
    int i;

    // Create the scaling factor
    float _min; float _max;
    min_max(_min, _max);
    double dmin = (double)_min;
    double conv = ((double)num_bins)/(double)(_max - _min);

    // initialize arrays
    VecD _bins(num_bins);
    VecI _freqs(num_bins, 0);
    int _len = this->dim();
    float *me = this->pointer();

    // Create the histogram:
    for (i = 0; i < _len; ++i) {
        int index = (int)((me[i]-_min)*conv);
        if (index == num_bins) {
            --index;
        }
        _freqs[index]++;
    }

    // Create the bins:
    double iconv = 1.0/conv;
    for (i = 0; i < num_bins; ++i) {
        _bins[i] = ((i+0.5) * iconv) + dmin;  // avg
        // _bins[i] = ((i+0.5) * iconv) + dmin; //min
    }

    bins.take(_bins);
    freqs.take(_freqs);
}

void VecF::logarithm(double base) {
    float *me = (float*)(*this);
    for (int i = 0; i < _n; ++i) {
        //printf("ME: %f\n", me[i]);
        me[i] = (float)(log((double)(me[i]))/log(base));
        //printf("MELOGGED: %f\n", me[i]);
    }
}

void VecF::min_max(float &mn, float &mx) {
    float *me = (float*)(*this);
    mn = me[0];
    mx = me[0];
    for (int n = 0; n < _n; ++n) {
        mn = min(mn, me[n]);
        mx = max(mx, me[n]);
    }
}



void VecF::print(bool without_length ) {
    if (!without_length) {
        //std::cout << _n << std::endl;
    }
    int i;
    for (i = 0; i < _n - 1; ++i) {
        //std::cout << _dat[i] << " ";
    }
    //std::cout << _dat[i]; // the last one
    //std::cout << std::endl;
}

void VecF::print(const char *filename, bool without_length) {
    std::ofstream fh(filename);
    if (!fh) {
        //std::cout << "Error opening file " << filename << std::endl;
    }
    this->print(fh, without_length);
    fh.close();
}

void VecF::print_tm() {
    int i;

    //std::cout << _n << std::endl;

    for (i = 0; i < _n - 1; ++i) {
        //std::cout << _dat[i] << " ";
    }
    //std::cout << _dat[i];
    //std::cout << std::endl;
}

void VecF::print(std::ostream &fout, bool without_length) {
    int i;
    if (!without_length) {
    fout <<"_n"<< _n << std::endl;
    }
    for (i = 0; i < _n - 1; ++i) {
        fout << _dat[i] << " ";
    }
    fout << _dat[i];
    fout << std::endl;
}




// Class functions:
// THIS MUST BE FOR float AND DOUBLE ONLY!!!
// This is a fairly precise Fortran->C translation of the SLATEC chim code
// Evaluate the deriv at each x point
// return 1 if less than 2 data points
// return 0 if no errors
// ASSUMES monotonicity of the X data points !!!!!
// ASSUMES that this->length() >= 2
// If length == 1 then derivs[0] is set to 0
// If length == 0 then prints message to STDERR and returns;
void VecF::chim(VecF &x, VecF &y, VecF &out_derivs) {
    float *tmp_derivs = new float[x.length()];

#ifdef JTP_BOUNDS_CHECK
    if (x.length() != y.length()) { Rprintf("x.length() != y.length()"); R_ShowMessage("Serious error in obiwarp.");}
#endif
    int length = x.length();
    float del1;
    float del2;
    float h1;
    float h2;
    float hsum;
    float w1;
    float w2;
    float dmax;
    float dmin;
    float three = (float)3.0;
    float dsave;
    float drat1;
    float drat2;
    float hsumt3;

    int ierr = 0;
    int lengthLess1 = length - 1;

    if (length < 2) {
        if (length == 1) {
            tmp_derivs[0] = 0;
            return;
        }
        else {
	  Rprintf("trying to chim with 0 data points!\n");
        }
    }

// THIS can be done BEFORE this routine if someone cares to...
//    // Check monotonicity
//    for (int i = 2; i < length; i++) {
//        if (x[i] <= x[i-1]) {
//            return 2;
//        }
//    }

    h1 = x[1] - x[0];
    del1 = (y[1] - y[0]) / h1;
    dsave = del1;

    // special case length=2 --use linear interpolation
    if (lengthLess1 < 2) {
        tmp_derivs[0] = del1;
        tmp_derivs[1] = del1;
        out_derivs.take(3, tmp_derivs);
        return;
    }

    // Normal case (length >= 3)
// 10

    h2 = x[2] - x[1];
    del2 = (y[2] - y[1]) / h2;

// SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
//     SHAPE-PRESERVING.

    hsum = h1 + h2;
    w1 = (h1 + hsum)/hsum;
    w2 = -h1/hsum;
    tmp_derivs[0] = (w1*del1) + (w2*del2);
    if ( pchst(tmp_derivs[0],del1) <= 0 ) {
        tmp_derivs[0] = (float)0;
    }
    else if ( pchst(del1,del2) < 0 ) {
        // need to do this check only if monotonicity switches
        dmax = three * del1;
        if (fabs(tmp_derivs[0]) > fabs(dmax)) {
            tmp_derivs[0] = dmax;
        }
    }

    int pchstval;
    int ind;

    for (ind = 1; ind < lengthLess1; ind++) {
        if (ind != 1) {
            h1 = h2;
            h2 = x[ind+1] - x[ind];
            hsum = h1 + h2;
            del1 = del2;
            del2 = (y[ind+1] - y[ind])/h2;
        }
// 40
        tmp_derivs[ind] = (float)0;

        pchstval = pchst(del1,del2);

// 45
        if (pchstval > 0) {
            hsumt3 = hsum+hsum+hsum;
            w1 = (hsum + h1)/hsumt3;
            w2 = (hsum + h2)/hsumt3;
            dmax = (float)max( fabs(del1), fabs(del2) );
            dmin = (float)min( fabs(del1), fabs(del2) );
            drat1 = del1/dmax;
            drat2 = del2/dmax;
            tmp_derivs[ind] = dmin/(w1*drat1 + w2*drat2);
        }
// 42
        else if (pchstval < 0 ) {
            ierr = ierr + 1;
            dsave = del2;
            continue;
        }
// 41
        else {  // equal to zero
            if (del2 == (float)0) { continue; }
            if (VecF::pchst(dsave,del2) < 0) { ierr = ierr + 1; }
            dsave = del2;
            continue;
        }
    }

// 50
    w1 = -h2/hsum;
    w2 = (h2 + hsum)/hsum;
    tmp_derivs[ind] = w1*del1 + w2*del2;
    if ( VecF::pchst(tmp_derivs[ind],del2) <= 0 ) {
        tmp_derivs[ind] = (float)0;
    }
    else if ( VecF::pchst(del1, del2) < 0) {
        // NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
        dmax = three*del2;
        if (fabs(tmp_derivs[ind]) > fabs(dmax)) {
            tmp_derivs[ind] = dmax;
        }
    }
    out_derivs.take(length, tmp_derivs);
    return;
}


void VecF::xy_to_x(VecF &x, VecF &y) {
    float *_x = (float*)x;
    float *_y = (float*)y;
    for (int i = 0; i < x.length(); i++) {
        _y[i] = _y[i] - _x[i];
    }
}

void VecF::x_to_xy(VecF &x, VecF &y) {
    float *_x = (float*)x;
    float *_y = (float*)y;
    for (int i = 0; i < x.length(); i++) {
        _y[i] = _y[i] + _x[i];
    }
}


void VecF::linear_derivs(VecF &x, VecF &y, VecF &out_derivs) {
    VecF tmp_d(x.size());
    for (int i = 0; i < x.size(); ++i) {
        tmp_d[i] = (y[i+1] - y[i]) / (x[i+1]-x[i]);
    }
    out_derivs.take(tmp_d);
}


void VecF::linear_interp(VecF &xin, VecF &yin, VecF &xe, VecF &out_ye, int sorted) {

    if (out_ye.size() == 0) {
        float *to_take = new float[xe.size()];
        out_ye.take(xe.size(), to_take);
    }

    // Calc the derivs:
    VecF derivs;
    VecF::linear_derivs(xin,yin,derivs);
    int i,j,ir;  // i indexes xin, j indexes xnew
    int ifirst = 0;

    // find the bounding points in xin
    int istart;
    float dt;

    if (sorted) {
        istart = 0;
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            istart = i;
            dt = xe[j] - xin[ifirst];  // diff in x, eval to input
            out_ye[j] = yin[ifirst] + (dt*derivs[ifirst]);
        }
    }
    else {

        // find the bounding points in xin
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            istart = 0;
            // @TODO: This should be a binary search:
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            dt = xe[j] - xin[ifirst];  // diff in x, eval to input
            out_ye[j] = yin[ifirst] + (dt * ((yin[ir] - yin[ifirst]) / (xin[ir]-xin[ifirst])) );
        }
    }
}

float VecF::sum_of_sq() {
    float *me = this->pointer();
    float total = 0;
    for( int n = 0; n < this->size(); n++) {
        total += me[n]*me[n];
    }
    return total;
}


double VecF::pearsons_r(VecF &x, VecF &y) {

    // Preparation:
    double sum_xTy = VecF::dot_product(x,y);
    double sum_x = x.sum();
    double sum_y = y.sum();
    // Could this step be sped up?
    double sum_x2 = x.sum_of_sq();
    double sum_y2 = y.sum_of_sq();
    int N = x.dim();

    // Here it is:
    // 'E' is Capital Sigma
    // r = EXY - (EXEY/N)
    //    -----------------
    //    sqrt( (EX^2 - (EX)^2/N) * (EY^2 - (EY)^2/N) )

    double top = sum_xTy - ((sum_x * sum_y)/N);
    double fbot = sum_x2 - ((sum_x*sum_x)/N);  //first part of bottom
    double sbot = sum_y2 - ((sum_y*sum_y)/N);  //second part of bottom
    return top / sqrt(fbot * sbot);

}

double VecF::covariance(VecF &x, VecF &y) {
    int i;
    int len = x.size();
    double mean_x = 0;
    double mean_y = 0;
    // get the means and x * y
    for (i = 0; i < len; ++i) {
        mean_x += x[i];
        mean_y += y[i];
    }
    mean_x /= len;
    mean_y /= len;
    double cov = 0;
    for (i = 0; i < len; ++i) {
        cov += (x[i] - mean_x) * (y[i] - mean_y);
    }
    return cov/len;
}

double VecF::euclidean(VecF &x, VecF &y) {
    VecF diff(x.size());
    double sum_of_diffs = 0;
    for (int i = 0; i < x.size(); ++i) {
        sum_of_diffs += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sqrt(sum_of_diffs);
}

float VecF::dot_product(VecF &x, VecF &y) {
    //assert(x.dim() == y.dim());
    float sum = 0;
    for (int i = 0; i < x.dim(); i++) {
        sum += (x[i] * y[i]);
    }
    return sum;
}

void VecF::chfe(VecF &xin, VecF &yin, VecF &xe, VecF &out_ye, int sorted) {
    //xin.print(); yin.print();

    if (out_ye.size() == 0) {
        float *to_take = new float[xe.size()];
        out_ye.take(xe.size(), to_take);
    }

    // Calc the derivs:
    VecF derivs;
    VecF::chim(xin,yin,derivs);
    int i,j,ir;  // i indexes xin, j indexes xnew
    int ifirst = 0;
    //derivs.print();
    // find the bounding points in xin
    int istart;


    if (sorted) {
        VecF c2(xin.size());
        VecF c3(xin.size());
        calc_cubic_coeff(xin, yin, derivs, c2, c3);
     //c2.print(); c3.print();
     //xe.print();
        istart = 0;
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            istart = i;
            chfev(xin[ifirst], yin[ifirst], derivs[ifirst], c2[ifirst], c3[ifirst], xe[j], out_ye[j]);
        }
    }
    else {

        // find the bounding points in xin
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            istart = 0;
            // @TODO: This should be a binary search:
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }

            chfev_all(xin[ifirst], xin[ir], yin[ifirst], yin[ir], derivs[ifirst], derivs[ir], xe[j], out_ye[j]);
        }
    }
}




void VecF::calc_cubic_coeff(VecF &x, VecF &y, VecF &derivs, VecF &c2, VecF &c3) {

    //  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
    float DEL1, DEL2, DELTA, H;
    for (int i = 0; i < x.size() - 1; ++i) {
        H = x[i+1] - x[i];
        DELTA = (y[i+1] - y[i])/H;
        DEL1 = (derivs[i] - DELTA)/H;
        DEL2 = (derivs[i+1] - DELTA)/H;
        c2[i] = -(DEL1+DEL1 + DEL2);
        c3[i] = (DEL1 + DEL2)/H;
    }

}

void VecF::chfev_all(float X1, float X2, float F1, float F2, float D1, float D2, float XE, float &FE) {
    float C2, C3, DEL1, DEL2, DELTA, H, X;

    H = X2 - X1;

    //  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
    DELTA = (F2 - F1)/H;
    DEL1 = (D1 - DELTA)/H;
    DEL2 = (D2 - DELTA)/H;
    C2 = -(DEL1+DEL1 + DEL2);
    C3 = (DEL1 + DEL2)/H;

    X = XE - X1;

    FE = F1 + X*(D1 + X*(C2 + X*C3));
}


void VecF::chfev(float X1, float F1, float D1, float C2, float C3, float XE, float &FE) {
    float X;
    X = XE - X1;

    FE = F1 + X*(D1 + X*(C2 + X*C3));
}


void VecF::chfe_xy(VecF &x, VecF &y, VecF &new_x, VecF &out_new_y, int sorted) {
    VecF::xy_to_x(x,y);
    chfe(x,y,new_x, out_new_y, sorted);
    x_to_xy(new_x, out_new_y);
    VecF::x_to_xy(x,y);
}

double VecF::sum_sq_res_yeqx(VecF &x, VecF &y) {
    //// PDL way
    //return sum(0.5*(($y - $x)**2));
    double __sum = 0.0;
    for (int i = 0; i < x.length(); ++i) {
        float diff = x[i] - y[i];
        __sum += 0.5*(diff*diff);
    }
    return __sum;
}

double VecF::avg_sq_res_yeqx(VecF &x, VecF &y) {
    return (sum_sq_res_yeqx(x,y))/x.length();
}

double VecF::avg_abs_diff(VecF &x, VecF &y) {
    double sum = 0.0;
    for (int i = 0; i < x.length(); ++i) {
        sum += fabs(x[i] - y[i]);
    }
    return sum/x.length();
}


void VecF::rsq_slope_intercept(VecF &x, VecF &y, double &rsq, double &slope, double &y_intercept) {
    int i;
    double mean_x = x.avg();
    double mean_y = y.avg();
    double sum_sq_res_xx = 0.0;
    double sum_sq_res_yy = 0.0;
    double sum_sq_res_xy = 0.0;

    for (i = 0; i < x.length(); ++i) {
        double x_minus_mean_i, y_minus_mean_i;
        x_minus_mean_i = ( (double)(x[i]) ) - mean_x;
        y_minus_mean_i = ( (double)(y[i]) ) - mean_y;
        sum_sq_res_xx += x_minus_mean_i*x_minus_mean_i;
        sum_sq_res_yy += y_minus_mean_i*y_minus_mean_i;
        sum_sq_res_xy += x_minus_mean_i*y_minus_mean_i;
    }
    slope = sum_sq_res_xy/sum_sq_res_xx;
    y_intercept = mean_y - (slope * mean_x);
    rsq = (sum_sq_res_xy*sum_sq_res_xy)/(sum_sq_res_xx*sum_sq_res_yy);
}


// END TEMPLATE

} // End namespace VEC
