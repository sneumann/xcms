
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "mat.h"
#include "vec.h"

#include <R.h>

//#define JTP_BOUNDS_CHECK

//#define JTP_DEBUG

namespace VEC {

// BEGIN TEMPLATE

/****************************************************************
 * MatI
 ***************************************************************/

// Constructors:
MatI::MatI() : _m(0), _n(0), _dat(0) {
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatI()!");
#endif
}

MatI::MatI(int m, int n) : _m(m), _n(n), _dat(m*n) {
#ifdef JTP_BOUNDS_CHECK
    if (m < 0 || n < 0) { Rprintf("m or n < 0"); R_ShowMessage("Serious error in obiwarp.");
}
#endif
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatI(m,n)!");
#endif
}

MatI::MatI(int m, int n, const int &val) : _m(m), _n(n), _dat(m*n, val) {
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatI(m,n,val)!");
#endif
}

MatI::MatI(int m, int n, int *arr, bool shallow) : _m(m), _n(n), _dat(m*n,arr,shallow) {
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatI(m,n,*arr,shallow) shallow=%d!\n", this->shallow());
#endif
}

MatI::MatI(const MatI &A, bool shallow) : _m(A._m), _n(A._n), _dat(A._dat, shallow) {
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatI(MatI &A,shallow) shallow=%d!\n", this->shallow());
#endif
}

void MatI::to_vec(VecI &outvec, bool shallow) {
    if (shallow) {
        outvec.set(_dat);
    }
    else {
        _dat.copy(outvec);
    }
}

void MatI::set(int m, int n, int *arr) {
    _dat.set(m*n,arr);
    _m = m;
    _n = n;
}

void MatI::set(MatI &A) {
    _dat.set(A._dat);
    _m = A._m;
    _n = A._n;
#ifdef JTP_DEBUG
    Rprintf("set called!");
#endif
}


void MatI::take(int m, int n, int *arr) {
    _dat.take(m*n,arr);
    _m = m;
    _n = n;
}

void MatI::take(MatI &A) {
    // Checking is done in Vec to ensure we're not taking a shallow!
    _dat.take(A._dat);
    _m = A._m;
    _n = A._n;
#ifdef JTP_DEBUG
    Rprintf("take called!");
#endif
}

void MatI::row_vecs(int &cnt, VecI *vecs) {
    cnt = rows();
    int _cols = cols();
    for (int i = 0; i < cnt; ++i) {
        int *ptr = this->pointer(i);
        vecs[i].set(_cols, ptr);  // shallow allocation
    }
}



bool MatI::operator==(const MatI &A) {
    // We don't care if one is shallow and the A is not!
    if (A._n == _n && A._m == _m) { // Same size
        return _dat == A._dat;
    }
    else {
        return false;
    }
}

void MatI::copy(MatI &receiver, bool shallow) const {
    receiver._m = _m;
    receiver._n = _n;
    _dat.copy(receiver._dat, shallow);
#ifdef JTP_DEBUG
    Rprintf("copy called!");
#endif
}

void MatI::file_rows_cols(std::ifstream &stream, int &rows, int &cols) {
    rows = 0;
    cols = 0;
    int BIGGEST_LINE = 1000000;
    char line[1000000];  // windows doesn't like that variable there
    stream.getline(line, BIGGEST_LINE);
    ++rows;
    char *ptr = line;
    int linelength = 0;
    while(*ptr != '\0') {
        if (*ptr == ' ') {
            *ptr = '\0';  // keep track of spaces
            ++cols;
        }
        ++ptr;
        ++linelength;
    }
    ++cols; // for the last one
    // Check for spaces on the end:
    while(1) {
        --ptr;
        // char is !\n && \r && \0
        if (*ptr == '\n' || *ptr == '\r') {
            continue;
        }
        else if (*ptr == '\0') {
            --cols;  // decrement for each space at the end
        }
        else { break; }
    }
    // Count the rows:
    while( stream.getline(line, BIGGEST_LINE) ) {
        // if the line starts with a real char:
        if (line[0] != ' ' && line[0] != '\n' && line[0] != '\r' && line[0] != '\0') {
            ++rows;
        }
    }
}

void MatI::set_from_ascii(std::ifstream &stream, int m, int n, MatI &out) {
    MatI tmp(m,n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            stream >> tmp(i,j);
        }
    }
    out.take(tmp);
}

void MatI::set_from_ascii(std::ifstream &stream, MatI &out) {
    int m,n;
    stream >> m >> n;
    MatI tmp(m,n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            stream >> tmp(i,j);
        }
    }
    out.take(tmp);
}

void MatI::set_from_ascii(const char *file, bool without_axes) {
    std::ifstream fh(file);
    if (fh.is_open()) {
        if (without_axes) {
            int m,n;
            file_rows_cols(fh,m,n);
            // Rewind the stream to beginning
            fh.clear(); // forget we saw the eof
            fh.seekg(0,std::ios::beg);
            set_from_ascii(fh,m,n,(*this));
        }
        else {
            set_from_ascii(fh,(*this));
        }
        fh.close();
    }
    else {
        Rprintf("Couldn't open %s\n", file);
        R_ShowMessage("Serious error in obiwarp.");
    }
}



MatI & MatI::operator=(const int &val) {
    _dat = val;
    return *this;
}

MatI & MatI::operator=(MatI &A) {
#ifdef JTP_DEBUG
    Rprintf("IN ASSIGNMENT OP tOP");
#endif
    if (this != &A) {
#ifdef JTP_DEBUG
        Rprintf("IN ASSIGNMENT OP MID");
#endif
        _m = A._m;
        _n = A._n;
        _dat = A._dat;
    }
    return *this;
}

MatI::~MatI( ) {
#ifdef JTP_DEBUG
    Rprintf("DESTRUCTOR");
#endif
}

/*************************
 * MATH OPERATORS
 ************************/
void MatI::operator+=(const MatI &A) {
    if (A._n == _n && A._m == _m) {
        _dat += A._dat;
    }
}

void MatI::operator-=(const MatI &A) {
    if (A._n == _n && A._m == _m) {
        _dat -= A._dat;
    }
}

void MatI::operator*=(const MatI &A) {
    if (A._n == _n && A._m == _m) {
        _dat *= A._dat;
    }
}

void MatI::operator/=(const MatI &A) {
    if (A._n == _n && A._m == _m) {
        _dat /= A._dat;
    }
}

void MatI::add(const MatI &toadd, MatI &out) {
    if (_n == toadd._n && _m == toadd._m) {
        _dat.add(toadd._dat, out._dat);
    }
}

void MatI::sub(const MatI &tosub, MatI &out) {
    if (_n == tosub._n && _m == tosub._m) {
        _dat.sub(tosub._dat, out._dat);
    }
}

void MatI::mul(const MatI &tomul, MatI &out) {
    if (_n == tomul._n && _m == tomul._m) {
        _dat.mul(tomul._dat, out._dat);
    }
}

void MatI::div(const MatI &todiv, MatI &out) {
    if (_n == todiv._n && _m == todiv._m) {
        _dat.div(todiv._dat, out._dat);
    }
}

void MatI::transpose(MatI &out) {
    MatI me(*this, 1);
    MatI tmp(me.nlen(), me.mlen());  // reverse m,n
    for (int m = 0; m < mlen(); ++m) {
        for (int n = 0; n < nlen(); ++n) {
            tmp(n,m) = me(m,n);
        }
    }
    out.take(tmp);
}


void MatI::expand(MatI &result, int match, int expand_x_lt, int expand_x_rt, int expand_y_up, int expand_y_dn, int expand_diag_lt_up, int expand_diag_rt_up, int expand_diag_lt_dn, int expand_diag_rt_dn ) {
    int i;
    int m_len = this->dim1();
    int n_len = this->dim2();
    this->copy(result);
    for (int m = 0; m < m_len; ++m) {
        for (int n = 0; n < n_len; ++n) {
            if ((*this)(m,n) == match) {
                for (i = 1; i <= expand_x_lt; ++i) {
                    if (n-i >= 0) {
                        result(m,n-i) = match;
                    }
                }
                for (i = 1; i <= expand_x_rt; ++i) {
                    if (n+i < n_len) {
                        result(m,n+i) = match;
                    }
                }
                for (i = 1; i <= expand_y_up; ++i) {
                    if (m-i >= 0) {
                        result(m-i,n) = match;
                    }
                }
                for (i = 1; i <= expand_y_dn; ++i) {
                    if (m+i < m_len) {
                        result(m+i,n) = match;
                    }
                }
                for (i = 1; i <= expand_diag_lt_up; ++i) {
                    if (n-i >= 0 && m-i >=0) {
                        result(m-i,n-i) = match;
                    }
                }
                for (i = 1; i <= expand_diag_rt_up; ++i) {
                    if (n+i < n_len && m-i >= 0) {
                        result(m-i,n+i) = match;
                    }
                }
                for (i = 1; i <= expand_diag_lt_dn; ++i) {
                    if (n-i >= 0 && m+i < m_len) {
                        result(m+i,n-i) = match;
                    }
                }
                for (i = 1; i <= expand_diag_rt_dn; ++i) {
                    if (n+i < n_len && m+i < m_len) {
                        result(m+i,n+i) = match;
                    }
                }
            }
        }
    }
}


void MatI::mask_as_vec(int return_val, MatI &mask, VecI &out) {
    _dat.mask_as_vec(return_val, mask._dat, out);
}


int MatI::sum(int m) {
    int sum = 0;
    int *ptr = pointer(m);
    for (int i = 0; i < _n; ++i) {
        sum += ptr[i];
    }
    return sum;
}


void MatI::print(bool without_axes) {
    MatI tmp((*this),1);
    if (!without_axes) {
        //std::cout << _m << ' ' << _n << std::endl;
    }
    for (int m = 0; m < _m; ++m) {
        int n;
        for (n = 0; n < _n - 1; ++n) {
            //std::cout << tmp(m,n) << " ";
        }
        //std::cout << tmp(m,n);
        //std::cout << std::endl;
    }
}

void MatI::print(const char *filename, bool without_axes) {
    std::ofstream fh(filename);
    if (!fh) {
        //std::cout << "Error opening file " << filename << std::endl;
    }
    this->print(fh, without_axes);
    fh.close();
}

void MatI::print(std::ostream &fout, bool without_axes) {
    int m;
    if (!without_axes) {
        fout << _m << ' ' << _n << std::endl;
    }
    for (m = 0; m < _m; m++) {
        int n;
        for (n = 0; n < _n - 1; n++) {
            fout << _dat[(m*_n)+n] << " ";
        }
        fout << _dat[m*_n+n];
        fout << std::endl;
    }
}

void MatI::write(const char *file) {
    if (file != NULL) {
        FILE *fh = fopen(file, "wb");
        fwrite(&_m, sizeof(int), 1, fh);
        fwrite(&_n, sizeof(int), 1, fh);
        fwrite((int*)(_dat), sizeof(int), _m*_n, fh);
        fclose(fh);
    }
    // else {
    //     fwrite(&_m, sizeof(int), 1, stdout);
    //     fwrite(&_n, sizeof(int), 1, stdout);
    //     fwrite((int*)(_dat), sizeof(int), _m*_n, stdout);
    // }
}



/****************************************************************
 * MatD
 ***************************************************************/

// Constructors:
MatD::MatD() : _m(0), _n(0), _dat(0) {
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatD()!");
#endif
}

MatD::MatD(int m, int n) : _m(m), _n(n), _dat(m*n) {
#ifdef JTP_BOUNDS_CHECK
    if (m < 0 || n < 0) { Rprintf("m or n < 0"); R_ShowMessage("Serious error in obiwarp.");
#endif
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatD(m,n)!");
#endif
}

MatD::MatD(int m, int n, const double &val) : _m(m), _n(n), _dat(m*n, val) {
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatD(m,n,val)!");
#endif
}

MatD::MatD(int m, int n, double *arr, bool shallow) : _m(m), _n(n), _dat(m*n,arr,shallow) {
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatD(m,n,*arr,shallow) shallow=%d!\n", this->shallow());
#endif
}

MatD::MatD(const MatD &A, bool shallow) : _m(A._m), _n(A._n), _dat(A._dat, shallow) {
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatD(MatD &A,shallow) shallow=%d!\n", this->shallow());
#endif
}

void MatD::to_vec(VecD &outvec, bool shallow) {
    if (shallow) {
        outvec.set(_dat);
    }
    else {
        _dat.copy(outvec);
    }
}

void MatD::set(int m, int n, double *arr) {
    _dat.set(m*n,arr);
    _m = m;
    _n = n;
}

void MatD::set(MatD &A) {
    _dat.set(A._dat);
    _m = A._m;
    _n = A._n;
#ifdef JTP_DEBUG
    Rprintf("set called!");
#endif
}


void MatD::take(int m, int n, double *arr) {
    _dat.take(m*n,arr);
    _m = m;
    _n = n;
}

void MatD::take(MatD &A) {
    // Checking is done in Vec to ensure we're not taking a shallow!
    _dat.take(A._dat);
    _m = A._m;
    _n = A._n;
#ifdef JTP_DEBUG
    Rprintf("take called!");
#endif
}

void MatD::row_vecs(int &cnt, VecD *vecs) {
    cnt = rows();
    int _cols = cols();
    for (int i = 0; i < cnt; ++i) {
        double *ptr = this->pointer(i);
        vecs[i].set(_cols, ptr);  // shallow allocation
    }
}



bool MatD::operator==(const MatD &A) {
    // We don't care if one is shallow and the A is not!
    if (A._n == _n && A._m == _m) { // Same size
        return _dat == A._dat;
    }
    else {
        return false;
    }
}

void MatD::copy(MatD &receiver, bool shallow) const {
    receiver._m = _m;
    receiver._n = _n;
    _dat.copy(receiver._dat, shallow);
#ifdef JTP_DEBUG
    Rprintf("copy called!");
#endif
}

void MatD::file_rows_cols(std::ifstream &stream, int &rows, int &cols) {
    rows = 0;
    cols = 0;
    int BIGGEST_LINE = 1000000;
    char line[1000000];  // windows doesn't like that variable there
    stream.getline(line, BIGGEST_LINE);
    ++rows;
    char *ptr = line;
    int linelength = 0;
    while(*ptr != '\0') {
        if (*ptr == ' ') {
            *ptr = '\0';  // keep track of spaces
            ++cols;
        }
        ++ptr;
        ++linelength;
    }
    ++cols; // for the last one
    // Check for spaces on the end:
    while(1) {
        --ptr;
        // char is !\n && \r && \0
        if (*ptr == '\n' || *ptr == '\r') {
            continue;
        }
        else if (*ptr == '\0') {
            --cols;  // decrement for each space at the end
        }
        else { break; }
    }
    // Count the rows:
    while( stream.getline(line, BIGGEST_LINE) ) {
        // if the line starts with a real char:
        if (line[0] != ' ' && line[0] != '\n' && line[0] != '\r' && line[0] != '\0') {
            ++rows;
        }
    }
}

void MatD::set_from_ascii(std::ifstream &stream, int m, int n, MatD &out) {
    MatD tmp(m,n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            stream >> tmp(i,j);
        }
    }
    out.take(tmp);
}

void MatD::set_from_ascii(std::ifstream &stream, MatD &out) {
    int m,n;
    stream >> m >> n;
    MatD tmp(m,n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            stream >> tmp(i,j);
        }
    }
    out.take(tmp);
}

void MatD::set_from_ascii(const char *file, bool without_axes) {
    std::ifstream fh(file);
    if (fh.is_open()) {
        if (without_axes) {
            int m,n;
            file_rows_cols(fh,m,n);
            // Rewind the stream to beginning
            fh.clear(); // forget we saw the eof
            fh.seekg(0,std::ios::beg);
            set_from_ascii(fh,m,n,(*this));
        }
        else {
            set_from_ascii(fh,(*this));
        }
        fh.close();
    }
    else {
        Rprintf("Couldn't open %s\n", file);
        R_ShowMessage("Serious error in obiwarp.");
    }
}



MatD & MatD::operator=(const double &val) {
    _dat = val;
    return *this;
}

MatD & MatD::operator=(MatD &A) {
#ifdef JTP_DEBUG
    Rprintf("IN ASSIGNMENT OP tOP");
#endif
    if (this != &A) {
#ifdef JTP_DEBUG
        Rprintf("IN ASSIGNMENT OP MID");
#endif
        _m = A._m;
        _n = A._n;
        _dat = A._dat;
    }
    return *this;
}

MatD::~MatD( ) {
#ifdef JTP_DEBUG
    Rprintf("DESTRUCTOR");
#endif
}

/*************************
 * MATH OPERATORS
 ************************/
void MatD::operator+=(const MatD &A) {
    if (A._n == _n && A._m == _m) {
        _dat += A._dat;
    }
}

void MatD::operator-=(const MatD &A) {
    if (A._n == _n && A._m == _m) {
        _dat -= A._dat;
    }
}

void MatD::operator*=(const MatD &A) {
    if (A._n == _n && A._m == _m) {
        _dat *= A._dat;
    }
}

void MatD::operator/=(const MatD &A) {
    if (A._n == _n && A._m == _m) {
        _dat /= A._dat;
    }
}

void MatD::add(const MatD &toadd, MatD &out) {
    if (_n == toadd._n && _m == toadd._m) {
        _dat.add(toadd._dat, out._dat);
    }
}

void MatD::sub(const MatD &tosub, MatD &out) {
    if (_n == tosub._n && _m == tosub._m) {
        _dat.sub(tosub._dat, out._dat);
    }
}

void MatD::mul(const MatD &tomul, MatD &out) {
    if (_n == tomul._n && _m == tomul._m) {
        _dat.mul(tomul._dat, out._dat);
    }
}

void MatD::div(const MatD &todiv, MatD &out) {
    if (_n == todiv._n && _m == todiv._m) {
        _dat.div(todiv._dat, out._dat);
    }
}

void MatD::transpose(MatD &out) {
    MatD me(*this, 1);
    MatD tmp(me.nlen(), me.mlen());  // reverse m,n
    for (int m = 0; m < mlen(); ++m) {
        for (int n = 0; n < nlen(); ++n) {
            tmp(n,m) = me(m,n);
        }
    }
    out.take(tmp);
}


void MatD::expand(MatD &result, double match, int expand_x_lt, int expand_x_rt, int expand_y_up, int expand_y_dn, int expand_diag_lt_up, int expand_diag_rt_up, int expand_diag_lt_dn, int expand_diag_rt_dn ) {
    int i;
    int m_len = this->dim1();
    int n_len = this->dim2();
    this->copy(result);
    for (int m = 0; m < m_len; ++m) {
        for (int n = 0; n < n_len; ++n) {
            if ((*this)(m,n) == match) {
                for (i = 1; i <= expand_x_lt; ++i) {
                    if (n-i >= 0) {
                        result(m,n-i) = match;
                    }
                }
                for (i = 1; i <= expand_x_rt; ++i) {
                    if (n+i < n_len) {
                        result(m,n+i) = match;
                    }
                }
                for (i = 1; i <= expand_y_up; ++i) {
                    if (m-i >= 0) {
                        result(m-i,n) = match;
                    }
                }
                for (i = 1; i <= expand_y_dn; ++i) {
                    if (m+i < m_len) {
                        result(m+i,n) = match;
                    }
                }
                for (i = 1; i <= expand_diag_lt_up; ++i) {
                    if (n-i >= 0 && m-i >=0) {
                        result(m-i,n-i) = match;
                    }
                }
                for (i = 1; i <= expand_diag_rt_up; ++i) {
                    if (n+i < n_len && m-i >= 0) {
                        result(m-i,n+i) = match;
                    }
                }
                for (i = 1; i <= expand_diag_lt_dn; ++i) {
                    if (n-i >= 0 && m+i < m_len) {
                        result(m+i,n-i) = match;
                    }
                }
                for (i = 1; i <= expand_diag_rt_dn; ++i) {
                    if (n+i < n_len && m+i < m_len) {
                        result(m+i,n+i) = match;
                    }
                }
            }
        }
    }
}


void MatD::mask_as_vec(double return_val, MatI &mask, VecD &out) {
    _dat.mask_as_vec(return_val, mask._dat, out);
}


double MatD::sum(int m) {
    double sum = 0;
    double *ptr = pointer(m);
    for (int i = 0; i < _n; ++i) {
        sum += ptr[i];
    }
    return sum;
}


void MatD::print(bool without_axes) {
    MatD tmp((*this),1);
    if (!without_axes) {
        //std::cout << _m << ' ' << _n << std::endl;
    }
    for (int m = 0; m < _m; ++m) {
        int n;
        for (n = 0; n < _n - 1; ++n) {
            //std::cout << tmp(m,n) << " ";
        }
        //std::cout << tmp(m,n);
        //std::cout << std::endl;
    }
}

void MatD::print(const char *filename, bool without_axes) {
    std::ofstream fh(filename);
    if (!fh) {
        //std::cout << "Error opening file " << filename << std::endl;
    }
    this->print(fh, without_axes);
    fh.close();
}

void MatD::print(std::ostream &fout, bool without_axes) {
    int m;
    if (!without_axes) {
        fout << _m << ' ' << _n << std::endl;
    }
    for (m = 0; m < _m; m++) {
        int n;
        for (n = 0; n < _n - 1; n++) {
            fout << _dat[(m*_n)+n] << " ";
        }
        fout << _dat[m*_n+n];
        fout << std::endl;
    }
}

void MatD::write(const char *file) {
    if (file != NULL) {
        FILE *fh = fopen(file, "wb");
        fwrite(&_m, sizeof(int), 1, fh);
        fwrite(&_n, sizeof(int), 1, fh);
        fwrite((double*)(_dat), sizeof(double), _m*_n, fh);
        fclose(fh);
    }
    // else {
    //     fwrite(&_m, sizeof(int), 1, stdout);
    //     fwrite(&_n, sizeof(int), 1, stdout);
    //     fwrite((double*)(_dat), sizeof(double), _m*_n, stdout);
    // }
}



/****************************************************************
 * MatF
 ***************************************************************/

// Constructors:
MatF::MatF() : _m(0), _n(0), _dat(0) {
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatF()!");
#endif
}

MatF::MatF(int m, int n) : _m(m), _n(n), _dat(m*n) {
#ifdef JTP_BOUNDS_CHECK
    if (m < 0 || n < 0) { Rprintf("m or n < 0"); R_ShowMessage("Serious error in obiwarp.");
#endif
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatF(m,n)!");
#endif
}

MatF::MatF(int m, int n, const float &val) : _m(m), _n(n), _dat(m*n, val) {
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatF(m,n,val)!");
#endif
}

MatF::MatF(int m, int n, float *arr, bool shallow) : _m(m), _n(n), _dat(m*n,arr,shallow) {
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatF(m,n,*arr,shallow) shallow=%d!\n", this->shallow());
#endif
}

MatF::MatF(const MatF &A, bool shallow) : _m(A._m), _n(A._n), _dat(A._dat, shallow) {
#ifdef JTP_DEBUG
    Rprintf("CONSTRUCTOR MatF(MatF &A,shallow) shallow=%d!\n", this->shallow());
#endif
}

void MatF::to_vec(VecF &outvec, bool shallow) {
    if (shallow) {
        outvec.set(_dat);
    }
    else {
        _dat.copy(outvec);
    }
}

void MatF::set(int m, int n, float *arr) {
    _dat.set(m*n,arr);
    _m = m;
    _n = n;
}

void MatF::set(MatF &A) {
    _dat.set(A._dat);
    _m = A._m;
    _n = A._n;
#ifdef JTP_DEBUG
    Rprintf("set called!");
#endif
}


void MatF::take(int m, int n, float *arr) {
    _dat.take(m*n,arr);
    _m = m;
    _n = n;
}

void MatF::take(MatF &A) {
    // Checking is done in Vec to ensure we're not taking a shallow!
    _dat.take(A._dat);
    _m = A._m;
    _n = A._n;
#ifdef JTP_DEBUG
    Rprintf("take called!");
#endif
}

void MatF::row_vecs(int &cnt, VecF *vecs) {
    cnt = rows();
    int _cols = cols();
    for (int i = 0; i < cnt; ++i) {
        float *ptr = this->pointer(i);
        vecs[i].set(_cols, ptr);  // shallow allocation
    }
}



bool MatF::operator==(const MatF &A) {
    // We don't care if one is shallow and the A is not!
    if (A._n == _n && A._m == _m) { // Same size
        return _dat == A._dat;
    }
    else {
        return false;
    }
}

void MatF::copy(MatF &receiver, bool shallow) const {
    receiver._m = _m;
    receiver._n = _n;
    _dat.copy(receiver._dat, shallow);
#ifdef JTP_DEBUG
    Rprintf("copy called!");
#endif
}

void MatF::file_rows_cols(std::ifstream &stream, int &rows, int &cols) {
    rows = 0;
    cols = 0;
    int BIGGEST_LINE = 1000000;
    char line[1000000];  // windows doesn't like that variable there
    stream.getline(line, BIGGEST_LINE);
    ++rows;
    char *ptr = line;
    int linelength = 0;
    while(*ptr != '\0') {
        if (*ptr == ' ') {
            *ptr = '\0';  // keep track of spaces
            ++cols;
        }
        ++ptr;
        ++linelength;
    }
    ++cols; // for the last one
    // Check for spaces on the end:
    while(1) {
        --ptr;
        // char is !\n && \r && \0
        if (*ptr == '\n' || *ptr == '\r') {
            continue;
        }
        else if (*ptr == '\0') {
            --cols;  // decrement for each space at the end
        }
        else { break; }
    }
    // Count the rows:
    while( stream.getline(line, BIGGEST_LINE) ) {
        // if the line starts with a real char:
        if (line[0] != ' ' && line[0] != '\n' && line[0] != '\r' && line[0] != '\0') {
            ++rows;
        }
    }
}

void MatF::set_from_ascii(std::ifstream &stream, int m, int n, MatF &out) {
    MatF tmp(m,n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            stream >> tmp(i,j);
        }
    }
    out.take(tmp);
}

void MatF::set_from_ascii(std::ifstream &stream, MatF &out) {
    int m,n;
    stream >> m >> n;
    MatF tmp(m,n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            stream >> tmp(i,j);
        }
    }
    out.take(tmp);
}

void MatF::set_from_ascii(const char *file, bool without_axes) {
    std::ifstream fh(file);
    if (fh.is_open()) {
        if (without_axes) {
            int m,n;
            file_rows_cols(fh,m,n);
            // Rewind the stream to beginning
            fh.clear(); // forget we saw the eof
            fh.seekg(0,std::ios::beg);
            set_from_ascii(fh,m,n,(*this));
        }
        else {
            set_from_ascii(fh,(*this));
        }
        fh.close();
    }
    else {
        Rprintf("Couldn't open %s\n", file);
        R_ShowMessage("Serious error in obiwarp.");
    }
}



MatF & MatF::operator=(const float &val) {
    _dat = val;
    return *this;
}

MatF & MatF::operator=(MatF &A) {
#ifdef JTP_DEBUG
    Rprintf("IN ASSIGNMENT OP tOP");
#endif
    if (this != &A) {
#ifdef JTP_DEBUG
        Rprintf("IN ASSIGNMENT OP MID");
#endif
        _m = A._m;
        _n = A._n;
        _dat = A._dat;
    }
    return *this;
}

MatF::~MatF( ) {
#ifdef JTP_DEBUG
    Rprintf("DESTRUCTOR");
#endif
}

/*************************
 * MATH OPERATORS
 ************************/
void MatF::operator+=(const MatF &A) {
    if (A._n == _n && A._m == _m) {
        _dat += A._dat;
    }
}

void MatF::operator-=(const MatF &A) {
    if (A._n == _n && A._m == _m) {
        _dat -= A._dat;
    }
}

void MatF::operator*=(const MatF &A) {
    if (A._n == _n && A._m == _m) {
        _dat *= A._dat;
    }
}

void MatF::operator/=(const MatF &A) {
    if (A._n == _n && A._m == _m) {
        _dat /= A._dat;
    }
}

void MatF::add(const MatF &toadd, MatF &out) {
    if (_n == toadd._n && _m == toadd._m) {
        _dat.add(toadd._dat, out._dat);
    }
}

void MatF::sub(const MatF &tosub, MatF &out) {
    if (_n == tosub._n && _m == tosub._m) {
        _dat.sub(tosub._dat, out._dat);
    }
}

void MatF::mul(const MatF &tomul, MatF &out) {
    if (_n == tomul._n && _m == tomul._m) {
        _dat.mul(tomul._dat, out._dat);
    }
}

void MatF::div(const MatF &todiv, MatF &out) {
    if (_n == todiv._n && _m == todiv._m) {
        _dat.div(todiv._dat, out._dat);
    }
}

void MatF::transpose(MatF &out) {
    MatF me(*this, 1);
    MatF tmp(me.nlen(), me.mlen());  // reverse m,n
    for (int m = 0; m < mlen(); ++m) {
        for (int n = 0; n < nlen(); ++n) {
            tmp(n,m) = me(m,n);
        }
    }
    out.take(tmp);
}


void MatF::expand(MatF &result, float match, int expand_x_lt, int expand_x_rt, int expand_y_up, int expand_y_dn, int expand_diag_lt_up, int expand_diag_rt_up, int expand_diag_lt_dn, int expand_diag_rt_dn ) {
    int i;
    int m_len = this->dim1();
    int n_len = this->dim2();
    this->copy(result);
    for (int m = 0; m < m_len; ++m) {
        for (int n = 0; n < n_len; ++n) {
            if ((*this)(m,n) == match) {
                for (i = 1; i <= expand_x_lt; ++i) {
                    if (n-i >= 0) {
                        result(m,n-i) = match;
                    }
                }
                for (i = 1; i <= expand_x_rt; ++i) {
                    if (n+i < n_len) {
                        result(m,n+i) = match;
                    }
                }
                for (i = 1; i <= expand_y_up; ++i) {
                    if (m-i >= 0) {
                        result(m-i,n) = match;
                    }
                }
                for (i = 1; i <= expand_y_dn; ++i) {
                    if (m+i < m_len) {
                        result(m+i,n) = match;
                    }
                }
                for (i = 1; i <= expand_diag_lt_up; ++i) {
                    if (n-i >= 0 && m-i >=0) {
                        result(m-i,n-i) = match;
                    }
                }
                for (i = 1; i <= expand_diag_rt_up; ++i) {
                    if (n+i < n_len && m-i >= 0) {
                        result(m-i,n+i) = match;
                    }
                }
                for (i = 1; i <= expand_diag_lt_dn; ++i) {
                    if (n-i >= 0 && m+i < m_len) {
                        result(m+i,n-i) = match;
                    }
                }
                for (i = 1; i <= expand_diag_rt_dn; ++i) {
                    if (n+i < n_len && m+i < m_len) {
                        result(m+i,n+i) = match;
                    }
                }
            }
        }
    }
}


void MatF::mask_as_vec(float return_val, MatI &mask, VecF &out) {
    _dat.mask_as_vec(return_val, mask._dat, out);
}


float MatF::sum(int m) {
    float sum = 0;
    float *ptr = pointer(m);
    for (int i = 0; i < _n; ++i) {
        sum += ptr[i];
    }
    return sum;
}


void MatF::print(bool without_axes) {
    MatF tmp((*this),1);
    if (!without_axes) {
        //std::cout << _m << ' ' << _n << std::endl;
    }
    for (int m = 0; m < _m; ++m) {
        int n;
        for (n = 0; n < _n - 1; ++n) {
            //std::cout << tmp(m,n) << " ";
        }
        //std::cout << tmp(m,n);
        //std::cout << std::endl;
    }
}

void MatF::print(int __m, int __n,bool without_axes) {
    MatF tmp((*this),1);
    if (!without_axes) {
        //std::cout << _m << ' ' << _n << std::endl<< std::endl;
    }
    for (int m = 0; m < _m; ++m) {
      //std::cout << tmp(m,__n) << " ";
        }
}

void MatF::print(const char *filename, bool without_axes) {
    std::ofstream fh(filename);
    if (!fh) {
        //std::cout << "Error opening file " << filename << std::endl;
    }
    this->print(fh, without_axes);
    fh.close();
}

void MatF::print(std::ostream &fout, bool without_axes) {
    int m;
    if (!without_axes) {
        fout << _m << ' ' << _n << std::endl;
    }
    for (m = 0; m < _m; m++) {
        int n;
        for (n = 0; n < _n - 1; n++) {
            fout << _dat[(m*_n)+n] << " ";
        }
        fout << _dat[m*_n+n];
        fout << std::endl;
    }
}

void MatF::write(const char *file) {
    if (file != NULL) {
        FILE *fh = fopen(file, "wb");
        fwrite(&_m, sizeof(int), 1, fh);
        fwrite(&_n, sizeof(int), 1, fh);
        fwrite((float*)(_dat), sizeof(float), _m*_n, fh);
        fclose(fh);
    }
    // else {
    //     fwrite(&_m, sizeof(int), 1, stdout);
    //     fwrite(&_n, sizeof(int), 1, stdout);
    //     fwrite((float*)(_dat), sizeof(float), _m*_n, stdout);
    // }
}

// END TEMPLATE

} // End namespace VEC
