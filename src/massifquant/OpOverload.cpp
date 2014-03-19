//OpOverload.cpp

#include <iostream>
#include <stdio.h>
#include <math.h>
#include "OpOverload.h"

#include <R.h>

const int OPDIMM = 4;
const int OPDIMV = 2;

using namespace std;

std::vector<double> operator * (const std::vector<double> & A,
	       const std::vector<double> & B) {
  //C = AB;
  std::vector<double> C(OPDIMM, 0);
  C[0] = A[0]*B[0] + A[1]*B[2];
  C[1] = A[0]*B[1] + A[1]*B[3];
  C[2] = A[2]*B[0] + A[3]*B[2];
  C[3] = A[2]*B[1] + A[3]*B[3];

  return C;
}

std::vector<double> operator + (const std::vector<double> & A,
	       const std::vector<double> & B) {
  int m = A.size();
  std::vector<double> C(m,0);
  for (int j = 0; j < m; j++) {
    C[j] = A[j] + B[j];
  }
  return C;
}

std::vector<int> operator + (const std::vector<int> & A,
			    const std::vector<int> & B) {
  int m = A.size();
  std::vector<int> C(m, 0);
  for (int j = 0; j < m; j++) {
        C[j] = A[j] + B[j];
  }
  return C;
}

std::vector<double> operator - (const std::vector<double> & A,
	       const double & b) {
  int m = A.size();
  std::vector<double> C(m,0);
  for (int j = 0; j < m; j++) {
    C[j] = A[j] -  b;
  }
  return C;

}

std::vector<double> operator / (const std::vector<double> & A,
			       const double & b) {
  int m = A.size();
  std::vector<double>  C(m,0);
  for (int j = 0; j < m; j++) {
    C[j] = A[j] / b;
  }
  return C;

}

//integer logicals
std::vector<int> operator >= (const std::vector<int> & A,
			      const int & b) {
  int m = A.size();
  std::vector<int> C(m,0);
  for (int i = 0; i < m; i++) {
    if (A.at(i) >= b)
      C.at(i) = 1;
  }
  return C;
}

std::vector<int> operator <= (const std::vector<int> & A,
			      const int & b) {
  int m = A.size();
  std::vector<int> C(m,0);
  for (int i = 0; i < m; i++) {
    if (A.at(i) <= b)
      C.at(i) = 1;
  }
  return C;
}

std::vector<int> operator == (const std::vector<int> & A,
			      const int & b) {
  int m = A.size();
  std::vector<int> C;
  for (int i = 0; i < m; i++) {
    if (A.at(i) == b)
      C.push_back(i);//returns the index
  }
  return C;
}

std::vector<int> operator == (const std::list<int> & A,
			      const int & b) {
  std::vector<int> C;
  std::list<int>::const_iterator it;
  int i = 0;
  for (it = A.begin(); it != A.end(); ++it) {
    if (*it == b) {
      C.push_back(i);//returns the index
    }
    i++;
  }
   return C;
}

std::list<int> operator == (const std::list<int> & A,
	const std::list<int> & B) {

    if (A.size() != B.size()) Rf_error("assertion failled in massifquant\n");
    std::list<int> C;
    std::list<int>::const_iterator itA;
    std::list<int>::const_iterator itB = B.begin();
    int idx = 0;
    for (itA = A.begin(); itA != A.end(); ++itA) {
        if (*itA == *itB) { C.push_back(idx); }
        itB++;
        idx++;
    }
    return C;
}

std::vector<int> operator != (const std::list<int> & A,
			      const int & b) {
  std::vector<int> C;
  std::list<int>::const_iterator it;
  int i = 0;
  for (it = A.begin(); it != A.end(); ++it) {
    if (*it != b) {
      C.push_back(i);//returns the index
    }
    i++;
  }
   return C;
}

/*use for miss list */
std::list<int> operator != (const std::vector<int> & A,
            			      const int & b) {
  std::list<int> C;
  std::vector<int>::const_iterator it;
  int i = 0;
  for (it = A.begin(); it != A.end(); ++it) {
    if (*it != b) {
      C.push_back(i);//returns the index
    }
    i++;
  }
   return C;
}

//double logicals
std::vector<int> operator > (const std::vector<double> & A,
			      const double & b) {
  int m = A.size();
  std::vector<int> C(m,0);
  for (int i = 0; i < m; i++) {
    if (A.at(i) > b)
      C.at(i) = 1;
  }
  return C;
}

std::vector<int> operator >= (const std::vector<double> & A,
			      const double & b) {
  int m = A.size();
  std::vector<int> C(m,0);
  for (int i = 0; i < m; i++) {
    if (A.at(i) >= b)
      C.at(i) = 1;
  }
  return C;
}

std::vector<int> operator <= (const std::vector<double> & A,
			      const double & b) {
  int m = A.size();
  std::vector<int> C(m,0);
  for (int i = 0; i < m; i++) {
    if (A.at(i) <= b)
      C.at(i) = 1;
  }
  return C;
}

std::vector<double> multiplyMatVec(const std::vector<double>  & A,
		   const std::vector<double> & b) {
  //Ab = x;
  std::vector<double> x(OPDIMV,0);
  x[0] = A[0]*b[0] + A[1]*b[1];
  x[1] = A[2]*b[0] + A[3]*b[1];

  return x;
}

std::vector<double> dottimes (const std::vector<double> & A,
				const std::vector<double> & B) {
  int m = A.size();
  std::vector<double> C(m, 0);
  for (int j = 0; j < m; j++) {
    C[j] = A[j]*B[j];
  }
  return C;
}

std::vector<double> dotadd (const std::vector<double> & A,
				const std::vector<double> & B) {
  int m = A.size();
  std::vector<double> C(m, 0);
  for (int j = 0; j < m; j++) {
    C[j] = A[j] + B[j];
  }
  return C;
}

std::vector<double> copySubIdx(const std::vector<double> & A,
			       const std::vector<int> & subidx) {
  int m = subidx.size();
  std::vector<double> Asubset(m,0);
  std::vector<int>::const_iterator it;
  int i = 0;
  for (it = subidx.begin(); it != subidx.end(); ++it) {
        Asubset[i] = A.at(*it); //assign it the current value;
        i++;
  }
  return Asubset;
}

std::vector<int> copySubIdx(const std::vector<int> & A,
        const std::vector<int> & subidx) {
    int m = subidx.size();
    std::vector<int> Asubset(m,0);
    std::vector<int>::const_iterator it;
    int i = 0;
    for (it = subidx.begin(); it != subidx.end(); ++it) {
        Asubset[i] = A.at(*it); //assign it the current value;
        i++;
    }
    return Asubset;

}

std::vector<int> createSequence(const int start,
        const int stop, const int spacing) {

    int upperBound = stop + 1;
    std::vector<int> seq(stop - start + 1, 0);
    int count = 0;
    for(int i = start; i < upperBound; i = i + spacing) {
        seq[count] = i;
        count++;
    }
    return seq;
}


double computeAnyXbar(const std::list<double> & x) {

    double xbar = 0;
    std::list<double>::const_iterator it;
    for (it = x.begin(); it != x.end(); ++it) {
        xbar +=  *it;
    }
    xbar = xbar/x.size();
    return xbar;
}

double computeAnySampVar(const std::list<double> & x) {

    double s2 = 0;
    double xbar = computeAnyXbar(x);
    std::list<double>::const_iterator it;
    for (it = x.begin(); it != x.end(); ++it) {
        s2 += pow(*it - xbar, 2);
    }
    s2 = s2/(x.size() - 1); // n - 1 is unbiased
    return s2;
}

void printvector(std::vector<int> & myvec) {
      for (size_t i = 0; i < myvec.size(); i++) {
        Rprintf("%d ", myvec.at(i));
                }
        Rprintf("\n");
}

void printList(const std::list<int> & mylist) {
        std::list<int>::const_iterator it;
            for (it = mylist.begin(); it != mylist.end(); ++it) {
                        Rprintf("%d  ", *it);
                            }
                Rprintf("\n");
}

void printList(const std::list<double> & mylist) {
        std::list<double>::const_iterator it;
            for (it = mylist.begin(); it != mylist.end(); ++it) {
                        Rprintf("%f ", *it);
                            }
                Rprintf("\n");
}


void printvector(const std::vector<double> & myvec) {
      for (size_t i = 0; i < myvec.size(); i++) {
              Rprintf("%f", myvec.at(i));
              Rprintf(" \n");
                  }
        Rprintf("\n");
}

int lowerBound(double val, std::vector<double> mzvals, int first, int length){
    int half, mid;
    while (length > 0) {
        half = length >> 1;
        mid = first;
        mid += half;
        if ( mzvals.at(mid) < val){
            first = mid;
            first ++;
            length = length - half -1;
        }
        else length = half;
    }
    return(first);
}

int upperBound(double val, std::vector<double> mzvals, int first, int length){
    int half, mid;
    while (length > 0) {
        half = length >> 1;
        mid = first;
        mid += half;
        if (val < mzvals.at(mid) ){
         length = half;
        }
        else {
            first = mid;
            first ++;
            length = length - half -1;
        }
    }
    return(first);
}
