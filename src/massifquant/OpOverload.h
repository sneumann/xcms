#ifndef OP_h
#define OP_h

#include <list>
#include <vector>

/*
  Operator Overloading
*/

/*matrix multiplication of 2 X 2 matrices*/
std::vector<double> operator * (const std::vector<double> & A,
	       const std::vector<double> & B);

/*vector or matrix addition for type double*/
std::vector<double> operator+ (const std::vector<double> & A,
	       const std::vector<double> & B);
/*vector or matrix addition for type int*/
std::vector<int> operator + (const std::vector<int> & A,
			    const std::vector<int> & B);

/*vector or matrix subtract off scalar for type double*/
std::vector<double> operator - (const std::vector<double> & A,
			       const double & b);

/*vector or matrix division by scalar for type double*/
std::vector<double> operator / (const std::vector<double> & A,
			       const double & b);

/*logical indexing for type int*/
std::vector<int> operator >= (const std::vector<int> & A,
			      const int & b);
std::vector<int> operator <= (const std::vector<int> & A,
			      const int & b);
std::vector<int> operator == (const std::vector<int> & A,
			      const int & b);

std::vector<int> operator == (const std::list<int> & A,
			      const int & b);

/*made for init data of trackers (should be the same size) */
std::list<int> operator == (const std::list<int> & A,
	const std::list<int> & B);

std::vector<int> operator != (const std::list<int> & A,
			      const int & b);

std::list<int> operator != (const std::vector<int> & A,
			      const int & b);

/*logical indexing for type double*/

std::vector<int> operator > (const std::vector<double> & A,
			      const double & b);

std::vector<int> operator >= (const std::vector<double> & A,
			      const double & b);
std::vector<int> operator <= (const std::vector<double> & A,
			      const double & b);

/*Perform this specific operation: Ab = x;*/
std::vector<double> multiplyMatVec(const std::vector<double>  & A,
				   const std::vector<double> & b);

/*element-wise multiplication of vectors*/
std::vector<double> dottimes (const std::vector<double> & A,
			      const std::vector<double> & B);

/*element-wise addition of vectors*/
std::vector<double> dotadd (const std::vector<double> & A,
			    const std::vector<double> & B);

std::vector<double> copySubIdx(const std::vector<double> & A,
			       const std::vector<int> & subidx);

std::vector<int> copySubIdx(const std::vector<int> & A,
			       const std::vector<int> & subidx);

std::vector<int> createSequence(const int start,
        const int stop, const int spacing);

double computeAnyXbar(const std::list<double> & x);

double computeAnySampVar(const std::list<double> & x);

void printvector(const std::vector<int> & myvec);

 void printvector(const std::vector<double> & myvec);

void printList(const std::list<int> & mylist);
void printList(const std::list<double> & mylist);
bool myuniqcomp(int i, int j);
int lowerBound(double val, std::vector<double> mzvals, int first, int length);
int upperBound(double val, std::vector<double> mzvals, int first, int length);

struct scanBuf {
    std::vector<double> mz;
    std::vector<double> intensity;
};


#endif
