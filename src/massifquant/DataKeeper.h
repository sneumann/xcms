#ifndef DK_h
#define DK_h

#ifndef WIN32
//------------Unix------------
#include <sys/types.h>
#define MKDIR(path,mask) mkdir(path,mask)
#define GETCWD(buf,len)  getcwd(buf,len)
#define OPEN(buf,mode,access) open(buf,mode,access)
//------------Unix: 32-bit/64-bit compatibility------------
typedef __int64_t     int64;
typedef __uint64_t    uint64;
typedef __int32_t     int32;
typedef __uint32_t    uint32;
typedef __int16_t     int16;
typedef __uint16_t    uint16;
typedef char          int8;
typedef unsigned char uint8;
typedef double        float64;
#else

//------------Win: 32-bit/64-bit compatibility------------
typedef __int64            int64;
typedef unsigned __int64   uint64;
typedef __int32            int32;
typedef unsigned int       uint32;
typedef __int16            int16;
typedef unsigned __int16   uint16;
typedef __int8             int8;
typedef unsigned __int8    uint8;
//typedef __float64          float64;
typedef double 			   float64;
#endif

//outside libraries
#include <string>
#include <vector>

/* Not able to find these header files */
#include <R.h>
#include <Rdefines.h>


const int MZITR = 0;
const int IITR = 1;
const int FILECHARMAX = 300;
const int RANGEMAXNUM = 6;
class DataKeeper {

    private:
        uint32 num_scans;
        std::vector<int> scan_idx;
        std::vector<double> rt; //a single big array
        std::vector<double> mz; //an
        std::vector<double> intensity;
   
        struct scanBuf * scbuf;
        double * pmz;
        double * pinten; 
        int * pscanindex;
        int nmz;
        int lastScan;
        double * pscantime;
        char filename[FILECHARMAX];

        double initMZS2;
        double initIS2; //var
        double initIS;   //sd

        void printVec(const std::vector<double> & myvec);
        void printList(const std::list<int> & mylist);
void printList(const std::list<double> & mylist);
        void assign_values(float64* data, uint32 data_len, std::vector<double> & vec, int vec_len); 
        
        std::vector<double> privGetMZScan(int s);
        std::vector<double> privGetIScan(int s);
        
        void privGetScanXcms(int scan, std::vector<double> & mzScan,
                std::vector<double> & intenScan);

    //want to be able to change data
        std::vector<double> transformIntensity(std::vector<double> & A); 
        void transformIntensityR(); 

    public:

        DataKeeper(SEXP mz, SEXP inten, SEXP scanindex, SEXP ls, SEXP scantime);
        DataKeeper(const char* dotplms1);
        ~DataKeeper();
        
        uint32 getTotalScanNumbers();
        int getTotalCentroidCount();
        double getInitMZS2();
        double getInitIS2();
        double getInitIS();

        std::vector<double> getMZScan(int s);
        std::vector<double> getIScan(int s);

        void  getScanMQ(int s, std::vector<double> & mzScan, std::vector<double> & intenScan);
        void  getScanXcms(int scan, int nmz, int lastScan, std::vector<double> & mzScan, std::vector<double> & intenScan);

        double getScanTime(int s);

        void ghostScan();
        
        void ghostScanR();  

        };
#endif
