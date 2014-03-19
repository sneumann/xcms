#ifndef DK_h
#define DK_h

#include <stdint.h>

#ifndef WIN32
//------------Unix------------
#include <sys/types.h>
#define MKDIR(path,mask) mkdir(path,mask)
#define GETCWD(buf,len)  getcwd(buf,len)
#define OPEN(buf,mode,access) open(buf,mode,access)
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
        uint32_t num_scans;
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

/*not working on windows build of bioconductor*/
//void assign_values(float64* data, uint32_t data_len, std::vector<double> & vec, int vec_len);

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

        uint32_t getTotalScanNumbers();
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
