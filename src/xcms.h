void ProfBinLin(double *xvals, double *yvals, int *numin,
                double *xstart, double *xend, int *numout, double *out);

void ProfBinLinM(double *xvals, double *yvals, int *numin, int *mindex, int *nummi,
                 double *xstart, double *xend, int *numout, double *out);

void ProfBinLinBase(double *xvals, double *yvals, int *numin, double *baselevel, double *basespace,
                    double *xstart, double *xend, int *numout, double *out);

void ProfBinLinBaseM(double *xvals, double *yvals, int *numin, int *mindex, int *nummi,
                     double *baselevel, double *basespace, double *xstart, double *xend,
                     int *numout, double *out);

void ProfIntLin(double *xvals, double *yvals, int *numin,
                double *xstart, double *xend, int *numout, double *out);

void ProfIntLinM(double *xvals, double *yvals, int *numin, int *mindex, int *nummi,
                 double *xstart, double *xend, int *numout, double *out);

void ProfBin(double *xvals, double *yvals, int *numin,
             double *xstart, double *xend, int *numout, double *out);

void ProfBinM(double *xvals, double *yvals, int *numin, int *mindex, int *nummi,
              double *xstart, double *xend, int *numout, double *out);

void ProfMaxIdx(double *xvals, double *yvals, int *numin,
                double *xstart, double *xend, int *numout, int *out);

void ProfMaxIdxM(double *xvals, double *yvals, int *numin, int *mindex, int *nummi,
                 double *xstart, double *xend, int *numout, int *out);

void MedianFilter(double *inmat, int *m, int *n, int *mrad, int *nrad, double *outmat);

int CompareDouble(const void *a, const void *b);



