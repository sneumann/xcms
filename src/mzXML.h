void MzXMLOpen(const char *fileName[], int *mxid, int *status);

void MzXMLClose(const int *mxid, int *status);

void MzXMLNumScans(const int *mxid, int *numscans, int *status);

void MzXMLRTTICPeaks(const int *mxid, double *rt, double *tic, 
                     int *peaksTotal, int *status);

void MzXMLSIPeaks(const int *mxid, int *scanindex, double *mz, 
                  double *intensity, int *status);
