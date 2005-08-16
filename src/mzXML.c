#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include "mzXML.h"
#include "ramp.h"

void MzXMLOpen(const char *fileName[], int *mxid, int *status) {

    *mxid = open(fileName[0], O_RDONLY, 0);
	
	*status = (*mxid == -1) ? errno : 0;
}

void MzXMLClose(const int *mxid, int *status) {

    *status = (close(*mxid) == -1) ? errno : 0;
}

void MzXMLNumScans(const int *mxid, int *numscans, int *status) {

    FILE  *mxFile;
	int   i, totalScans;
	off_t indexOffset, *index;
	
	*numscans = 0;
	
	mxFile = fdopen(*mxid, "r");
	
	if (!mxFile) {
	    *status = errno;
		return;
	}
	
	indexOffset = getIndexOffset(mxFile);
	
	if (!indexOffset) {
	   *status = -1;
	   return;
	}
	
	index = readIndex(mxFile, indexOffset, &totalScans);
	
	for (i = 1; i <= totalScans; i++)
	    if (readMsLevel(mxFile, index[i]) == 1)
		    (*numscans)++;
	
	free(index);
	*status = 0;
}

void MzXMLRTTICPeaks(const int *mxid, double *rt, double *tic, 
                     int *peaksTotal, int *status) {

    FILE   *mxFile;
	int    i, totalScans, scani = 0;
	off_t  indexOffset, *index;
	struct ScanHeaderStruct scanHeader;
	
	*peaksTotal = 0;
	
	mxFile = fdopen(*mxid, "r");
	
	if (!mxFile) {
	    *status = errno;
		return;
	}
	
	indexOffset = getIndexOffset(mxFile);
	
	if (!indexOffset) {
	   *status = -1;
	   return;
	}
	
	index = readIndex(mxFile, indexOffset, &totalScans);
	
	for (i = 1; i <= totalScans; i++) {
	    readHeader(mxFile, index[i], &scanHeader);
	    if (scanHeader.msLevel == 1) {
		    rt[scani] = scanHeader.retentionTime;
		    tic[scani] = scanHeader.totIonCurrent;
		    *peaksTotal += scanHeader.peaksCount;
			scani++;
        }
    }
    
    free(index);
    *status = 0;
}

void MzXMLSIPeaks(const int *mxid, int *scanindex, double *mz, 
                  double *intensity, int *status) {

    FILE   *mxFile;
	int    i, totalScans, scani = 0, peaki = 0;
	off_t  indexOffset, *index;
	float  *peaks, *peaksPtr;
	
	mxFile = fdopen(*mxid, "r");
	
	if (!mxFile) {
	    *status = errno;
		return;
	}
	
	indexOffset = getIndexOffset(mxFile);
	
	if (!indexOffset) {
	   *status = -1;
	   return;
	}
	
	index = readIndex(mxFile, indexOffset, &totalScans);

    for (i = 1; i <= totalScans; i++)
	    if (readMsLevel(mxFile, index[i]) == 1) {
		    scanindex[scani++] = peaki;
		    peaks = readPeaks(mxFile, index[i]);
		    if (!peaks) {
		        *status = -1;
		        free(index);
		        return;
		    }
		    peaksPtr = peaks;
		    while (*peaksPtr >= 0) {
		        mz[peaki] = *(peaksPtr++);
		        intensity[peaki] = *(peaksPtr++);
		        peaki++;
		    }
		    free(peaks);
		}
    
    free(index);
    *status = 0;
}
