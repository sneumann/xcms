/***************************************************************************
                             RAMP


Non sequential parser for mzXML files

                             -------------------
    begin                : Wed Oct 10
    copyright            : (C) 2003 by Pedrioli Patrick, ISB, Proteomics
    email                : ppatrick@student.ethz.ch
    additional work for C++, >2GB files in WIN32, and portability (C) 2004 by Brian Pratt, Insilicos LLC 
 ***************************************************************************/

/***************************************************************************
*																								  *
*	 This program is free software; you can redistribute it and/or modify  *
*	 it under the terms of the GNU Library or "Lesser" General Public 	  *
*	 License (LGPL) as published by the Free Software Foundation;			  *
*	 either version 2 of the License, or (at your option) any later		  *
*	 version.																				  *
*																								  *
***************************************************************************/

#ifndef _RAMP_H
#define _RAMP_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
#define RAP_EXTERN_C extern "C"
#else 
#define RAP_EXTERN_C 
#endif

#ifdef _WIN32
#include <winsock2.h>
#include <sys/types.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
typedef unsigned int uint32_t;
typedef unsigned __int64 uint64_t;
#ifndef strcasecmp
#define strcasecmp stricmp
#endif
#else
#if !defined(sparc) && !defined(__sparc)
#include <stdint.h>
#endif
#include <netinet/in.h>
#endif

#ifndef _LARGEFILE_SOURCE 
typedef int ramp_filehandle_t; // use MSFT API for 64 bit file pointers
#else
typedef FILE * ramp_filehandle_t; // a real OS with real file handling
#endif

// set mz and intensity precision
#ifdef RAMPREAL_DOUBLE
typedef double RAMPREAL; 
#else
typedef float RAMPREAL; 
#endif

//
// we use this struct instead of FILE* so we can track what kind of files we're parsing
//
typedef struct {
   ramp_filehandle_t fileHandle;
   int bIsMzData; // mzXML or mzData?
} RAMPFILE;

#ifndef _LARGEFILE_SOURCE // use MSFT API for 64 bit file pointers
typedef __int64 ramp_fileoffset_t;
#define ramp_fseek(a,b,c) _lseeki64(a->fileHandle,b,c)
#define ramp_ftell(a) _lseeki64(a->fileHandle,0,SEEK_CUR)
#define ramp_fread(buf,len,handle) read(handle->fileHandle,buf,len)
RAP_EXTERN_C char *ramp_fgets(char *buf,int len,RAMPFILE *handle);
#define ramp_feof(handle) eof(handle->fileHandle)
#define atoll(a) _atoi64(a)

#else // a real OS with real file handling
typedef off_t ramp_fileoffset_t;
#define ramp_fseek(a,b,c) fseeko(a->fileHandle,b,c)
#define ramp_ftell(a) ftello(a->fileHandle)
#define ramp_fread(buf,len,handle) fread(buf,1,len,handle->fileHandle)
#define ramp_fgets(buf,len,handle) fgets(buf, len, handle->fileHandle)
#define ramp_feof(handle) feof(handle->fileHandle)
#endif

#include <string.h>
#include <ctype.h>


#ifdef __cplusplus
extern "C" {
#endif

#define INSTRUMENT_LENGTH 2000
#define SCANTYPE_LENGTH 32


struct ScanHeaderStruct
{
   int seqNum; // number in sequence observed file (1-based)
   int acquisitionNum; // scan number as declared in File (may be gaps)
   int  msLevel;
   int  peaksCount;
   double totIonCurrent;
   double retentionTime;        /* in seconds */
   double basePeakMZ;
   double basePeakIntensity;
   double collisionEnergy;
   double ionisationEnergy;
   double lowMZ;
   double highMZ;
   int precursorScanNum; /* only if MS level > 1 */
   double precursorMZ;  /* only if MS level > 1 */
   int precursorCharge;  /* only if MS level > 1 */
   char scanType[SCANTYPE_LENGTH];
   ramp_fileoffset_t filePosition; /* where in the file is this header? */
};

struct RunHeaderStruct
{
  int scanCount;
  double lowMZ;
  double highMZ;
  double startMZ;
  double endMZ;
  double dStartTime;
  double dEndTime;
};

typedef struct InstrumentStruct
{
   char manufacturer[INSTRUMENT_LENGTH];
   char model[INSTRUMENT_LENGTH];
   char ionisation[INSTRUMENT_LENGTH];
   char analyzer[INSTRUMENT_LENGTH];
   char detector[INSTRUMENT_LENGTH];
   //char msType[INSTRUMENT_LENGTH];
} InstrumentStruct;

// file open/close
RAMPFILE *rampOpenFile(const char *filename);
void rampCloseFile(RAMPFILE *pFI);

// construct a filename in buf from a basename, adding .mzXML or .mzData as exists
// returns buf, or NULL if neither .mzXML or .mzData file exists
char *rampConstructInputFileName(char *buf,int buflen,const char *basename);

// trim a filename of its .mzData or .mzXML extension
// return trimmed buffer, or null if no proper .ext found
char *rampTrimBaseName(char *buf);

// locate the .mzData or .mzXML extension in the buffer
// return pointer to extension, or NULL if not found
char *rampValidFileType(char *buf);

// exercise at least some of the ramp interface - return non-0 on failure
int rampSelfTest(char *filename); // if filename is non-null we'll exercise reader with it

ramp_fileoffset_t getIndexOffset(RAMPFILE *pFI);
ramp_fileoffset_t *readIndex(RAMPFILE *pFI,
                ramp_fileoffset_t indexOffset,
                int *iLastScan);
void readHeader(RAMPFILE *pFI,
                ramp_fileoffset_t lScanIndex, // read from this file position
                struct ScanHeaderStruct *scanHeader);
int  readMsLevel(RAMPFILE *pFI,
                 ramp_fileoffset_t lScanIndex);
double readStartMz(RAMPFILE *pFI,
		   ramp_fileoffset_t lScanIndex);
double readEndMz(RAMPFILE *pFI,
		   ramp_fileoffset_t lScanIndex);
int readPeaksCount(RAMPFILE *pFI,
                 ramp_fileoffset_t lScanIndex);
RAMPREAL *readPeaks(RAMPFILE *pFI,
                 ramp_fileoffset_t lScanIndex);
void readRunHeader(RAMPFILE *pFI,
                   ramp_fileoffset_t *pScanIndex,
                   struct RunHeaderStruct *runHeader,
                   int iLastScan);
void readMSRun(RAMPFILE *pFI,
                   struct RunHeaderStruct *runHeader);

int setTagValue(const char* text,
      char* storage,
      int maxlen,
      const char* lead,
      const char* tail);

InstrumentStruct* getInstrumentStruct(RAMPFILE *pFI);

#ifdef __cplusplus
}
#endif

#endif
