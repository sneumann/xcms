/***************************************************************************
                             RAMP


Non sequential parser for mzXML files
and mzData files, too!

                             -------------------
    begin                : Wed Oct 10
    copyright            : (C) 2003 by Pedrioli Patrick, ISB, Proteomics
    email                : ppatrick@student.ethz.ch
    additional work for C++, >2GB files in WIN32, and portability (C) 2004 by Brian Pratt, Insilicos LLC 
    additional work for mzData input (C) 2005 Brian Pratt Insilicos LLC
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

#include "ramp.h"

#include "base64.h"

#include "time.h"

#define SIZE_BUF 512

#if defined __LITTLE_ENDIAN
#define swapbytes(x) ntohl(x)  /* use system byteswap (ntohl is a noop on bigendian systems) */
#else
uint32_t swapbytes(uint32_t x) {
  return ((x & 0x000000ffU) << 24) | 
         ((x & 0x0000ff00U) <<  8) | 
         ((x & 0x00ff0000U) >>  8) | 
         ((x & 0xff000000U) >> 24);
}
#endif

uint64_t swapbytes64(uint64_t x) {
  return ((((uint64_t)swapbytes((uint32_t)(x & 0xffffffffU)) << 32) | 
            (uint64_t)swapbytes((uint32_t)(x >> 32))));
}

/****************************************************************
 * Utility functions					*
 ***************************************************************/

const char* skipspace(const char* pStr)
{
	while (isspace(*pStr))
		pStr++;
	if (*pStr == '\0')
		return NULL;
	return pStr;
}

static void getIsLittleEndian(const char *buf, int *result) {
   const char *p = strstr(buf,"byteOrder");
   if (p) {
      p = strchr(p,'\"');
      if (p++) {
         *result = (0!=strncmp(p,"network",7));
      }
   }
   // leave *result alone if we don't see anything!
}

/**************************************************
* open and close files *
**************************************************/
RAMPFILE *rampOpenFile(const char *filename) {
   RAMPFILE *result = (RAMPFILE *)calloc(1,sizeof(RAMPFILE));
   if (result) {
      int bOK;
#ifndef _LARGEFILE_SOURCE
     result->fileHandle = open(filename,_O_BINARY|_O_RDONLY);
     bOK = (result->fileHandle >= 0);
#else
     result->fileHandle = fopen(filename,"rb");
     bOK = (result->fileHandle != NULL);
#endif
     if (!bOK) {
        free(result);
        result = NULL;
     } else {
        char buf[1024];
        buf[sizeof(buf)-1] = 0;
        while (!ramp_feof(result)) {
           ramp_fgets(buf,sizeof(buf)-1,result);
           if (strstr(buf,"<msRun")) {
              result->bIsMzData = 0;
              break;
           } else if (strstr(buf,"<mzData")) {
              result->bIsMzData = 1;
              break;
           }
        }
        ramp_fseek(result,0,SEEK_SET); // rewind
        // now set up the element tags we seek
     }
   }
   return result;
}

void rampCloseFile(RAMPFILE *pFI) {
   if (pFI) {
#ifndef _LARGEFILE_SOURCE
      close(pFI->fileHandle);
#else
      fclose(pFI->fileHandle);
#endif
      free(pFI);
   }
}

/**************************************************
* fgets() for win32 long files *
* TODO: this could be a LOT more efficient...
**************************************************/
#ifndef _LARGEFILE_SOURCE
char *ramp_fgets(char *buf,int len,RAMPFILE *handle) {
   int nread=0;
   int chunk;
   ramp_fileoffset_t pos = ramp_ftell(handle);
   buf[--len]=0; // nullterm for safety
   chunk = max(len/4,1); // usually all that's needed is a short read
   while (nread <= len) {
      char *newline;
      int nread_now = ramp_fread(buf+nread,chunk,handle); 
      if (!nread_now) {
         return nread?buf:NULL;
      }
      newline = strchr(buf+nread,'\n');
      if (newline) {
         *newline = 0;
         ramp_fseek(handle,pos+(newline-buf)+1,SEEK_SET); // so next read is at next line
         break;
      }
      nread+=nread_now;
      if (nread >= len) {
         break;
      }
      // apparently we need bigger reads
      chunk = len-nread;
   }
   return buf;
}
#endif

/****************************************************************
 * Find the Offset of the index					*
 ***************************************************************/
ramp_fileoffset_t getIndexOffset(RAMPFILE *pFI)
{
   int  i;
   ramp_fileoffset_t  indexOffset, indexOffsetOffset;
   char indexOffsetTemp[SIZE_BUF], buf;

   if (pFI->bIsMzData) {
      return -1; // no index in mzData
   }

   for (indexOffsetOffset = -120;  indexOffsetOffset++ < 0 ;)
   {
      char seekbuf[SIZE_BUF];
      char *target = "<indexOffset>";
      int  nread;

      ramp_fseek(pFI, indexOffsetOffset, SEEK_END);
      nread = ramp_fread(seekbuf, strlen(target), pFI);
      seekbuf[nread] = '\0';
      
      if (!strcmp(seekbuf, target))
      {
         break;
      }
   }

   if (indexOffsetOffset >= 0) {
      return -1; // no answer
   }

   indexOffset = ramp_ftell(pFI);

   i = 0;
   while (ramp_fread(&buf, 1, pFI) && buf != '<')
   {
      indexOffsetTemp[i] = buf;
      i++;
   }
   indexOffsetTemp[i] = '\0';

   if (sizeof(ramp_fileoffset_t)==8)
      indexOffset = (atoll(indexOffsetTemp));
   else
      indexOffset = (atol(indexOffsetTemp));
   // now test this
   ramp_fseek(pFI, indexOffset, SEEK_SET);
   ramp_fread(indexOffsetTemp, sizeof(indexOffsetTemp), pFI);
   indexOffsetTemp[sizeof(indexOffsetTemp)-1] = 0;
   if (!strstr(indexOffsetTemp,"<index")) {
      indexOffset = -1; // broken somehow
   }
   return indexOffset;
}


/****************************************************************
 * Reads the Scan index in a list				*
 * Returns pScanIndex which becomes property of the caller	*
 * pScanIndex is -1 terminated					*
 ***************************************************************/
char buf[SIZE_BUF*16];

ramp_fileoffset_t *readIndex(RAMPFILE *pFI,
                ramp_fileoffset_t indexOffset,
                int *iLastScan)
{
   int  n;
   int  reallocSize = 8000;    /* initial # of scan indexes to expect */
   char *beginScanOffset;
   ramp_fileoffset_t *pScanIndex=NULL;
   int retryLoop;
   
   for (retryLoop = 2;retryLoop--;) {
     n = 1; // ramp is one based
     *iLastScan = 0;
     free(pScanIndex);
      if ((indexOffset < 0) || (retryLoop==0)) { // derive the index by inspection
         // no index found, derive it
         const char *scantag = pFI->bIsMzData?"<spectrum id":"<scan";
         int taglen = (int)strlen(scantag);
         ramp_fileoffset_t index = 0;
         pScanIndex = (ramp_fileoffset_t *)malloc( sizeof(ramp_fileoffset_t)*reallocSize); // allocate space for the scan index info
         if (!pScanIndex) {
            printf("Cannot allocate memory\n");
            return NULL;
         }
         ramp_fseek(pFI,0,SEEK_SET);
         buf[sizeof(buf)-1] = 0;
         while ((int)ramp_fread(buf,sizeof(buf)-1,pFI)>taglen) {
            char *find;
            char *look=buf;
            int nread;
            while (NULL != (find = strstr(look,scantag))) {
               pScanIndex[n++] = index+(find-buf); // ramp is 1-based
               (*iLastScan)++;
               look = find+taglen;
               if (reallocSize<=n) {
                  reallocSize+=500;
                  pScanIndex = (ramp_fileoffset_t *)realloc(pScanIndex, sizeof(ramp_fileoffset_t)*reallocSize); // allocate space for the scan index info
                  if (!pScanIndex) {
                     printf("Cannot allocate memory\n");
                     return NULL;
                  }
               }
            }
            nread = strlen(look)+(look-buf);
            if (*look && strchr(scantag,buf[nread-1])) { // check last char of buffer
               // possible that next scantag overhangs end of buffer
               ramp_fseek(pFI,-taglen,SEEK_CUR); // so next get includes it
            } 
            index = ramp_ftell(pFI);
         }
         break; // no need to retry
      } else {  // read the index out of the file (then check it!)
         struct ScanHeaderStruct scanHeader; // for test of index integrity
         int indexOK=1; // until we show otherwise

         if ((pScanIndex = (ramp_fileoffset_t *) malloc(reallocSize * sizeof(ramp_fileoffset_t))) == NULL) {
            printf("Cannot allocate memory\n");
            return NULL;
         }
         
         ramp_fseek(pFI, indexOffset, SEEK_SET);
         
         ramp_fgets(buf, SIZE_BUF, pFI);
         while( !strstr( buf , "<offset" ) ) {
            ramp_fgets(buf, SIZE_BUF, pFI);
         }
         
         while (!strstr(buf, "/index")) {
            if ((beginScanOffset = (char *) (strstr(buf, ">"))) == NULL) {
               ramp_fgets(buf, SIZE_BUF, pFI);
               continue;
            }
            beginScanOffset++;
            
            if (sizeof(ramp_fileoffset_t)==8) {
               pScanIndex[n] = atoll(beginScanOffset);
            } else {
               pScanIndex[n] = atol(beginScanOffset);
            }
            n++;
            (*iLastScan)++;
            
            if (n == reallocSize) {
               ramp_fileoffset_t *pTmp;
               
               reallocSize = reallocSize + 500;
               
               pTmp = (ramp_fileoffset_t*) realloc(pScanIndex, reallocSize*sizeof(ramp_fileoffset_t));
               
               if (pTmp==NULL) {
                  printf("Cannot allocate memory\n");
                  return NULL;
               } else {
                  pScanIndex=pTmp;
               }
            }
            
            ramp_fgets(buf, SIZE_BUF, pFI);
         }
         if (n) {
            // OK, now test that to see if it's a good index or not
            readHeader(pFI,pScanIndex[1],&scanHeader);
            if (scanHeader.acquisitionNum == -1) { // first
               indexOK = 0; // bogus index
               free(pScanIndex);
               pScanIndex = NULL;
            } 
            if (indexOK && (n>1)) { // middle
               readHeader(pFI,pScanIndex[n/2],&scanHeader);
               if (scanHeader.acquisitionNum == -1) {
                  indexOK = 0; // bogus index
                  free(pScanIndex);
                  pScanIndex = NULL;
               }
            }
            if (indexOK && (n>1)) { // last
               readHeader(pFI,pScanIndex[n-1],&scanHeader);
               if (scanHeader.acquisitionNum == -1) {
                  indexOK = 0; // bogus index
                  free(pScanIndex);
                  pScanIndex = NULL;
               }
            }
         }
         if (indexOK) {
            break; // no retry
         }
      } // end if we claim to have an index
   } // end for retryloop
   pScanIndex[n] = -1;

   return (pScanIndex);
}

// helper func for reading mzData
const char *findMzDataTagValue(const char *pStr, const char *tag) {
   const char *find = strstr(pStr,tag);
   if (find) {
      find = strstr(find+1,"value=");
      if (find) {
         find = strchr(find,'\"');
         if (find) {
            find++; // pointing at value string
         }
      }
   }
   return find;
}

/*
 * Reads a time string, returns time in seconds.
 */
static double rampReadTime(RAMPFILE *pFI,const char *pStr) {
   double t=0;
   if (pFI->bIsMzData) {
      const char *tag = findMzDataTagValue(pStr, "TimeInMinutes");
      if (tag) {
         t = 60.0*atof(tag);
      } else if (NULL!=(tag = findMzDataTagValue(pStr, "TimeInSeconds"))) { // von Steffan Neumann
         t = atof(tag);
      }
   } else if (!sscanf(pStr, "PT%lfS", &t)) {  // usually this is elapsed run time
      /* but could be stored in for PxYxMxDTxHxMxS */
      struct tm fullTime; // apologies to those working after January 18, 19:14:07, 2038
      double secondsFrac=0;
      int bDate = 1; 
      while (*++pStr != '\"') {
         double val;
         if ('T'==*pStr) {
            pStr++;
            bDate = 0; // we're into the minutes:seconds portion
         }
         val = atof(pStr);
         while (('.'==*pStr)||isdigit(*pStr)) {
            pStr++;
         }
         switch (*pStr) {
         case 'Y':
            fullTime.tm_year = (int)val-1900; // years since 1900
            break;
         case 'M':
            if (bDate) {
               fullTime.tm_mon = (int)val-1; // range 0-11
            } else {
               fullTime.tm_min = (int)val; // range 0-59
            }            
            break;
         case 'D':
            fullTime.tm_mday = (int)val; // range 1-31
            break;
         case 'H':
            fullTime.tm_hour = (int)val; // range 0-23
            break;
         case 'S':
            fullTime.tm_sec = (int)val;
            secondsFrac = val-(double)fullTime.tm_sec;
            break;
         }
      }
      t = (double)mktime(&fullTime)+secondsFrac;
   }
   return t;
}

/*
 * helper func for faster parsing
 */
const char *matchAttr(const char *where,const char *attr,int len) {
   const char *look = where; // we assume this is pointed at '=', look back at attr
   while (len--) {
      if (*--look != attr[len]) {
         return NULL; // no match
      }
   }
   return where+2; // point past ="
}


/*
 * Reads scan header information.
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE
 *    RETURNING !
 */
void readHeader(RAMPFILE *pFI,
                ramp_fileoffset_t lScanIndex, // look here
                struct ScanHeaderStruct *scanHeader)
{
   char stringBuf[SIZE_BUF];
   char *pStr2;

   ramp_fseek(pFI, lScanIndex, SEEK_SET);

   /*
    * initialize defaults
    */
   memset(scanHeader,0,sizeof(struct ScanHeaderStruct)); // mostly we want 0's
#define LOWMZ_UNINIT 1.111E6
   scanHeader->lowMZ =  LOWMZ_UNINIT;
   scanHeader->acquisitionNum = -1;
   scanHeader->seqNum = -1;


   if (pFI->bIsMzData) {
      int bHasPrecursor = 0;
      while (ramp_fgets(stringBuf, SIZE_BUF, pFI))
      {
         const char *attrib;
         const char *pStr;
         // find each attribute in stringBuf
         for (attrib=stringBuf-1;NULL!=(attrib=strchr(attrib+1,'='));) {
            if ((pStr = matchAttr(attrib, "spectrum id",11))) {
               sscanf(pStr, "%d\"", &(scanHeader->acquisitionNum));
            //} else if ((pStr = matchAttr(attrib, "basePeakMz",10)))  {
            //   sscanf(pStr, "%lf\"", &(scanHeader->basePeakMZ));      
            //} else if ((pStr = matchAttr(attrib, "totIonCurrent",13)))  {
            //   sscanf(pStr, "%lf\"", &(scanHeader->totIonCurrent));      
            //} else if ((pStr = matchAttr(attrib, "basePeakIntensity",17)))  {
            //   sscanf(pStr, "%lf\"", &(scanHeader->basePeakIntensity));      
            } else if ((pStr = matchAttr(attrib, "msLevel",7)))  {
               sscanf(pStr, "%d\"", &(scanHeader->msLevel));
            //} else if ((pStr = matchAttr(attrib, "length",6)))  { get this from array element
            //   sscanf(pStr, "%d\"", &(scanHeader->peaksCount));
            } else if ((pStr = findMzDataTagValue(attrib,"TimeInMinutes")))  {
               scanHeader->retentionTime = rampReadTime(pFI,stringBuf);
            } else if ((pStr = matchAttr(attrib, "mzRangeStart",12)))  {
               sscanf(pStr, "%lf\"", &(scanHeader->lowMZ));
            } else if ((pStr = matchAttr(attrib, "mzRangeStop",11)))  {
               sscanf(pStr, "%lf\"", &(scanHeader->highMZ));
            } else if ((pStr = findMzDataTagValue(attrib, "ScanMode"))) { 
               if ((pStr2 = (char *) strchr(pStr, '\"'))) {
                  memcpy(&(scanHeader->scanType), pStr, pStr2-pStr);
                  scanHeader->scanType[pStr2-pStr] = '\0';
               }
            }
         }
         if (strstr(stringBuf, "</spectrumSettings>")) {
            break; // into data territory now
         }
      }
      do {
         
         /*
         * read precursor info
         */
         const char *pStr,*pStr2;
         if ((pStr = (char *) strstr(stringBuf, "<precursorList")))
         {
            bHasPrecursor = 1;
         }
         if ((pStr = (char *) strstr(stringBuf, "</precursorList")))
         {
            bHasPrecursor = 0;
         }
         if (bHasPrecursor) { // in precursor section
            if (NULL!=(pStr2 = (char *) strstr(stringBuf, "spectrumRef=\""))) 
            {
               sscanf(pStr2 + 13, "%d\"", &(scanHeader->precursorScanNum));
            }
            if (NULL!=(pStr2 = findMzDataTagValue(stringBuf,"ChargeState"))) 
            {
               scanHeader->precursorCharge = atoi(pStr2);
            }
            
            if (NULL!=(pStr2 = findMzDataTagValue(stringBuf,"MassToChargeRatio"))) 
            {
               scanHeader->precursorMZ = atoi(pStr2);
            }
            if (NULL!=(pStr2 = findMzDataTagValue(stringBuf,"mz"))) 
            {
               scanHeader->precursorMZ = atoi(pStr2);
            }
         }
         if (strstr(stringBuf, "</spectrumDesc>")) {
            break; // into data territory now
         }
         if (strstr(stringBuf, "</precursorList>")) {
            break; // into data territory now
         }
      } while (ramp_fgets(stringBuf, SIZE_BUF, pFI));
      // now read peaks count
      do {
         if (strstr(stringBuf, "ArrayBinary>")) {
            do {
               char *cp=strstr(stringBuf,"length=");
               if (cp) {
                  scanHeader->peaksCount = atoi(cp+8);
                  break;
               }
            } while (ramp_fgets(stringBuf, SIZE_BUF, pFI));
            break;
         }
      } while (ramp_fgets(stringBuf, SIZE_BUF, pFI));
   } else { // mzXML
      while (ramp_fgets(stringBuf, SIZE_BUF, pFI))
      {
         const char *attrib;
         const char *pStr;
         // find each attribute in stringBuf
         for (attrib=stringBuf-1;NULL!=(attrib=strchr(attrib+1,'='));) {
            if ((pStr = matchAttr(attrib, "num",3))) {
               sscanf(pStr, "%d\"", &(scanHeader->acquisitionNum));
            } else if ((pStr = matchAttr(attrib, "basePeakMz",10)))  {
               sscanf(pStr, "%lf\"", &(scanHeader->basePeakMZ));      
            } else if ((pStr = matchAttr(attrib, "totIonCurrent",13)))  {
               sscanf(pStr, "%lf\"", &(scanHeader->totIonCurrent));      
            } else if ((pStr = matchAttr(attrib, "basePeakIntensity",17)))  {
               sscanf(pStr, "%lf\"", &(scanHeader->basePeakIntensity));      
            } else if ((pStr = matchAttr(attrib, "collisionEnergy",15)))  {
               sscanf(pStr, "%lf\"", &(scanHeader->collisionEnergy));
            } else if ((pStr = matchAttr(attrib, "msLevel",7)))  {
               sscanf(pStr, "%d\"", &(scanHeader->msLevel));
            } else if ((pStr = matchAttr(attrib, "peaksCount",10)))  {
               sscanf(pStr, "%d\"", &(scanHeader->peaksCount));
            } else if ((pStr = matchAttr(attrib, "retentionTime",13)))  {
               scanHeader->retentionTime = rampReadTime(pFI,pStr);
            } else if ((pStr = matchAttr(attrib, "lowMz",5)))  {
               sscanf(pStr, "%lf\"", &(scanHeader->lowMZ));
            } else if ((pStr = matchAttr(attrib, "highMz",6)))  {
               sscanf(pStr, "%lf\"", &(scanHeader->highMZ));
            } else if ((scanHeader->lowMZ==LOWMZ_UNINIT) &&  
               ((pStr = matchAttr(attrib, "startMz",7)))) {
               sscanf(pStr, "%lf\"", &(scanHeader->lowMZ));  
            } else if ((!scanHeader->highMZ) &&
               ((pStr = matchAttr(attrib, "endMz",5)))) {
               sscanf(pStr, "%lf\"", &(scanHeader->highMZ));
            } else if ((pStr = matchAttr(attrib, "scanType", 8))) { 
               if ((pStr2 = (char *) strchr(pStr, '\"'))) {
                  memcpy(&(scanHeader->scanType), pStr, sizeof(char)*((pStr2-pStr)));
                  scanHeader->scanType[pStr2-pStr] = '\0';
               }
            } else if ((pStr = matchAttr(attrib, "collisionEnergy", 15)))  {
                 sscanf(pStr, "%lf\"", &(scanHeader->collisionEnergy));
            }
         }
         
         /*
         * read precursor mass
         */
         if ((pStr = (char *) strstr(stringBuf, "<precursorMz ")))
         {
            if ((pStr2 = (char *) strstr(stringBuf, "precursorScanNum=\""))) 
            {
               sscanf(pStr2 + 18, "%d\"", &(scanHeader->precursorScanNum));
            }
            
            /*
            * Check for precursor charge.
            */
            if ((pStr2 = (char *) strstr(pStr, "precursorCharge=\""))) 
            {
               sscanf(pStr2 + 17, "%d\"", &(scanHeader->precursorCharge));
            }
            
            /*
            * Find end of tag.
            */
            while (!(pStr = strchr(pStr, '>')))
            {      
               ramp_fgets(stringBuf, SIZE_BUF, pFI);
               pStr = stringBuf;
               if ((pStr2 = (char *) strstr(stringBuf, "precursorScanNum=\"")))
                  sscanf(pStr2 + 18, "%d\"", &(scanHeader->precursorScanNum));
               if ((pStr2 = (char *) strstr(stringBuf, "precursorCharge=\""))) 
               {
                  sscanf(pStr2 + 17, "%d\"", &(scanHeader->precursorCharge));
               }
            }
            pStr++;	// Skip >
            
                     /*
                     * Skip past white space.
            */
            while (!(pStr = skipspace(pStr)))
            {
               ramp_fgets(stringBuf, SIZE_BUF, pFI);
               pStr = stringBuf;
               if ((pStr2 = (char *) strstr(stringBuf, "precursorScanNum=\""))) {
                  sscanf(pStr2 + 18, "%d\"", &(scanHeader->precursorScanNum));
               }
               if ((pStr2 = (char *) strstr(stringBuf, "precursorCharge=\""))) 
               {
                  sscanf(pStr2 + 17, "%d\"", &(scanHeader->precursorCharge));
               }
            }
            
            sscanf(pStr, "%lf<", &(scanHeader->precursorMZ));
            //         printf("precursorMass = %lf\n", scanHeader->precursorMZ);
         }
         if (strstr(stringBuf, "<peaks")) {
            break; // into data territory now
         }
      }
   } // end else mzXML
   if (!scanHeader->retentionTime) { 
      scanHeader->retentionTime = scanHeader->acquisitionNum; // just some unique nonzero value
   }
   scanHeader->seqNum = scanHeader->acquisitionNum; //  default sequence number
   scanHeader->filePosition = lScanIndex;
}

/****************************************************************
 * Reads the MS level of the scan.				*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

int readMsLevel(RAMPFILE *pFI,
      ramp_fileoffset_t lScanIndex)
{
   int  msLevelLen;
   char stringBuf[SIZE_BUF];
   char szLevel[12];
   char *beginMsLevel, *endMsLevel;

   ramp_fseek(pFI, lScanIndex, SEEK_SET);

   ramp_fgets(stringBuf, SIZE_BUF, pFI);

   while (!(beginMsLevel = (char *) strstr(stringBuf, "msLevel=\"")))
   {
      ramp_fgets(stringBuf, SIZE_BUF, pFI);
   }

   beginMsLevel += 9;           // We need to move the length of msLevel="
   endMsLevel = (char *) strstr(beginMsLevel, "\"");
   msLevelLen = endMsLevel - beginMsLevel;

   strncpy(szLevel, beginMsLevel, msLevelLen);
   szLevel[msLevelLen] = '\0';

   return atoi(szLevel);
}


/****************************************************************
 * Reads startMz and endMz of the scan.				*
 * Returns 1.E6 if startMz was not set.*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

double readStartMz(RAMPFILE *pFI,
		   ramp_fileoffset_t lScanIndex)
{
  char stringBuf[SIZE_BUF];
  double startMz = 1.E6;
  char *pStr;
  const char *tag = pFI->bIsMzData?"mzRangeStart":"startMz";

  ramp_fseek(pFI, lScanIndex, SEEK_SET);
   
  while (ramp_fgets(stringBuf, SIZE_BUF, pFI))
  {
     if (strstr(stringBuf, pFI->bIsMzData?"</spectrumDesc":"<peaks"))
        break; // ran to end
      if ((pStr = strstr(stringBuf, tag))){
        sscanf(pStr + strlen(tag)+2, "%lf\"", &startMz);
        break;
      }
  }

  return startMz;
}


/****************************************************************
 * Reads startMz and endMz of the scan.				*
 * Returns 1.E6 if startMz was not set.	*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

double readEndMz(RAMPFILE *pFI,
		   ramp_fileoffset_t lScanIndex)
{
  char stringBuf[SIZE_BUF];
  double endMz = 0.0;
  char *pStr;
  const char *tag = pFI->bIsMzData?"mzRangeStop":"endMz";

  ramp_fseek(pFI, lScanIndex, SEEK_SET);
   
  while (ramp_fgets(stringBuf, SIZE_BUF, pFI))
  {
     if (strstr(stringBuf, pFI->bIsMzData?"</spectrumDesc":"<peaks"))
        break; // ran to end
      if ((pStr = strstr(stringBuf, tag))){
        sscanf(pStr + strlen(tag)+2, "%lf\"", &endMz);
        break;
      }
  }

  return endMz;
}

#define N_EXT_TYPES 4
static char *data_ext[N_EXT_TYPES]={".mzXML",".mzData",".mzxml",".mzdata"};


// construct a filename in inbuf from a basename, adding .mzXML or .mzData as exists
// returns inbuf, or NULL if neither .mzXML or .mzData file exists
char *rampConstructInputFileName(char *inbuf,int inbuflen,const char *basename_in) {
   char *result = NULL;
   int i;
   char *basename = (char *)basename_in;
   char *tmpbuf = (char *)malloc(inbuflen);
   if (basename_in==inbuf) { // same pointer
      basename = (char *)malloc(inbuflen);
      strncpy(basename,basename_in,inbuflen);
   }

   for (i=0;i<N_EXT_TYPES;i++) {
      FILE *test;
      strncpy(tmpbuf,basename,inbuflen);
      tmpbuf[inbuflen-1]='\0';
      strncat(tmpbuf,data_ext[i],inbuflen-strlen(tmpbuf));
      test = fopen(tmpbuf,"r");
      if (test != NULL) {
         if (result) { // conflict! both mzXML and mzData are present
            if (strcasecmp(tmpbuf,result)) { // win32 isn't case sensitive
               printf("found both %s and %s, using %s\n",
                  tmpbuf,result,result);
            }
         } else { // found file, copy the constructed name
            strncpy(inbuf,tmpbuf,inbuflen);
            result = inbuf;
         }
         fclose(test);
      }
   }
   if (!result) { // failed - caller can complain about lack of .mzXML
      strncpy(inbuf,basename,inbuflen);
      strncat(inbuf,data_ext[0],inbuflen-strlen(tmpbuf));
   }
   if (basename_in==inbuf) { // same pointer
      free(basename);
   }
   free(tmpbuf);
   return result;
}

// return NULL if fname is not of extension type we handle,
// otherwise return pointer to .ext
char *rampValidFileType(char *fname) {
   char *result;
   int i;
   for (i = N_EXT_TYPES;i--;) {
      result = strrchr(fname,'.');
      if (!strcasecmp(result,data_ext[i])) {
         break;
      }
      result = NULL; // try again
   }
   return result;
}

// remove the filename .ext, if found
// return NULL if no .ext found, else return fname
char *rampTrimBaseName(char *fname) {
   char *ext = rampValidFileType(fname);
   if (ext) {
      *ext = 0; // trim the extension
   }
   return ext?fname:NULL;
}


/****************************************************************
 * READS the number of peaks.			*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/
int readPeaksCount(RAMPFILE *pFI,
      ramp_fileoffset_t lScanIndex)
{
   char stringBuf[SIZE_BUF];
   char *beginPeaksCount, *peaks;
   int result = 0;
   const char *tag = pFI->bIsMzData?"length=\"":"peaksCount=\"";
   ramp_fileoffset_t in_lScanIndex = lScanIndex;

   ramp_fseek(pFI, lScanIndex, SEEK_SET);

   // Get the num of peaks in the scan and allocate the space we need
   ramp_fgets(stringBuf, SIZE_BUF, pFI);
   while (!(beginPeaksCount = (char *) strstr(stringBuf, tag)))
   {
      lScanIndex = ramp_ftell(pFI);
      ramp_fgets(stringBuf, SIZE_BUF, pFI);
   }

   // We need to move forward the length of the tag
   beginPeaksCount += strlen(tag);
   result = atoi(beginPeaksCount);

   // mext call is probably to read the <peaks> section, position there if needed
   if (pFI->bIsMzData) {
      ramp_fseek(pFI, in_lScanIndex, SEEK_SET);
   } else {
      peaks = strstr(stringBuf,"<peaks");
      if (peaks) {
         ramp_fseek(pFI, lScanIndex+(peaks-stringBuf), SEEK_SET);
      }
   }
   return result;
}


/****************************************************************
 * READS the base64 encoded list of peaks.			*
 * Return a RAMPREAL* that becomes property of the caller!		*
 * The list is terminated by -1					*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

RAMPREAL *readPeaks(RAMPFILE *pFI,
      ramp_fileoffset_t lScanIndex)
{
   int  n=0;
   int  peaksCount=0;
   int  peaksLen;       // The length of the base64 section
   int precision = 0;
   RAMPREAL *pPeaks = NULL;

   int  endtest = 1;
   int  weAreLittleEndian = *((char *)&endtest);

   char *pData = NULL;
   char *pBeginData;
   char *pDecoded = NULL;

   char buf[1000];
   buf[sizeof(buf)-1] = 0;

   if (pFI->bIsMzData) {
      // intensity and mz are written in two different arrays
      int bGotInten = 0;
      int bGotMZ = 0;
      ramp_fseek(pFI,lScanIndex,SEEK_SET);
      while ((!(bGotInten && bGotMZ)) &&
           ramp_fgets(buf,sizeof(buf)-1,pFI)) {
         int isArray = 0;
         int isInten = 0;
         int isLittleEndian = 0;
         int byteOrderOK;
         if (strstr(buf,"<mzArrayBinary")) {
            isArray = bGotMZ = 1;
         } else if (strstr(buf,"<intenArrayBinary")) {
            isArray = isInten = bGotInten = 1;
         }
         if (isArray) {
            int partial,triplets,bytes;
            // now determine peaks count, precision
            while (!(pBeginData = (char *) strstr(buf, "<data")))
            {
               ramp_fgets(buf, sizeof(buf)-1, pFI);
            }
            
            // find precision="xx"
            if( !(pBeginData = strstr( buf , "precision=" )))
            {
               precision = 32; // default value
            } else { // we found declaration
               precision = atoi(strchr(pBeginData,'\"')+1);
            }
            
            // find length="xx"
            if((!peaksCount) && (pBeginData = strstr( buf , "length=" )))
            {
               peaksCount = atoi(strchr(pBeginData,'\"')+1);
            }

            if (peaksCount <= 0)
            { // No peaks in this scan!!
               return NULL;
            }

            // find endian="xx"
            if((pBeginData = strstr( buf , "endian=" )))
            {
               isLittleEndian = !strncmp("little",strchr(pBeginData,'\"')+1,6);
            }
            
            // find close of <data>
            while( !(pBeginData = strstr( buf , ">" )))
            {
               ramp_fgets(buf, sizeof(buf)-1 , pFI);
            }
            pBeginData++;	// skip the >
            
            // base64 has 4:3 bloat, precision/8 bytes per value
            bytes = (peaksCount*(precision/8));
            // for every 3 bytes base64 emits 4 characters - 1, 2 or 3 byte input emits 4 bytes
            triplets = (bytes/3)+((bytes%3)!=0);
            peaksLen = 4*triplets;
            
            if ((pData = (char *) realloc(pData,1 + peaksLen)) == NULL)
            {
               printf("Cannot allocate memory\n");
               return NULL;
            }
            
            // copy in any partial read of peak data, and complete the read
            strncpy(pData,pBeginData,peaksLen);
            pData[peaksLen] = 0;
            partial = strlen(pData);
            if (partial < peaksLen) {              
               ramp_fread(pData+partial,peaksLen-partial, pFI);
            }
                        
            if ((pDecoded = (char *) realloc(pDecoded,peaksCount * (precision/8) + 1)) == NULL)
            {
               printf("Cannot allocate memory\n");
               return NULL;
            }
            // Base64 decoding
            b64_decode_mio(pDecoded, pData);
            
            if ((!pPeaks) && ((pPeaks = (RAMPREAL *) malloc((peaksCount+1) * 2 * sizeof(RAMPREAL) + 1)) == NULL))
            {
               printf("Cannot allocate memory\n");
               return NULL;
            }
            
            // And byte order correction
            byteOrderOK = (isLittleEndian==weAreLittleEndian);
            if (32==precision) { // floats
               if (byteOrderOK) {
                  float *f = (float *) pDecoded;
                  for (n = 0; n < peaksCount; n++) {
                     pPeaks[isInten+(2*n)] = (RAMPREAL)*f++;
                  }
               } else {
                  uint32_t *u = (uint32_t *) pDecoded;
                  for (n = 0; n < peaksCount; n++) {
                     uint32_t tmp = swapbytes( *u++ );
                     pPeaks[isInten+(2*n)] = (RAMPREAL) (* (float *)(&tmp));
                  }
               }
            } else { // doubles
               if (byteOrderOK) {
                  double *d = (double *)pDecoded;
                  for (n = 0; n < peaksCount; n++) {
                     pPeaks[isInten+(2*n)] = (RAMPREAL)*d++;
                  }
               } else {
                  uint64_t *u = (uint64_t *) pDecoded;
                  for (n = 0; n < peaksCount; n++) {
                     uint64_t tmp = swapbytes64( *u++ );
                     pPeaks[isInten+(2*n)] = (RAMPREAL) (* (double *)(&tmp));
                  }
               }
            }
            if (bGotInten && bGotMZ) {
               break;
            }
         } // end if isArray
      } // end while we haven't got both inten and mz
      free(pData);
      free(pDecoded);
      pPeaks[peaksCount*2] = -1; // some callers want a terminator
   } else { // mzXML
      int partial,bytes,triplets;
      int isLittleEndian = 0; // default is network byte order (Big endian)
      int byteOrderOK;
      peaksCount = readPeaksCount(pFI, lScanIndex);
      if (peaksCount <= 0)
      { // No peaks in this scan!!
         return NULL;
      }
      
      // now determine peaks precision
      ramp_fgets(buf, sizeof(buf)-1, pFI);
      while (!(pBeginData = (char *) strstr(buf, "peaks")))
      {
         ramp_fgets(buf, sizeof(buf)-1, pFI);
      }
      getIsLittleEndian(buf,&isLittleEndian);

      // find precision="xx"
      while( !(pBeginData = strstr( buf , "precision=" )))
      {        
         if (NULL!=(pBeginData = strstr( buf , ">" ))) {
            precision = 32; // default value
            break; // no precision declared
         }
         ramp_fgets(buf, sizeof(buf)-1 , pFI);
         getIsLittleEndian(buf,&isLittleEndian); // in case it wasn't in previous line
      }
      if (!precision) { // we didn't hit end of element, ie we found declaration
         precision = atoi(strchr(pBeginData,'\"')+1);
      }
 
      // find close of <peaks precision="xx">
      while( !(pBeginData = strstr( buf , ">" )))
      {
         ramp_fgets(buf, sizeof(buf)-1 , pFI);
         getIsLittleEndian(buf,&isLittleEndian); // in case it wasn't in previous line
      }
      pBeginData++;	// skip the >
      
      // base64 has 4:3 bloat, precision/8 bytes per value, 2 values per peak
      bytes = (peaksCount*(precision/4));
      // for every 3 bytes base64 emits 4 characters - 1, 2 or 3 byte input emits 4 bytes
      triplets = (bytes/3)+((bytes%3)!=0);
      peaksLen = 4*triplets;
      
      if ((pData = (char *) malloc(1 + peaksLen)) == NULL)
      {
         printf("Cannot allocate memory\n");
         return NULL;
      }
      pData[peaksLen] = 0;
      
      // copy in any partial read of peak data, and complete the read
      strncpy(pData,pBeginData,peaksLen);
      partial = strlen(pData);
      if (partial < peaksLen) {
         ramp_fread(pData+partial,peaksLen-partial, pFI);
      }
      
      // 2 values per peak, precision/8 bytes per value
      if ((pDecoded = (char *) malloc(peaksCount * (precision/4) + 1)) == NULL)
      {
         printf("Cannot allocate memory\n");
         return NULL;
      }
      // Base64 decoding
      b64_decode_mio(pDecoded, pData);
      free(pData);
      
      if ((pPeaks = (RAMPREAL *) malloc((peaksCount+1) * 2 * sizeof(RAMPREAL) + 1)) == NULL)
      {
         printf("Cannot allocate memory\n");
         return NULL;
      }
      
      // And byte order correction
      byteOrderOK = (isLittleEndian==weAreLittleEndian);
      if (32==precision) { // floats
         if (byteOrderOK) {
            for (n = 0; n < (2 * peaksCount); n++) {
               pPeaks[n] = (RAMPREAL) ((float *) pDecoded)[n];
            } 
         } else {
            for (n = 0; n < (2 * peaksCount); n++) {
               uint32_t tmp = swapbytes(((uint32_t *) pDecoded)[n]);
               pPeaks[n] = (RAMPREAL) * (float *)(&tmp);
            } 
         }
      } else { // doubles
         if (byteOrderOK) {
            for (n = 0; n < (2 * peaksCount); n++) {
               pPeaks[n] = (RAMPREAL)((double *) pDecoded)[n];
            }
         } else {
            for (n = 0; n < (2 * peaksCount); n++) {
               uint64_t tmp = swapbytes64((uint64_t) ((uint64_t *) pDecoded)[n]);
               pPeaks[n] = (RAMPREAL) * (double *)(&tmp);
            }
         }
      }
      free(pDecoded);
      pPeaks[n] = -1;
   }

   return (pPeaks); // caller must free this pointer
}


/*
 * read just the info available in the msRun element
 */
void readMSRun(RAMPFILE *pFI,
                   struct RunHeaderStruct *runHeader)
{
   char stringBuf[SIZE_BUF];
   ramp_fseek(pFI, 0 , SEEK_SET); // rewind
   ramp_fgets(stringBuf, SIZE_BUF, pFI);

   while((!strstr( stringBuf , pFI->bIsMzData?"<mzData":"<msRun" )) && !ramp_feof(pFI))  /* this should not be needed if index offset points to correct location */
   {
      ramp_fgets(stringBuf, SIZE_BUF, pFI);
   }
   while(!ramp_feof(pFI))
   {
      const char *cp;
      if (NULL != (cp=strstr( stringBuf , pFI->bIsMzData?"spectrumList count":"scanCount" ))) {
         cp = strchr(cp,'\"');
         runHeader->scanCount = atoi(cp+1);
      }
      if (NULL != (cp=strstr( stringBuf , "startTime" ))) {
         cp = strchr(cp,'\"');
         runHeader->dStartTime = rampReadTime(pFI,cp+1);
      } 
      if (NULL != (cp=strstr( stringBuf , "endTime" ))) {
         cp = strchr(cp,'\"');
         runHeader->dEndTime = rampReadTime(pFI,cp+1);
      } 
      if (NULL != (cp=strstr( stringBuf , pFI->bIsMzData?"<spectrumDesc":"<scan" ))) {
         break; /* we're into data territory now */
      } 
      ramp_fgets(stringBuf, SIZE_BUF, pFI);
   }

}

/*
 * walk through each scan to find overall lowMZ, highMZ
 * sets overall start and end times also
 */
void readRunHeader(RAMPFILE *pFI,
                   ramp_fileoffset_t *pScanIndex,
                   struct RunHeaderStruct *runHeader,
                   int iLastScan)
{

   int i;
   struct ScanHeaderStruct scanHeader;
   double startMz = 0.0;
   double endMz = 0.0;

   readHeader(pFI, pScanIndex[1], &scanHeader);

   /*
    * initialize values to first scan
    */
   runHeader->lowMZ = scanHeader.lowMZ;
   runHeader->highMZ = scanHeader.highMZ;
   runHeader->dStartTime = scanHeader.retentionTime;
   runHeader->startMZ = readStartMz( pFI , pScanIndex[1] );
   runHeader->endMZ = readEndMz( pFI , pScanIndex[1] );
   for (i = 2; i <= iLastScan; i++)
   {
      readHeader(pFI, pScanIndex[i], &scanHeader);

      if (scanHeader.lowMZ < runHeader->lowMZ)
         runHeader->lowMZ = scanHeader.lowMZ;
      if (scanHeader.highMZ > runHeader->highMZ)
         runHeader->highMZ = scanHeader.highMZ;
      if( (startMz = readStartMz( pFI , pScanIndex[i] )) < runHeader->startMZ )
	runHeader->startMZ = startMz;
      if( (endMz = readEndMz( pFI , pScanIndex[i] )) > runHeader->endMZ )
	runHeader->endMZ = endMz;   
   }

   runHeader->dEndTime = scanHeader.retentionTime;
}



int setTagValue(const char* text,
   char* storage,
   int maxlen,
   const char* lead,
   const char* tail)
{
  char* result = NULL;
  char* term = NULL;
  int len = maxlen - 1;

  result = strstr(text, lead);
  if(result != NULL)
  {
    term = strstr(result + strlen(lead), tail);
    if(term != NULL)
    {
      if((int)(strlen(result) - strlen(term) - strlen(lead)) < len)
	len = strlen(result) - strlen(term) - strlen(lead);

      strncpy(storage, result + strlen(lead), len);
      storage[len] = 0;
      return 1;
    } // if term
  }
  return 0;
}


InstrumentStruct* getInstrumentStruct(RAMPFILE *pFI)
{
  InstrumentStruct* output = NULL;
  char* result = NULL;
  int found[] = {0, 0, 0, 0, 0};
  char stringBuf[SIZE_BUF];

   if ((output = (InstrumentStruct *) calloc(1,sizeof(InstrumentStruct))) == NULL)
   {
      printf("Cannot allocate memory\n");
      return NULL;
   } else {
      const char *cpUnknown="UNKNOWN";
      strncpy(output->analyzer,cpUnknown,sizeof(output->analyzer));
      strncpy(output->detector,cpUnknown,sizeof(output->analyzer));
      strncpy(output->ionisation,cpUnknown,sizeof(output->detector));
      strncpy(output->manufacturer,cpUnknown,sizeof(output->ionisation));
      strncpy(output->model,cpUnknown,sizeof(output->manufacturer));
   }

   ramp_fgets(stringBuf, SIZE_BUF, pFI);

   if (pFI->bIsMzData) {
   } else {
      int isAncient=0;
      while( !strstr( stringBuf , "<msInstrument" ) && 
            !(isAncient=(NULL!=strstr( stringBuf , "<instrument" ))) && 
            !strstr(stringBuf, "<dataProcessing") && 
            !ramp_feof(pFI))  /* this should not be needed if index offset points to correct location */
      {
         ramp_fgets(stringBuf, SIZE_BUF, pFI);
      }
      while(! strstr(stringBuf, isAncient?"</instrument":"</msInstrument") &&  ! strstr(stringBuf, "</dataProcessing") && !ramp_feof(pFI))
      {
        if(! found[0])
        {
          result = strstr(stringBuf, isAncient?"manufacturer=":"<msManufacturer");
          if(result != NULL && setTagValue(result, output->manufacturer, INSTRUMENT_LENGTH, isAncient?"manufacturer=\"":"value=\"", "\""))
	      found[0] = 1;
        }
        if(! found[1])
        {
           result = strstr(stringBuf, isAncient?"model=":"<msModel");
          if(result != NULL && setTagValue(result, output->model, INSTRUMENT_LENGTH, isAncient?"model=\"":"value=\"", "\""))
	      found[1] = 1;
        }
        if(! found[2])
        {
          result = strstr(stringBuf, isAncient?"ionisation=":"<msIonisation");
          if(result != NULL && setTagValue(result, output->ionisation, INSTRUMENT_LENGTH, isAncient?"ionisation=\"":"value=\"", "\""))
	      found[2] = 1;
        }
        if(! found[3])
        {
           result = strstr(stringBuf, isAncient?"msType=":"<msMassAnalyzer");
          if(result != NULL && setTagValue(result, output->analyzer, INSTRUMENT_LENGTH, isAncient?"msType=\"":"value=\"", "\""))
	      found[3] = 1;
        }
        if(! found[4])
        {
          result = strstr(stringBuf, "<msDetector");
          if(result != NULL && setTagValue(result, output->detector, INSTRUMENT_LENGTH, "value=\"", "\""))
	      found[4] = 1;
        }
        ramp_fgets(stringBuf, SIZE_BUF, pFI);

      } // while
   }

   if(found[0] || found[1] || found[2] || found[3] || found[4])
     return output;

   return NULL; // no data
}

// exercise at least some of the ramp interface - return non-0 on failure
int rampSelfTest(char *filename) { // if filename is non-null we'll exercise reader with it
#define N_TEST_NAME 5
   int result = 0; // assume success
   char buf[256];
   char buf2[256];
   int i;

   char *testname[N_TEST_NAME] = 
   {"foo.bar","foo.mzxml","foo.mzdata","foo.mzXML","foo.mzData"};

   // locate the .mzData or .mzXML extension in the buffer
   // return pointer to extension, or NULL if not found
   for (i=N_TEST_NAME;i--;) {
      result |= (!i) != !rampValidFileType(testname[i]); // 0th one in not a valid file type
   }

   // trim a filename of its .mzData or .mzXML extension
   // return trimmed buffer, or null if no proper .ext found
   for (i=N_TEST_NAME;i--;) {
      strncpy(buf,testname[i],sizeof(buf));
      result |= ((!i) != !rampTrimBaseName(buf));
      if (i) {
         result |= (strcmp(buf,"foo")!=0);
      }
   }

   if (filename && rampValidFileType(filename)) {
      // construct a filename in buf from a basename, adding .mzXML or .mzData as exists
      // returns buf, or NULL if neither .mzXML or .mzData file exists
      char *name;
      strncpy(buf,filename,sizeof(buf));
      rampTrimBaseName(buf);
      name = rampConstructInputFileName(buf,sizeof(buf),buf); // basename is in buf
      result |= (name==NULL);
      strncpy(buf,filename,sizeof(buf));
      rampTrimBaseName(buf);
      name = rampConstructInputFileName(buf2,sizeof(buf2),buf); // different buf, basename
      result |= (name==NULL);
   }
   return result;
}

