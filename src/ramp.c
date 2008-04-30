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
*    This program is free software; you can redistribute it and/or modify  *
*    it under the terms of the GNU Library or "Lesser" General Public      *
*    License (LGPL) as published by the Free Software Foundation;          *
*    either version 2 of the License, or (at your option) any later        *
*    version.                                                              *
***************************************************************************/

#include "ramp.h"

#include "base64.h"

#include "time.h"

#include "zlib.h"

#define SIZE_BUF 512

#ifdef _MSC_VER
#define file_sep '\\'
#else
#define file_sep '/'
#endif

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

// 
// do casts through unions to avoid running afoul of gcc strict aliasing
//
typedef union {
   uint32_t u32;
   float flt;
} U32;

typedef union {
   uint64_t u64;
   double dbl;
} U64;

/****************************************************************
 * Utility functions					*
 ***************************************************************/

static char *findquot(const char *cp) { /* " and ' are both valid attribute delimiters */
   char *result = strchr(cp,'\"');
   if (!result) {
      result = strchr(cp,'\'');
   }
   return result;
}

static int isquot(const char c) { /* " and ' are both valid attribute delimiters */
      return ('\"'==c)||('\''==c);
}

static int setTagValue(const char* text,
      char* storage,
      int maxlen,
      const char* lead);

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
      p = findquot(p);
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
      buf[nread+nread_now] = 0;
      if (!nread_now) {
         return nread?buf:NULL;
      }
      newline = strchr(buf+nread,'\n');
      if (newline) {
         *(newline+1) = 0; // real fgets includes the newline
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
   char indexOffsetTemp[SIZE_BUF+1], buf;

   if (pFI->bIsMzData) {
      return -1; // no index in mzData
   }

   for (indexOffsetOffset = -120;  indexOffsetOffset++ < 0 ;)
   {
      char seekbuf[SIZE_BUF+1];
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
   int  n, nread;
   int  reallocSize = 8000;    /* initial # of scan indexes to expect */
   char *beginScanOffset;
   
   int newN;
   char *beginOffsetId;
   
   ramp_fileoffset_t *pScanIndex=NULL;
   int retryLoop;
   char* s;
   
   for (retryLoop = 2;retryLoop--;) {
     n = 1; // ramp is one based
     *iLastScan = 0;
     free(pScanIndex);
      if ((indexOffset < 0) || (retryLoop==0)) { // derive the index by inspection
         // no index found, derive it
        
        // HENRY - look for <scan num" instead of just <scan - so that we can more easily access the actual scan number. same for mzData 
        const char *scantag = pFI->bIsMzData?"<spectrum id=\"":"<scan num=\"";
         int taglen = (int)strlen(scantag);
         ramp_fileoffset_t index = 0;
        // HENRY - in this new implementation, n should start at zero
        n = 0;
         pScanIndex = (ramp_fileoffset_t *)malloc( sizeof(ramp_fileoffset_t)*reallocSize); // allocate space for the scan index info
         if (!pScanIndex) {
            printf("Cannot allocate memory\n");
            return NULL;
         }
         ramp_fseek(pFI,0,SEEK_SET);
         buf[sizeof(buf)-1] = 0;
         while ((nread = (int)ramp_fread(buf,sizeof(buf)-1,pFI))>taglen) {
            char *find;
            int truncated = 0;
            char *look=buf;
            buf[nread] = 0;
            while (NULL != (find = strstr(look,scantag))) {
              int k,newN; 
              // HENRY - needs to read ahead a few chars to make sure the scan num is complete in this buf
              char *scanNumStr = find + taglen; // pointing to the first digit of the scan num
               while (++scanNumStr < buf + sizeof(buf) - 1 && *scanNumStr != '\"'); // increment until it hits the end quote or the end of buffer 
               if (scanNumStr >= buf + sizeof(buf) - 1) { 
                  // hitting the end of buffer, let's not read this scan; remember the length of the truncated piece
                  truncated = scanNumStr - find;
                 break;
               }
               
               // HENRY - reset scanNumStr to start of scan num
               scanNumStr = find + taglen;
               
               // HENRY - atoi will read until the end quote
               newN = atoi(scanNumStr);
 
               //              printf("newN = %d, offset = %lld\n", newN, index + (find - buf));
               // HENRY - realloc needs to make sure newN has a spot in pScanIndex
               if (reallocSize <= newN) {
                 reallocSize = newN + 500; 
                 pScanIndex = (ramp_fileoffset_t *)realloc(pScanIndex, sizeof(ramp_fileoffset_t)*reallocSize);
                 if (!pScanIndex) {
                   printf("Cannot allocate memory\n");
                   return NULL;
                 }
               }               
               
               // HENRY - sets all the skipped scans to offset -1 (here you see why I set n = 0 to begin, rather than n = 1)
               for (k = n + 1; k < newN; k++) {
                 pScanIndex[k] = -1;
               }
               
               // HENRY - puts the offset at pScanIndex[newN]
               pScanIndex[newN] = index+(find-buf); // ramp is 1-based
               n = newN;
               (*iLastScan) = newN;
               
               // HENRY - we can start looking from the end quote of the last scan number.
               look = scanNumStr;
               
               // HENRY - reallocation needs to happen earlier, before we set pScanIndex[newN], in case newN is already past the old alloc size
               /*
               if (reallocSize<=n) {
                  reallocSize+=500;
                  pScanIndex = (ramp_fileoffset_t *)realloc(pScanIndex, sizeof(ramp_fileoffset_t)*reallocSize); // allocate space for the scan index info
                  if (!pScanIndex) {
                     printf("Cannot allocate memory\n");
                     return NULL;
                  }
               }
               */
            }
            nread = strlen(look)+(look-buf);
            if (*look && strchr(scantag,buf[nread-1]) && !ramp_feof(pFI)) { // check last char of buffer
               // possible that next scantag overhangs end of buffer
               ramp_fseek(pFI,-taglen,SEEK_CUR); // so next get includes it
            
            // HENRY - if the scan number is truncated, we go back a few chars so that the next get will include it
            } else if (truncated != 0 && !ramp_feof(pFI)) {
               ramp_fseek(pFI, -truncated, SEEK_CUR);
            } 
            index = ramp_ftell(pFI);
         }
         break; // no need to retry
      } else {  // read the index out of the file (then check it!)
         struct ScanHeaderStruct scanHeader; // for test of index integrity
         int indexOK=1; // until we show otherwise

         // HENRY -- reset n to zero. Note that it should be zero here, not one -- as n points to the previous record in 
         // my nomenclature (and newN to the newly read record).
         n = 0;
         
         if ((pScanIndex = (ramp_fileoffset_t *) malloc(reallocSize * sizeof(ramp_fileoffset_t))) == NULL) {
            printf("Cannot allocate memory\n");
            return NULL;
         }
         
         ramp_fseek(pFI, indexOffset, SEEK_SET);
         
         s = ramp_fgets(buf, SIZE_BUF, pFI);
         while( s!=NULL && !strstr( buf , "<offset" ) ) {
            s = ramp_fgets(buf, SIZE_BUF, pFI);
         }
         
         if (s == NULL)
            break;   // end of file reached.
         
         while (!strstr(buf, "/index")) {
            int k;
            // HENRY -- also reads the "id" field, which is the scan num
            if ((beginOffsetId = (char *)(strstr(buf, "id=\""))) == NULL) {
               ramp_fgets(buf, SIZE_BUF, pFI);
               continue;
            }
            beginOffsetId += 4;
            
            newN = atol(beginOffsetId);
            
            // HENRY -- check if the new id is past the max size of the pScanIndex array
            // Note that it should be reallocSize - 1, because the very last record is set to offset=-1 
            // (see below)! In case newN is the very last record, we need to prepare the space for the offset=-1 thingy.
            if (newN >= reallocSize - 1) {
              ramp_fileoffset_t *pTmp;
              
              // HENRY -- we don't know how much newN is bigger than the old realloc size. In case it is more than 500 bigger,
              // then the old way of always reallocating for 500 more will break. Instead we jump to newN + 500.
              reallocSize = newN + 500;
              pTmp = (ramp_fileoffset_t*)realloc(pScanIndex, reallocSize * sizeof(ramp_fileoffset_t));
              if (pTmp == NULL) { 
                printf("Cannot allocate memory\n");
                return NULL;
              } else {
                pScanIndex=pTmp;
              }
            }
            
            // HENRY -- any scan number skipped between the last record and the new record will be assumed "missing"
            // the offset will be set to -1 for these scan numbers
            for (k = n + 1; k < newN; k++) {
              pScanIndex[k] = -1;
            }
            // HENRY -- this new record becomes the current one, and iLastScan gets the scan number of this new record
            n = newN;
            (*iLastScan) = n;
            
            // HENRY -- using merely the ">" as the beginning of the offset is somewhat scary, but I'm not changing it now.
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

            // HENRY -- I have moved the following realloc piece earlier. The reason is:
            // In the old way, the scan numbers are assumed to be consecutive, so one can expect the next scan number
            // to be 1 + the current one. In this case, you only have to make sure space is allocated for one more record.
            // In the new way, the scan numbers are not consecutive, so how many 
            // more spaces we need here is unpredictable. (If the next scan number is 1000 + this current one, then we could be in 
            // trouble.) It then makes sense to realloc AFTER we read the next scan number, which is what I am doing.
            
            //            printf ("(%d, %ld) ", n, pScanIndex[n]);
            //            n++;
//            (*iLastScan)++;
/*            
            if (n >= reallocSize) {
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
  */          
            ramp_fgets(buf, SIZE_BUF, pFI);
         }
         
         // HENRY -- We have no idea whether scan number 1, n/2 or n-1 is a missing scan or not. So we cannot just blindly test them.
         // Instead, we start from 1, n/2 and n to find a valid offset to test. (we can test n because n still points to the
         // last scan number (never n++ in this implementation). 
         
         if (n > 0) {
            // OK, now test that to see if it's a good index or not
            int testIndex = 1;   
            // HENRY -- iteratively finds the next valid offset
            while (testIndex <= n && pScanIndex[testIndex] <= 0) testIndex++;
            if (testIndex <= n) {
              readHeader(pFI,pScanIndex[testIndex],&scanHeader);
              if (scanHeader.acquisitionNum == -1) { // first
                 indexOK = 0; // bogus index
                 free(pScanIndex);
                 pScanIndex = NULL;
              }
            } 
            // HENRY -- n>3 is better. if n=2 or 3, n/2 is 1, which we just tested. 
            if (indexOK && (n>3)) { // middle
               testIndex = n/2;
               while (testIndex <= n && pScanIndex[testIndex] <= 0) testIndex++;
               if (testIndex <= n) {
                 readHeader(pFI,pScanIndex[testIndex],&scanHeader);
                 if (scanHeader.acquisitionNum == -1) {
                    indexOK = 0; // bogus index
                    free(pScanIndex);
                    pScanIndex = NULL;
                 }
              }
            }
            
            if (indexOK && (n>1)) { // last
               testIndex = n;
               while (testIndex >= 1 && pScanIndex[testIndex] <= 0) testIndex--;
               if (testIndex >= 1) {
                 readHeader(pFI,pScanIndex[testIndex],&scanHeader);
                 if (scanHeader.acquisitionNum == -1) {
                   indexOK = 0; // bogus index
                   free(pScanIndex);
                   pScanIndex = NULL;
                 }
               }
            }
     //       HENRY - Uncomment following to activate index creation from scratch
     //       indexOK = 0;
         }
         
         if (indexOK) {
            break; // no retry
         }
      } // end if we claim to have an index
   } // end for retryloop
   
   // HENRY -- Here we set the n+1 record to offset=-1. (Note that unlike in the old implementation, I have never n++,
   // so n still is the scan number of the last record.
   pScanIndex[n + 1] = -1;

   return (pScanIndex);
}

// helper func for reading mzData
const char *findMzDataTagValue(const char *pStr, const char *tag) {
   const char *find = strstr(pStr,tag);
   if (find) {
      find = strstr(find+1,"value=");
      if (find) {
         find = findquot(find);
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
      while (!isquot(*++pStr)) {
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


//
// helper function to deal with mzdata with no newlines - breaks
// lines up at </...>
//
static char *ramp_nextTag(char *buf, int buflen, RAMPFILE *pFI) {
   char *result;
   result = ramp_fgets(buf,buflen,pFI);
   if (result && !strchr(buf,'\n')) { // no newline found
      char *closer = strstr(buf+1,"</");
      if (closer) {
         *closer = 0; // temp. nullterm
         ramp_fseek(pFI,(1+closer-buf)-buflen,SEEK_CUR); // reposition for next read
      }
   }
   return result;
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

   /*
    * initialize defaults
    */
   memset(scanHeader,0,sizeof(struct ScanHeaderStruct)); // mostly we want 0's
#define LOWMZ_UNINIT 1.111E6
   scanHeader->lowMZ =  LOWMZ_UNINIT;
   scanHeader->acquisitionNum = -1;
   scanHeader->seqNum = -1;

   // HENRY - missing scans due to dta2mzXML get offset of zero
   // missing scans without index entries get offset of -1
   // check for those cases and populate an empty scanHeader and return
   // if we fseek(-1) unexpected behavior will result!
   if (lScanIndex <= 0) { 
     return;
   }

   ramp_fseek(pFI, lScanIndex, SEEK_SET);



   if (pFI->bIsMzData) {
      int bHasPrecursor = 0;
      while (ramp_nextTag(stringBuf, SIZE_BUF, pFI))
      {
         const char *attrib;
         const char *pStr;
         const char *closeTag = strstr(stringBuf, "</spectrumSettings>");
         // find each attribute in stringBuf
         for (attrib=stringBuf-1;NULL!=(attrib=strchr(attrib+1,'='));) {
            if (closeTag && (closeTag < attrib)) {
               break; // into data territory now
            }
            if ((pStr = matchAttr(attrib, "spectrum id",11))) {
               sscanf(pStr, "%d", &(scanHeader->acquisitionNum));
            //} else if ((pStr = matchAttr(attrib, "basePeakMz",10)))  {
            //   sscanf(pStr, "%lf", &(scanHeader->basePeakMZ));      
            //} else if ((pStr = matchAttr(attrib, "totIonCurrent",13)))  {
            //   sscanf(pStr, "%lf", &(scanHeader->totIonCurrent));      
            //} else if ((pStr = matchAttr(attrib, "basePeakIntensity",17)))  {
            //   sscanf(pStr, "%lf", &(scanHeader->basePeakIntensity));      
            } else if ((pStr = matchAttr(attrib, "msLevel",7)))  {
               sscanf(pStr, "%d", &(scanHeader->msLevel));
            //} else if ((pStr = matchAttr(attrib, "length",6)))  { get this from array element
            //   sscanf(pStr, "%d", &(scanHeader->peaksCount));
            } else if ((pStr = findMzDataTagValue(attrib,"TimeInMinutes")))  {
               scanHeader->retentionTime = rampReadTime(pFI,stringBuf);
            } else if ((pStr = findMzDataTagValue(attrib,"TimeInSeconds")))  {
               scanHeader->retentionTime = rampReadTime(pFI,stringBuf);
            } else if ((pStr = matchAttr(attrib, "mzRangeStart",12)))  {
               sscanf(pStr, "%lf", &(scanHeader->lowMZ));
            } else if ((pStr = matchAttr(attrib, "mzRangeStop",11)))  {
               sscanf(pStr, "%lf", &(scanHeader->highMZ));
            } else if ((pStr = findMzDataTagValue(attrib, "ScanMode"))) { 
               if ((pStr2 = (char *) findquot(pStr))) {
                  memcpy(&(scanHeader->scanType), pStr, pStr2-pStr);
                  scanHeader->scanType[pStr2-pStr] = '\0';
               }
            }
         }
         if (closeTag) {
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
         } else if ((pStr = (char *) strstr(stringBuf, "</precursorList")))
         {
            bHasPrecursor = 0;
         }
         if (bHasPrecursor) { // in precursor section
            if (NULL!=(pStr2 = (char *) strstr(stringBuf, "spectrumRef="))) 
            {
               sscanf(pStr2 + 13, "%d", &(scanHeader->precursorScanNum));
            }
            if (NULL!=(pStr2 = findMzDataTagValue(stringBuf,"ChargeState"))) 
            {
               scanHeader->precursorCharge = atoi(pStr2);
            }
            //Paul : added support for Collision Energy 25-01-08
            if (NULL!=(pStr2 = findMzDataTagValue(stringBuf,"CollisionEnergy"))) 
            {
               scanHeader->collisionEnergy = atoi(pStr2);
            }
            
            if (NULL!=(pStr2 = findMzDataTagValue(stringBuf,"MassToChargeRatio"))) 
            {
               scanHeader->precursorMZ = atof(pStr2);
//Paul altered atoi to atof so we get the decimals :D 04.04.08
            }
            if (NULL!=(pStr2 = findMzDataTagValue(stringBuf,"mz"))) 
            {
               scanHeader->precursorMZ = atof(pStr2);
            }
	    // STN
            if (NULL!=(pStr2 = findMzDataTagValue(stringBuf,"Intensity"))) 
            {
               scanHeader->precursorIntensity = atof(pStr2); //PB: atoi -> atof 
            }
         }
         if (strstr(stringBuf, "</spectrumDesc>")) {
            break; // into data territory now
         }
         if (strstr(stringBuf, "</precursorList>")) {
            break; // into data territory now
         }
      } while (ramp_nextTag(stringBuf, SIZE_BUF, pFI));
      // now read peaks count
      do {
         if (strstr(stringBuf, "ArrayBinary>")) {
            do {
               char *cp=strstr(stringBuf,"length=");
               if (cp) {
                  scanHeader->peaksCount = atoi(cp+8);
                  break;
               }
            } while (ramp_nextTag(stringBuf, SIZE_BUF, pFI));
            break;
         }
      } while (ramp_nextTag(stringBuf, SIZE_BUF, pFI));
   } else { // mzXML
      while (ramp_fgets(stringBuf, SIZE_BUF, pFI))
      {
         const char *attrib;
         const char *pStr;
         // find each attribute in stringBuf
         for (attrib=stringBuf-1;NULL!=(attrib=strchr(attrib+1,'='));) {
            if ((pStr = matchAttr(attrib, "num",3))) {
               sscanf(pStr, "%d", &(scanHeader->acquisitionNum));
            } else if ((pStr = matchAttr(attrib, "basePeakMz",10)))  {
               sscanf(pStr, "%lf", &(scanHeader->basePeakMZ));      
            } else if ((pStr = matchAttr(attrib, "totIonCurrent",13)))  {
               sscanf(pStr, "%lf", &(scanHeader->totIonCurrent));      
            } else if ((pStr = matchAttr(attrib, "basePeakIntensity",17)))  {
               sscanf(pStr, "%lf", &(scanHeader->basePeakIntensity));      
            } else if ((pStr = matchAttr(attrib, "msLevel",7)))  {
               sscanf(pStr, "%d", &(scanHeader->msLevel));
            } else if ((pStr = matchAttr(attrib, "peaksCount",10)))  {
               sscanf(pStr, "%d", &(scanHeader->peaksCount));
            } else if ((pStr = matchAttr(attrib, "retentionTime",13)))  {
               scanHeader->retentionTime = rampReadTime(pFI,pStr);
            } else if ((pStr = matchAttr(attrib, "lowMz",5)))  {
               sscanf(pStr, "%lf", &(scanHeader->lowMZ));
            } else if ((pStr = matchAttr(attrib, "highMz",6)))  {
               sscanf(pStr, "%lf", &(scanHeader->highMZ));
            } else if ((scanHeader->lowMZ==LOWMZ_UNINIT) &&  
               ((pStr = matchAttr(attrib, "startMz",7)))) {
               sscanf(pStr, "%lf", &(scanHeader->lowMZ));  
            } else if ((!scanHeader->highMZ) &&
               ((pStr = matchAttr(attrib, "endMz",5)))) {
               sscanf(pStr, "%lf", &(scanHeader->highMZ));
            } else if ((pStr = matchAttr(attrib, "scanType", 8))) { 
               if ((pStr2 = (char *) findquot(pStr))) {
                  memcpy(&(scanHeader->scanType), pStr, sizeof(char)*((pStr2-pStr)));
                  scanHeader->scanType[pStr2-pStr] = '\0';
               }
            } else if ((pStr = matchAttr(attrib, "collisionEnergy", 15)))  {
                 sscanf(pStr, "%lf", &(scanHeader->collisionEnergy));
            }
         }
         
         /*
         * read precursor mass
         */
         if ((pStr = (char *) strstr(stringBuf, "<precursorMz ")))
         {
            if ((pStr2 = (char *) strstr(stringBuf, "precursorScanNum="))) 
            {
               sscanf(pStr2 + 18, "%d", &(scanHeader->precursorScanNum));
            }
            
            /*
            * Check for precursor charge.
            */
            if ((pStr2 = (char *) strstr(pStr, "precursorCharge="))) 
            {
               sscanf(pStr2 + 17, "%d", &(scanHeader->precursorCharge));
            }
            if ((pStr2 = (char *) strstr(pStr, "precursorIntensity="))) {
               sscanf(pStr2 + 20, "%lf", &(scanHeader->precursorIntensity));
            }
            
            /*
            * Find end of tag.
            */
            while (!(pStr = strchr(pStr, '>')))
            {      
               ramp_fgets(stringBuf, SIZE_BUF, pFI);
               pStr = stringBuf;

               if ((pStr2 = (char *) strstr(stringBuf, "precursorScanNum="))) 
                  sscanf(pStr2 + 18, "%d", &(scanHeader->precursorScanNum));
               if ((pStr2 = (char *) strstr(stringBuf, "precursorCharge="))) 
               {
                  sscanf(pStr2 + 17, "%d", &(scanHeader->precursorCharge));
               }
               if ((pStr2 = (char *) strstr(pStr, "precursorIntensity="))) {
                  sscanf(pStr2 + 20, "%lf", &(scanHeader->precursorIntensity));
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
            }
            
            sscanf(pStr, "%lf<", &(scanHeader->precursorMZ));
            //         printf("precursorMass = %lf\n", scanHeader->precursorMZ);
         }
         if (strstr(stringBuf, "<peaks")) {
            break; // into data territory now
         }
         if ((-1==scanHeader->acquisitionNum) &&
             ((strstr(stringBuf, "</dataProcessing>")||
               strstr(stringBuf, "</msInstrument>")))) {
            break; // ??? we're before scan 1 - indicates a broken index
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
   char stringBuf[SIZE_BUF+1];
   char szLevel[12];
   char *beginMsLevel, *endMsLevel;

   // HENRY - check if index is valid. the uninitialized value of ms level is probably zero.
   if (lScanIndex <= 0) {
     return (0);
   }

   ramp_fseek(pFI, lScanIndex, SEEK_SET);

   ramp_fgets(stringBuf, SIZE_BUF, pFI);

   while (!(beginMsLevel = (char *) strstr(stringBuf, "msLevel=")))
   {
      ramp_fgets(stringBuf, SIZE_BUF, pFI);
   }

   beginMsLevel += 9;           // We need to move the length of msLevel="
   endMsLevel = (char *) findquot(beginMsLevel);
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
  char stringBuf[SIZE_BUF+1];
  double startMz = 1.E6;
  char *pStr;
  const char *tag = pFI->bIsMzData?"mzRangeStart":"startMz";

  // HENRY -- again, check for invalid offset first. Is startMz = 1.E6 a good uninitialized value?
  if (lScanIndex <= 0) {
    return (startMz);
  }

  ramp_fseek(pFI, lScanIndex, SEEK_SET);
   
  while (ramp_fgets(stringBuf, SIZE_BUF, pFI))
  {
     if (strstr(stringBuf, pFI->bIsMzData?"</spectrumDesc":"<peaks"))
        break; // ran to end
      if ((pStr = strstr(stringBuf, tag))){
        sscanf(pStr + strlen(tag)+2, "%lf", &startMz);
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
  char stringBuf[SIZE_BUF+1];
  double endMz = 0.0;
  char *pStr;
  const char *tag = pFI->bIsMzData?"mzRangeStop":"endMz";

  // HENRY -- again, check for invalid offset first. is endMz = 0 a good uninitialized value?
  if (lScanIndex <= 0) {
    return (endMz);
  }

  ramp_fseek(pFI, lScanIndex, SEEK_SET);
   
  while (ramp_fgets(stringBuf, SIZE_BUF, pFI))
  {
     if (strstr(stringBuf, pFI->bIsMzData?"</spectrumDesc":"<peaks"))
        break; // ran to end
      if ((pStr = strstr(stringBuf, tag))){
        sscanf(pStr + strlen(tag)+2, "%lf", &endMz);
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
   return rampConstructInputPath(inbuf, inbuflen, "", basename_in);
}

char *rampConstructInputPath(char *inbuf,int inbuflen,const char *dir,const char *basename_in) {
   char *result = NULL;
   int i;
   char *basename = (char *) basename_in;
   char *tmpbuf = (char *)malloc(strlen(dir) + strlen(basename) + 20);
   char *append;
   if (dir != NULL && *dir != '\0')
   {
       // If directory for mzXML was supplied, strip off directory part.
       char *basename_sep = strrchr(basename, '/');
       if (basename_sep != NULL)
           basename = basename_sep + 1;
       basename_sep = strrchr(basename, '\\');
       if (basename_sep != NULL)
           basename = basename_sep + 1;
   }

   if (basename_in==inbuf) { // same pointer
      char *basename_buff = (char *)malloc(inbuflen);
      strncpy(basename_buff,basename,inbuflen);
      basename = basename_buff;
   }

   *tmpbuf= 0;
   if (dir != NULL && *dir != '\0')
   {
       int len_dir = strlen(dir);
   strcpy(tmpbuf, dir);
       if (tmpbuf[len_dir - 1] != file_sep)
       {
           tmpbuf[len_dir] = file_sep;
           tmpbuf[len_dir+1] = 0;
       }
   }
   strcat(tmpbuf, basename);

   append = tmpbuf + strlen(tmpbuf);

   for (i=0;i<N_EXT_TYPES;i++) {
      FILE *test;
      strcpy(append,data_ext[i]);
      test = fopen(tmpbuf,"r");
      if (test != NULL) {
         if (result) { // conflict! both mzXML and mzData are present
            if (strcasecmp(tmpbuf,result)) { // win32 isn't case sensitive
               printf("found both %s and %s, using %s\n",
                  tmpbuf,result,result);
            }
         } else { // found file, copy the constructed name
            result = strdup(tmpbuf);
         }
         fclose(test);
      }
   }
   if (!result) { // failed - caller can complain about lack of .mzXML
      strcpy(append,data_ext[0]);
      result = strdup(tmpbuf);
   }
   if (basename_in==inbuf) { // same pointer
      free(basename);
   }
   free(tmpbuf);

   if ((int) strlen(result) < inbuflen) {
      strcpy(inbuf, result);
      free(result);
      result = inbuf;
   } else {
      printf("buffer too small for file %s\n",
         result);
      free(result);
      result = NULL;
   }
   return result;
}

// return NULL if fname is not of extension type we handle,
// otherwise return pointer to .ext
char *rampValidFileType(const char *fname) {
   char *result;
   int i;
   for (i = N_EXT_TYPES;i--;) {
      result = strrchr(fname,'.');
      if (result && !strcasecmp(result,data_ext[i])) {
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
   char *stringBuf=(char *)malloc(SIZE_BUF+1);
   char *beginPeaksCount, *peaks;
   int result = 0;
   const char *tag = pFI->bIsMzData?"length=":"peaksCount=";
   ramp_fileoffset_t in_lScanIndex = lScanIndex;

   // HENRY -- check invalid offset. is 0 a good uninitialized value for peakCount?
   if (lScanIndex <= 0) {
     return (0);
   }

   ramp_fseek(pFI, lScanIndex, SEEK_SET);

   // Get the num of peaks in the scan and allocate the space we need
   ramp_nextTag(stringBuf, SIZE_BUF, pFI);
   while (!(beginPeaksCount = (char *) strstr(stringBuf, tag)))
   {
      lScanIndex = ramp_ftell(pFI);
      ramp_nextTag(stringBuf, SIZE_BUF, pFI);
   }

   // We need to move forward the length of the tag
   beginPeaksCount += (strlen(tag)+1);
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
   free(stringBuf);
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
   RAMPREAL *pPeaksDeRuled = NULL;
   
   int  endtest = 1;
   int  weAreLittleEndian = *((char *)&endtest);

   char *pData = NULL;
   char *pBeginData;
   char *pDecoded = NULL;

   char buf[1000];
   buf[sizeof(buf)-1] = 0;

   // HENRY - check invalid offset... is returning NULL here okay? I think it should be because
   // NULL is also returned later in this function when there is no peak found.
   if (lScanIndex <= 0) {
     return (NULL);
   }


   if (pFI->bIsMzData) {
      // intensity and mz are written in two different arrays
      int bGotInten = 0;
      int bGotMZ = 0;
      ramp_fseek(pFI,lScanIndex,SEEK_SET);
      while ((!(bGotInten && bGotMZ)) &&
           ramp_nextTag(buf,sizeof(buf)-1,pFI)) {
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
            const char *datastart;
            // now determine peaks count, precision
            while (!(datastart= (char *) strstr(buf, "<data")))
            {
               ramp_nextTag(buf, sizeof(buf)-1, pFI);
            }
            
            // find precision="xx"
            if( !(pBeginData = strstr( buf , "precision=" )))
            {
               precision = 32; // default value
            } else { // we found declaration
               precision = atoi(findquot(pBeginData)+1);
            }
            
            // find length="xx"
            if((!peaksCount) && (pBeginData = strstr( buf , "length=" )))
            {
               peaksCount = atoi(findquot(pBeginData)+1);
            }

            if (peaksCount <= 0)
            { // No peaks in this scan!!
               return NULL;
            }

            // find endian="xx"
            if((pBeginData = strstr( buf , "endian=" )))
            {
               isLittleEndian = !strncmp("little",findquot(pBeginData)+1,6);
            }
            
            // find close of <data>
            while( !(pBeginData = strstr( datastart , ">" )))
            {
               ramp_nextTag(buf, sizeof(buf)-1 , pFI);
               datastart = buf;
            }
            pBeginData++;	// skip the >
            
            // base64 has 4:3 bloat, precision/8 bytes per value
            bytes = (peaksCount*(precision/8));
            // for every 3 bytes base64 emits 4 characters - 1, 2 or 3 byte input emits 4 bytes
            triplets = (bytes/3)+((bytes%3)!=0);
            peaksLen = (4*triplets)+1; // read the "<" from </data> too, to confirm lack of whitespace
            
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

            // whitespace may be present in base64 char stream
            while (pData[peaksLen-1]!='<') {
               char *cp;
               partial = 0;
               // didn't read all the peak info - must be whitespace
               for (cp=pData;*cp;) {
                  if (strchr("\t\n\r ",*cp)) {
                     memmove(cp,cp+1,peaksLen-(partial+cp-pData));
                     partial++;
                  } else {
                     cp++;
                  }
               }
               if (!ramp_fread(pData+peaksLen-partial,partial, pFI)) {
                  break;
               }
            }
            pData[peaksLen-1] = 0; // pure base64 now
                        
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
                  U32 tmp;
                  for (n = 0; n < peaksCount; n++) {
                     tmp.u32 = swapbytes( *u++ );
                     pPeaks[isInten+(2*n)] = (RAMPREAL) tmp.flt;
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
                  U64 tmp;
                  for (n = 0; n < peaksCount; n++) {
                     tmp.u64 = swapbytes64( *u++ );
                     pPeaks[isInten+(2*n)] = (RAMPREAL) tmp.dbl;
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
       Byte *pUncompr;
       int isCompressed = 0;
      int partial,bytes,triplets;
      int isLittleEndian = 0; // default is network byte order (Big endian)
      int byteOrderOK;
      int       compressedLen = 0;
      int       decodedSize;
      char      *pToBeCorrected;
      e_contentType contType = mzInt; // default to m/z-int

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

          // TODO ALL OF THE FOLLOWING CHECKS ASSUME THAT THE NAME AND THE VALUE OF THE
          // ATTRIBUTE ARE PRESENT AT THE SAME TIME IN THE BUFFER.
          // ADD A CHECK FOR THAT!
      while( 1 )
      { // Untill the end of the peaks element
          if( (pBeginData = strstr( buf , "precision=" )))
          { // read the precision attribute
              precision = atoi(strchr(pBeginData,'\"')+1);
          }
          if( (pBeginData = strstr( buf , "contentType=" )))
          { // read the contentType attribute
                  // we are only supporting m/z-int for the moment > return if it is something else
                  // TODO add support for the other content types
              if( (pBeginData = strstr( buf , "m/z-int" )))
              {
                  contType = mzInt;
              }
              else if( (pBeginData = strstr( buf , "m/z ruler" )))
              {
                  contType = mzRuler;
              }
              else
              {
                  char* pEndAttrValue;
                  pEndAttrValue = strchr( pBeginData + strlen( "contentType=\"") + 1 , '\"' );
                  pEndAttrValue  = '\0';
                  fprintf(stderr, "%s Unsupported content type\n" , pBeginData ); 
                  return NULL;
              }
          }
          if( (pBeginData = strstr( buf , "compressionType=" )))
          { // read the compressionType attribute.
              if( (pBeginData = strstr( buf , "zlib" )))
              {
                  isCompressed = 1;
              }
              else if( (pBeginData = strstr( buf , "none" )))
              {
                  isCompressed = 0;
              }
              else
              {
                  char* pEndAttrValue;
                  pEndAttrValue = strchr( pBeginData + strlen( "compressionType=\"") + 1 , '\"' );
                  pEndAttrValue = '\0';
                  fprintf(stderr, "%s Unsupported compression type\n" , pBeginData ); 
                  return NULL;
              }
          }
          if( (pBeginData = strstr( buf , "compressedLen=\"")))
          {
              compressedLen = atoi( pBeginData + strlen( "compressedLen=\"" ) );
          }
          if( !(pBeginData = strstr( buf , ">" )))
          { // There is more to read
              ramp_fgets(buf, sizeof(buf)-1 , pFI);
              getIsLittleEndian(buf,&isLittleEndian);
          }
          else
          {
              pBeginData++;	// skip the >
              break;
          }
      }
      if( !precision )
      { // precision attribute was not defined assume 32 by default
          precision = 32;
      }

      if( isCompressed )
      {
          bytes = compressedLen;
      }
      else
      {
          bytes = (peaksCount*(precision/4));
      }
      
      // base64 has 4:3 bloat, precision/8 bytes per value, 2 values per peak
      // for every 3 bytes base64 emits 4 characters - 1, 2 or 3 byte input emits 4 bytes
      triplets = (bytes/3)+((bytes%3)!=0);
      peaksLen = (4*triplets)+1; // read the "<" from </data> too, to confirm lack of whitespace
      
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
      // whitespace may be present in base64 char stream
      while (pData[peaksLen-1]!='<') {
         char *cp;
         partial = 0;
         // didn't read all the peak info - must be whitespace
         for (cp=pData;*cp;) {
            if (strchr("\t\n\r ",*cp)) {
               memmove(cp,cp+1,peaksLen+1-(partial+cp-pData));
               partial++;
            } else {
               cp++;
            }
         }
         if (!ramp_fread(pData+peaksLen-partial,partial, pFI)) {
            break;
         }
      }
      if( isCompressed )
      {
          decodedSize = compressedLen + 1;
      }
      else
      {
          // 2 values per peak, precision/8 bytes per value
          decodedSize = peaksCount * (precision/4) + 1;
      }
      pData[peaksLen-1] = 0; // pure base64 now
      
      if ((pDecoded = (char *) malloc( decodedSize )) == NULL)
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

          //Zlib decompression 
      if( isCompressed )
      {
          int err;
//        printf("Decompressing data\n");
          uLong uncomprLen = (peaksCount * precision/4 + 1);
			
          pUncompr = (Byte*)calloc((uInt) uncomprLen , 1);
			
          err = uncompress( pUncompr , &uncomprLen , (const Bytef*)pDecoded , decodedSize );
          free( pDecoded );
          pToBeCorrected = (char *)pUncompr;
      }
      else
      {
          pToBeCorrected = pDecoded;
      }
      
      // And byte order correction
      byteOrderOK = (isLittleEndian==weAreLittleEndian);
      if (32==precision) { // floats
         if (byteOrderOK) {
            for (n = 0; n < (2 * peaksCount); n++) {
               pPeaks[n] = (RAMPREAL) ((float *) pToBeCorrected)[n];
            } 
         } else {
            U32 tmp;
            for (n = 0; n < (2 * peaksCount); n++) {
               tmp.u32 = swapbytes(((uint32_t *) pToBeCorrected)[n]);
               pPeaks[n] = (RAMPREAL) tmp.flt;
            } 
         }
      } else { // doubles
         if (byteOrderOK) {
            for (n = 0; n < (2 * peaksCount); n++) {
               pPeaks[n] = (RAMPREAL)((double *) pToBeCorrected)[n];
            }
         } else {
            U64 tmp;
            for (n = 0; n < (2 * peaksCount); n++) {
               tmp.u64 = swapbytes64((uint64_t) ((uint64_t *) pToBeCorrected)[n]);
               pPeaks[n] = (RAMPREAL) tmp.dbl;
            }
         }
      }

      if( contType == mzRuler )
      { // Convert back from m/z ruler contentType into m/z - int pairs
		RAMPREAL lastMass;
		RAMPREAL  deltaMass;
		  int multiplier;
          int j = 0;

          if ((pPeaksDeRuled = (RAMPREAL *) malloc((peaksCount+1) * 2 * sizeof(RAMPREAL) + 1)) == NULL)
          {
              printf("Cannot allocate memory\n");
              return NULL;
          }
         
          for (n = 0; n < (2 * peaksCount); )
          {
// printf("%f\n" , pPeaks[j] );
              if( (int) pPeaks[j] == -1 )
              { // Change in delta m/z
				++j;
				lastMass = (RAMPREAL) pPeaks[j++];
                deltaMass = pPeaks[j++];
				multiplier = 0;
//printf("%f %f\n" , lastMass , deltaMass );
              }
   		    pPeaksDeRuled[n++] = lastMass + (RAMPREAL) multiplier * deltaMass;
			++multiplier;
            pPeaksDeRuled[n++] = pPeaks[j++];
          }
          
          free(pToBeCorrected);
          pPeaksDeRuled[n] = -1;

          free( pPeaks );
          return (pPeaksDeRuled); // caller must free this pointer
      }
      
      free(pToBeCorrected);
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
   char stringBuf[SIZE_BUF+1];
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
         cp = findquot(cp);
         runHeader->scanCount = atoi(cp+1);
      }
      if (NULL != (cp=strstr( stringBuf , "startTime" ))) {
         cp = findquot(cp);
         runHeader->dStartTime = rampReadTime(pFI,cp+1);
      } 
      if (NULL != (cp=strstr( stringBuf , "endTime" ))) {
         cp = findquot(cp);
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
   int firstScan = 1;

   // HENRY -- initialize runHeader to some uninitialized values in case of failure
   runHeader->lowMZ = 0;
   runHeader->highMZ = 0;
   runHeader->dStartTime = 0;
   runHeader->startMZ = 1.E6;
   runHeader->endMZ = 0;
   
   // HENRY -- skipping over all the "missing scans"
   while (firstScan <= iLastScan && pScanIndex[firstScan] <= 0) { 
     firstScan++;
   }
   if (firstScan > iLastScan) {
     // HENRY -- this means there is no scan! do we need to initialize runHeader to something so that the caller
     // can check for this?
     return;
   }
   
   readHeader(pFI, pScanIndex[firstScan], &scanHeader);

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
     // HENRY -- skipping over all the missing scans
     if (pScanIndex[i] <= 0) {
       continue;
     }
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
   const char* lead)
{
  char* result = NULL;
  char* term = NULL;
  int len = maxlen - 1;
  int leadlen = strlen(lead)+1; // include the opening quote

  result = strstr(text, lead);
  if(result != NULL)
  {
    char tail = *(result+leadlen-1); // nab the quote char (is it single or double quote?)
    term = strchr(result + leadlen, tail);
    if(term != NULL)
    {
      if((int)(strlen(result) - strlen(term) - leadlen) < len)
        len = strlen(result) - strlen(term) - leadlen;

      strncpy(storage, result + leadlen , len);
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
  char stringBuf[SIZE_BUF+1];

  // HENRY - need to rewind to get instrument info
  ramp_fseek(pFI, 0 , SEEK_SET);
  
   if ((output = (InstrumentStruct *) calloc(1,sizeof(InstrumentStruct))) == NULL)
   {
      printf("Cannot allocate memory\n");
      return NULL;
   } else {
      const char *cpUnknown="UNKNOWN";
      strncpy(output->analyzer,cpUnknown,sizeof(output->analyzer));
      strncpy(output->detector,cpUnknown,sizeof(output->detector));
      strncpy(output->ionisation,cpUnknown,sizeof(output->ionisation));
      strncpy(output->manufacturer,cpUnknown,sizeof(output->manufacturer));
      strncpy(output->model,cpUnknown,sizeof(output->model));
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
          if(result != NULL && setTagValue(result, output->manufacturer, INSTRUMENT_LENGTH, isAncient?"manufacturer=":"value="))
	      found[0] = 1;
      
        }
        if(! found[1])
        {
           result = strstr(stringBuf, isAncient?"model=":"<msModel");
          if(result != NULL && setTagValue(result, output->model, INSTRUMENT_LENGTH, isAncient?"model=":"value="))
	      found[1] = 1;
        }
        if(! found[2])
        {
          result = strstr(stringBuf, isAncient?"ionisation=":"<msIonisation");
          if(result != NULL && setTagValue(result, output->ionisation, INSTRUMENT_LENGTH, isAncient?"ionisation=":"value="))
	      found[2] = 1;
        }
        if(! found[3])
        {
           result = strstr(stringBuf, isAncient?"msType=":"<msMassAnalyzer");
          if(result != NULL && setTagValue(result, output->analyzer, INSTRUMENT_LENGTH, isAncient?"msType=":"value="))
	      found[3] = 1;
        }
        if(! found[4])
        {
          result = strstr(stringBuf, "<msDetector");
          if(result != NULL && setTagValue(result, output->detector, INSTRUMENT_LENGTH, "value="))
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

// Cache support

// Get a new cache instance with a specified window size.  A larger window
// requires more memory, obviously.  Too small a window for the required
// function can lead perf comparable to no caching at all.

// Iterating over a range of 200 scans with a cache that contains 100 or
// fewer scans will yield no cache hits on each iteration.  Pick a cache
// size slightly larger than the window you expect to cover.

struct ScanCacheStruct *getScanCache(int size)
{
    struct ScanCacheStruct* cache = (struct ScanCacheStruct*) malloc(sizeof(struct ScanCacheStruct));
    cache->seqNumStart = 0;
    cache->size = size;
    cache->headers = (struct ScanHeaderStruct*) calloc(size, sizeof(struct ScanHeaderStruct));
    cache->peaks = (RAMPREAL**) calloc(size, sizeof(RAMPREAL*));
    return cache;
}

// Free all memory associated with a cache struct.
void freeScanCache(struct ScanCacheStruct* cache)
{
   if (cache) {
    int i;
    for (i = 0; i < cache->size; i++)
    {
        if (cache->peaks[i] != NULL)
            free(cache->peaks[i]);
    }
    free(cache->peaks);
    free(cache->headers);
    free(cache);
}
}

// Clear all cached values, freeing peaks, but not the cache arrays themselves.
void clearScanCache(struct ScanCacheStruct* cache)
{
    int i;
    for (i = 0; i < cache->size; i++)
    {
        if (cache->peaks[i] == NULL)
            continue;

        free(cache->peaks[i]);
        cache->peaks[i] = NULL;
    }
    memset(cache->headers, 0, cache->size * sizeof(struct ScanHeaderStruct));
}

// Shift the cache start index by a number of scans.  This moves the cache
// window left (negative) or right (positive).
void shiftScanCache(struct ScanCacheStruct* cache, int nScans)
{
    int i;
    cache->seqNumStart += nScans;
    if (abs(nScans) > cache->size)
    {
        // If the shift is larger than the size of the cache window,
        // just clear the whole cache.
        clearScanCache(cache);
    }
    else if (nScans > 0)
    {
        // Shifting window to the right.  Memory moves right, with new
        // empty scans on the end.

        // Free the peaks that memmove will overwrite.
        for (i = 0; i < nScans; i++)
        {
            if (cache->peaks[i] != NULL)
                free(cache->peaks[i]);
        }
        memmove(cache->peaks, cache->peaks + nScans,
            (cache->size - nScans) * sizeof(RAMPREAL*));
        memset(cache->peaks + cache->size - nScans, 0, nScans * sizeof(RAMPREAL*));
        memmove(cache->headers, cache->headers + nScans,
            (cache->size - nScans) * sizeof(struct ScanHeaderStruct));
        memset(cache->headers + cache->size - nScans, 0, nScans * sizeof(struct ScanHeaderStruct));
    }
    else if (nScans < 0)
    {
        // Shifting window to the left.  Memory moves right, with new
        // empty scans at the beginning.
        nScans = -nScans;

        // Free the peaks that memmove will overwrite.
        for (i = 0; i < nScans; i++)
        {
            if (cache->peaks[cache->size - 1 - i] != NULL)
                free(cache->peaks[cache->size - 1 - i]);
        }
        memmove(cache->peaks + nScans, cache->peaks,
            (cache->size - nScans) * sizeof(RAMPREAL*));
        memset(cache->peaks, 0, nScans * sizeof(RAMPREAL*));
        memmove(cache->headers  + nScans, cache->headers,
            (cache->size - nScans) * sizeof(struct ScanHeaderStruct));
        memset(cache->headers, 0, nScans * sizeof(struct ScanHeaderStruct));
    }
}

// Convert a scan index into a cache index, adjusting the cache window
// if necessary.
int getCacheIndex(struct ScanCacheStruct* cache, int seqNum)
{
    int seqNumStart = cache->seqNumStart;
    int size = cache->size;

    // First access, just set the start to seqNum.
    if (seqNumStart == 0)
        cache->seqNumStart = seqNum;
    // If requested scan is less than cache start, shift cache window
    // left to start at requested scan.
    else if (seqNum < seqNumStart)
        shiftScanCache(cache, (int) (seqNum - seqNumStart));
    // If requested scan is greater than cache end, shift cache window
    // right so last entry is requested scan.
    else if (seqNum >= seqNumStart + size)
        shiftScanCache(cache, (int) (seqNum - (seqNumStart + size - 1)));

    return (int) (seqNum - cache->seqNumStart);
}

struct ScanHeaderStruct* readHeaderCached(struct ScanCacheStruct* cache, int seqNum, RAMPFILE* pFI, ramp_fileoffset_t lScanIndex)
{
    int i = getCacheIndex(cache, seqNum);
    if (cache->headers[i].msLevel == 0)
        readHeader(pFI, lScanIndex, cache->headers + i);
    return cache->headers + i;
}

int readMsLevelCached(struct ScanCacheStruct* cache, int seqNum, RAMPFILE* pFI, ramp_fileoffset_t lScanIndex)
{
    struct ScanHeaderStruct* header = readHeaderCached(cache, seqNum, pFI, lScanIndex);
    return header->msLevel;
}

RAMPREAL *readPeaksCached(struct ScanCacheStruct* cache, int seqNum, RAMPFILE* pFI, ramp_fileoffset_t lScanIndex)
{
    int i = getCacheIndex(cache, seqNum);
    if (cache->peaks[i] == NULL)
        cache->peaks[i] = readPeaks(pFI, lScanIndex);
    return cache->peaks[i];
}

