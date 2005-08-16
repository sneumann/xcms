/***************************************************************************
                             RAMP


Non sequential parser for mzXML files

                             -------------------
    begin                : Wed Oct 10
    copyright            : (C) 2003 by Pedrioli Patrick, ISB, Proteomics
    email                : ppatrick@student.ethz.ch
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#include "ramp.h"

#define SIZE_BUF 4096
#ifdef __MINGW32__
#define fseeko fseek
#endif

/****************************************************************
 * Utility functions					*
 ***************************************************************/

char* skipspace(char* pStr)
{
	while (isspace(*pStr))
		pStr++;
	if (*pStr == '\0')
		return NULL;
	return pStr;
}

/****************************************************************
 * Find the Offset of the index					*
 ***************************************************************/
off_t getIndexOffset(FILE * pFI)
{
   char indexOffsetTemp[SIZE_BUF];
   char *target = "<indexOffset>";
   char *start, *end;

   fseeko(pFI, (off_t)-120, SEEK_END);
   while( fgets(indexOffsetTemp, SIZE_BUF, pFI) )
   {
      if ( !( start = strstr( indexOffsetTemp, "<indexOffset>" ) ) )
         continue;
      start = strstr( start, ">" ) + 1;
      if ( !( end = strstr( start, "<" ) ) )
         break;
      *end = '\0';

      if (sizeof(off_t)==8)
         return (atoll(start));
      else
         return (atol(start));
   }

   return 0; /* Error: never found the index offset */
}

/****************************************************************
 * Reads the Scan index in a list				*
 * Returns pScanIndex which becomes property of the caller	*
 * pScanIndex is -1 terminated					*
 ***************************************************************/

off_t *readIndex(FILE * pFI,
                off_t indexOffset,
                int *iLastScan)
{
   int  n;
   int  reallocSize = 8000;    /* initial # of scan indexes to expect */
   int  scanOffsetLength;
   char *beginScanOffset, *endScanOffset;
   char stringBuf[SIZE_BUF], scanOffset[SIZE_BUF];
   off_t *pScanIndex;

   *iLastScan = 0;

   if ((pScanIndex = (off_t *) malloc(reallocSize * sizeof(off_t))) == NULL)
   {
      printf("Cannot allocate memory\n");
      exit(1);
   }

   fseeko(pFI, indexOffset, SEEK_SET);

   fgets(stringBuf, SIZE_BUF, pFI);
   while( !strstr( stringBuf , "<offset " ) && !strstr( stringBuf, "/index") )
      fgets(stringBuf, SIZE_BUF, pFI);

   n = 1;

   while (!strstr(stringBuf, "/index"))
   {
      (*iLastScan)++;

      if ((beginScanOffset = (char *) (strstr(stringBuf, ">") + 1)) == NULL)
      {
         fgets(stringBuf, SIZE_BUF, pFI);
         continue;
      }
      if ((endScanOffset = (char *) (strstr(beginScanOffset, "<"))) == NULL)
      {
         fgets(stringBuf, SIZE_BUF, pFI);
         continue;
      }

      scanOffsetLength = endScanOffset - beginScanOffset;
      strncpy(scanOffset, beginScanOffset, scanOffsetLength);
      scanOffset[scanOffsetLength] = '\0';

      if (sizeof(off_t)==8)
         pScanIndex[n] = atoll(scanOffset);
      else
         pScanIndex[n] = atol(scanOffset);

      n++;

      if (n == reallocSize)
      {
         off_t *pTmp;

         reallocSize = reallocSize + 500;

         pTmp = (off_t*) realloc(pScanIndex, reallocSize*sizeof(off_t));
         
         if (pTmp==NULL)
         {
            printf("Cannot allocate memory\n");
            exit(1);
         }
         else
            pScanIndex=pTmp;
      }

      fgets(stringBuf, SIZE_BUF, pFI);
   }

   pScanIndex[n] = -1;

   return (pScanIndex);
}


/*
 * Reads scan header information.
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE
 *    RETURNING !
 */
void readHeader(FILE * pFI,
                off_t lScanIndex,
                struct ScanHeaderStruct *scanHeader)
{
   char stringBuf[SIZE_BUF];
   char *pStr, *pStr2;

   fseeko(pFI, lScanIndex, SEEK_SET);

   /*
    * initialize defaults
    */
   scanHeader->msLevel = 0;
   scanHeader->peaksCount = 0;
   scanHeader->retentionTime = 0;
   scanHeader->lowMZ = 1.E6;
   scanHeader->highMZ = 0.0;
   scanHeader->basePeakMz = 0.0;
   scanHeader->basePeakIntensity = 0.0;
   scanHeader->totIonCurrent = 0.0;
   scanHeader->precursorMZ = 0.0;
   scanHeader->precursorCharge = 0;
   scanHeader->scanType[0] = '\0';
   
   while (fgets(stringBuf, SIZE_BUF, pFI))
   {
      if (strstr(stringBuf, "<peaks "))
         return;
      if ((pStr = (char *) strstr(stringBuf, "msLevel=\"")))
         sscanf(pStr + 9, "%d\"", &(scanHeader->msLevel));      /* +9 moves the length of msLevel=" */
      if ((pStr = (char *) strstr(stringBuf, "peaksCount=\"")))
         sscanf(pStr + 12, "%d\"", &(scanHeader->peaksCount));
      if ((pStr = (char *) strstr(stringBuf, "retentionTime=\"")))
         sscanf(pStr + 15, "PT%lfS\"", &(scanHeader->retentionTime));   /* stored in for PT????S */
      if ((pStr = (char *) strstr(stringBuf, "lowMz=\"")))
         sscanf(pStr + 7, "%lf\"", &(scanHeader->lowMZ));
      if ((pStr = (char *) strstr(stringBuf, "highMz=\"")))
         sscanf(pStr + 8, "%lf\"", &(scanHeader->highMZ));
      if ((pStr = (char *) strstr(stringBuf, "startMz=\"")) && 
	  scanHeader->lowMZ == 1.E6)
	sscanf(pStr + 9, "%lf\"", &(scanHeader->lowMZ));
      if ((pStr = (char *) strstr(stringBuf, "endMz=\"")) && 
	  scanHeader->highMZ == 0.0)
	sscanf(pStr + 7, "%lf\"", &(scanHeader->highMZ));
      if ((pStr = (char *) strstr(stringBuf, "basePeakMz=\"")))
         sscanf(pStr + 12, "%lf\"", &(scanHeader->basePeakMz));
      if ((pStr = (char *) strstr(stringBuf, "basePeakIntensity=\"")))
         sscanf(pStr + 19, "%lf\"", &(scanHeader->basePeakIntensity));
      if ((pStr = (char *) strstr(stringBuf, "totIonCurrent=\"")))
         sscanf(pStr + 15, "%lf\"", &(scanHeader->totIonCurrent));

      if ((pStr = (char *) strstr(stringBuf, "scanType=\""))) {
	pStr = pStr+10;
	if ((pStr2 = (char *) strstr(pStr, "\""))) {
	  memcpy(&(scanHeader->scanType), pStr, sizeof(char)*(pStr2-pStr));
	  scanHeader->scanType[pStr2-pStr] = '\0';
	}
      }

      /*
       * read precursor mass
       */
      if ((pStr = (char *) strstr(stringBuf, "<precursorMz ")))
      {
         for (;;)
         {
            /*
             * Check for precursor charge.
             */
            if ((pStr2 = (char *) strstr(pStr, "precursorCharge=\"")))
            {
               sscanf(pStr2 + 17, "%d\"", &(scanHeader->precursorCharge));
//               printf("precursorCharge = %d\n", scanHeader->precursorCharge);
	    }
            /*
             * Find end of tag.
             */
            if ((pStr = strchr(pStr, '>')))
		break;

            fgets(stringBuf, SIZE_BUF, pFI);
            pStr = stringBuf;
         }

         pStr++;      // Skip >
            
         /*
          * Skip past white space.
          */
         while (!(pStr = skipspace(pStr)))
         {
            fgets(stringBuf, SIZE_BUF, pFI);
            pStr = stringBuf;
         }

         sscanf(pStr, "%lf<", &(scanHeader->precursorMZ));
//         printf("precursorMass = %lf\n", scanHeader->precursorMZ);
      }
   }
}


/****************************************************************
 * Reads the MS level of the scan.				*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

int readMsLevel(FILE * pFI,
      off_t lScanIndex)
{
   int  level;
   int  msLevelLen;
   char stringBuf[SIZE_BUF];
   char szLevel[12];
   char *beginMsLevel, *endMsLevel;

   fseeko(pFI, lScanIndex, SEEK_SET);

   fgets(stringBuf, SIZE_BUF, pFI);

   while (!(beginMsLevel = (char *) strstr(stringBuf, "msLevel=\"")))
   {
      fgets(stringBuf, SIZE_BUF, pFI);
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
 * Returns NULL if startMz was not set. Don't free the memory!	*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

double readStartMz(FILE * pFI,
		   off_t lScanIndex)
{
  char stringBuf[SIZE_BUF];
  double startMz = 1.E6;
  char *pStr;

  fseeko(pFI, lScanIndex, SEEK_SET);
   
  while (fgets(stringBuf, SIZE_BUF, pFI))
  {
      if (strstr(stringBuf, "<peaks "))
        return startMz;
      if ((pStr = strstr(stringBuf, "startMz=\""))){
        sscanf(pStr + 9, "%lf\"", &startMz);
        break;
      }
  }

  return startMz;
}


/****************************************************************
 * Reads startMz and endMz of the scan.				*
 * Returns NULL if startMz was not set. Don't free the memory!	*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

double readEndMz(FILE * pFI,
		   off_t lScanIndex)
{
  char stringBuf[SIZE_BUF];
  double endMz = 0.0;
  char *pStr;

  fseeko(pFI, lScanIndex, SEEK_SET);
   
  while (fgets(stringBuf, SIZE_BUF, pFI))
  {
      if (strstr(stringBuf, "<peaks "))
        return endMz;
      if ((pStr = strstr(stringBuf, "endMz=\""))){
        sscanf(pStr + 9, "%lf\"", &endMz);
        break;
      }
  }

  return endMz;
}


/****************************************************************
 * Reads RT of the scan.					*
 * Return a char* that becomes property of the caller!		*
 * Returns NULL if the RT was not set. Don't free the memory!	*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

char *readRT(FILE * pFI,
      off_t lScanIndex)
{
   int  RTLen;
   char stringBuf[SIZE_BUF];
   char *beginRT, *endRT;
   char *pRT;

   fseeko(pFI, lScanIndex, SEEK_SET);
 
   fgets(stringBuf, SIZE_BUF, pFI);

   while (!(beginRT = (char *) strstr(stringBuf, "retentionTime=\"")))
   {
      if (strstr(stringBuf, "retentionTime=\""))
         return NULL;
      fgets(stringBuf, SIZE_BUF, pFI);
   }

   beginRT += 15;
   endRT = (char *) strstr(beginRT, "\"");
   RTLen = endRT - beginRT;

   if ((pRT = (char *) malloc(RTLen * sizeof(char) + 1)) == NULL)
   {
      printf("Cannot allocate memory\n");
      exit(1);
   }

   strncpy(pRT, beginRT, RTLen);

   pRT[RTLen] = '\0';

   return pRT;
}


/****************************************************************
 * READS the base64 encoded list of peaks.			*
 * Return a float* that becomes property of the caller!		*
 * The list is terminated by -1					*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

float *readPeaks(FILE * pFI,
      off_t lScanIndex)
{
   int  n;
   int  peaksCount;
   int  peaksCountLength;
   int  peaksLen;       // The length of the base64 section
   int  trail;
   char szPeaksCount[50];
   char stringBuf[SIZE_BUF];
   char *beginPeaksCount, *endPeaksCount;
   float *pPeaks;

   char *pData;
   char *pBeginData;
   char *pDecoded;

   fseeko(pFI, lScanIndex, SEEK_SET);

   // Get the num of peaks in the scan and allocate the space we need
   fgets(stringBuf, SIZE_BUF, pFI);
   while (!(beginPeaksCount = (char *) strstr(stringBuf, "peaksCount=\"")))
   {
      fgets(stringBuf, SIZE_BUF, pFI);
   }

   // We need to move forward the length of peaksCount="
   beginPeaksCount += 12;

   endPeaksCount = (char *) strstr(beginPeaksCount, "\"");
   peaksCountLength = endPeaksCount - beginPeaksCount;
   strncpy(szPeaksCount, beginPeaksCount, peaksCountLength);
   szPeaksCount[peaksCountLength] = '\0';
   peaksCount = atoi(szPeaksCount);

   // We add 100 to compensate for initial white spaces and the opening of the element
   peaksLen = peaksCount * (64 / 3) * sizeof(char) + 100;
   if (peaksCount == 0)
     { // No peaks in this scan!!
       return NULL;
     }

   if ((pData = (char *) malloc(1 + peaksLen)) == NULL)
   {
      printf("Cannot allocate memory\n");
      exit(1);
   }

   fgets(pData, peaksLen, pFI);
   while (!(pBeginData = (char *) strstr(pData, "peaks")))
   {
      fgets(pData, peaksLen, pFI);
   }

   // skip <peaks precision="xx">
   //pBeginData += 21;
   while( !(pBeginData = strstr( pData , ">" )))
   {
      fgets(pData , peaksLen , pFI);
   }
   pBeginData++;	// skip the >

   // Chop after the end of the actual peak list
   trail = (peaksCount % 3);
   if( !trail )
   {
      pBeginData[ ((peaksCount * 32/3)* sizeof(char)) ] = '\0';
   }
   else
   {
      pBeginData[ (( (peaksCount * 32/3) + trail + 1)* sizeof(char)) ] = '\0';
   }
   //  printf( "%s\n" , pData );
   //printf( "%s\n" , pBeginData );

   if ((pPeaks = (float *) malloc(peaksCount * 8 + 1)) == NULL ||
       (pDecoded = (char *) malloc(peaksCount * 8 + 1)) == NULL)
   {
      printf("Cannot allocate memory\n");
      exit(1);
   }

   // Base64 decoding
   b64_decode_mio(pDecoded, pBeginData);


   // And byte order correction

   for (n = 0; n < (2 * peaksCount); n++)
   {
      ((uint32_t *) pPeaks)[n] =
         ntohl((uint32_t) ((uint32_t *) pDecoded)[n]);
   }
   pPeaks[n] = -1;

   free(pData);
   free(pDecoded);

   return (pPeaks);
}


/*
 * walk through each scan to find overall lowMZ, highMZ
 * sets overall start and end times also
 */
void readRunHeader(FILE * pFI,
                   off_t *pScanIndex,
                   struct RunHeaderStruct *runHeader,
                   int iLastScan)
{
   int  i;
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
   runHeader->startMZ = readStartMz( pFI , pScanIndex[i] );
   runHeader->endMZ = readEndMz( pFI , pScanIndex[i] );

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



int setTagValue(char* text,
   char* storage,
   int maxlen,
   char* lead,
   char* tail)
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
      if(strlen(result) - strlen(term) - strlen(lead) < len)
	len = strlen(result) - strlen(term) - strlen(lead);

      strncpy(storage, result + strlen(lead), len);
      storage[len] = 0;
      return 1;
    } // if term
  }
  return 0;
}


InstrumentStruct* getInstrumentStruct(FILE* pFI)
{
  InstrumentStruct* output = NULL;
  char* result = NULL;
  char* term = NULL;
  int found[] = {0, 0, 0, 0, 0};
  char stringBuf[SIZE_BUF];

   if ((output = (InstrumentStruct *) malloc(sizeof(InstrumentStruct))) == NULL)
   {
      printf("Cannot allocate memory\n");
      exit(1);
   }
   output->manufacturer[0] = 0;
   output->model[0] = 0;
   output->ionisation[0] = 0;
   output->analyzer[0] = 0;
   output->detector[0] = 0;

   fgets(stringBuf, SIZE_BUF, pFI);


   while( !strstr( stringBuf , "<msInstrument" ) &&  ! strstr(stringBuf, "<dataProcessing") && !feof(pFI))  /* this should not be needed if index offset points to correct location */
   {
      fgets(stringBuf, SIZE_BUF, pFI);
   }
   while(! strstr(stringBuf, "</msInstrument") &&  ! strstr(stringBuf, "</dataProcessing") && !feof(pFI))
   {
     if(! found[0])
     {
       result = strstr(stringBuf, "<msManufacturer");
       if(result != NULL && setTagValue(result, output->manufacturer, INSTRUMENT_LENGTH, "value=\"", "\""))
	   found[0] = 1;
     }
     if(! found[1])
     {
       result = strstr(stringBuf, "<msModel");
       if(result != NULL && setTagValue(result, output->model, INSTRUMENT_LENGTH, "value=\"", "\""))
	   found[1] = 1;
     }
     if(! found[2])
     {
       result = strstr(stringBuf, "<msIonisation");
       if(result != NULL && setTagValue(result, output->ionisation, INSTRUMENT_LENGTH, "value=\"", "\""))
	   found[2] = 1;
     }
     if(! found[3])
     {
       result = strstr(stringBuf, "<msMassAnalyzer");
       if(result != NULL && setTagValue(result, output->analyzer, INSTRUMENT_LENGTH, "value=\"", "\""))
	   found[3] = 1;
     }
     if(! found[4])
     {
       result = strstr(stringBuf, "<msDetector");
       if(result != NULL && setTagValue(result, output->detector, INSTRUMENT_LENGTH, "value=\"", "\""))
	   found[4] = 1;
     }
     fgets(stringBuf, SIZE_BUF, pFI);

   } // while

   if(found[0] || found[1] || found[2] || found[3] || found[4])
     return output;

   return NULL; // no data
}
