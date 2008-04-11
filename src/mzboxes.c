#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "R.h"
#include "Rdefines.h"

#define DEBUG
#undef DEBUG 

#define RAMPREAL double

#undef TRUE
#define TRUE    1
#undef FALSE
#define FALSE    0

#define SIZE_PEAKBUFS 200000   // max. # of "short" boxes
#define SIZE_PEAKBUFL 5000    // max. # of "long" boxes
#define DIM_PEAKBUFS 50       // max length of "short" boxes
#define DIM_PEAKBUFL 2500     // max length of "long" boxes = max number of spectra
#define DIM_MZBUF    10000     // max # m/z values to hold
#define DIM_SCANBUF 10000      // max # m/z values per scan

#define UNDEF_BUF      0
#define IN_SHORT_BUF  -1
#define IN_LONG_BUF    3

typedef unsigned int bufnrType;

struct scanStruct
{
   double  mz         [DIM_SCANBUF] ;
   double  intensity  [DIM_SCANBUF] ;
   int     length ;
} scanbuf;

struct pickOptionsStruct
{
   int    minEntries;
   float dev;
} pickOptions;

struct peakbufStruct {
  bufnrType snum  [SIZE_PEAKBUFS]; // actually used length of each buffer
  int       sscan [SIZE_PEAKBUFS][DIM_PEAKBUFS];
  RAMPREAL  smz   [SIZE_PEAKBUFS][DIM_PEAKBUFS];
  RAMPREAL  sint  [SIZE_PEAKBUFS][DIM_PEAKBUFS];
  
  bufnrType lnum  [SIZE_PEAKBUFL]; // actually used length of each buffer
  int       lscan [SIZE_PEAKBUFL][DIM_PEAKBUFL];
  RAMPREAL  lmz   [SIZE_PEAKBUFL][DIM_PEAKBUFL];
  RAMPREAL        lint  [SIZE_PEAKBUFL][DIM_PEAKBUFL];
  
  int             PeaksInBuf;
  unsigned char   freelists [SIZE_PEAKBUFS];
  unsigned char   freelistl [SIZE_PEAKBUFL];
} peakbuf;

struct mzvalStruct {
  RAMPREAL          mz    [DIM_MZBUF];
  bufnrType         bufnr [DIM_MZBUF];  // address in peakbuf
  int               slbuf [DIM_MZBUF];  // is in small OR long buffer
  int               length;             // actual # of m/z values in mz
} mzval;

void insertMZ(RAMPREAL val, int i, bufnrType bufpos, int which_buf, struct mzvalStruct *psbuf)
{
   if ((psbuf->length) >= DIM_MZBUF) error("MZ BUF SIZE too small ! \n"); 
   int n = psbuf->length - i;
   if (n>0) {
    memmove(psbuf->mz + i +1, psbuf->mz + i, n*sizeof(RAMPREAL));
    memmove(psbuf->bufnr + i +1, psbuf->bufnr + i, n*sizeof(bufnrType));
    memmove(psbuf->slbuf + i +1, psbuf->slbuf + i, n*sizeof(int));
   }  // else { printf("psbuf->length: %d , n: %d, insertMZ failed ! \n",psbuf->length,n); exit(-5);}
   
   psbuf->mz[i] = val; 
   psbuf->bufnr[i] = bufpos;
   psbuf->slbuf[i] = which_buf;  
   psbuf->length++;
}

void deleteMZ(int i, struct mzvalStruct *psbuf, struct peakbufStruct *peakbuf, int removeValues) 
{
  if (removeValues == TRUE) {
    bufnrType bufpos   = psbuf->bufnr[i];
    int       whichbuf = psbuf->slbuf[i];
    if (whichbuf == UNDEF_BUF) error("UNDEF_BUF error! \n"); 
    
    if (whichbuf == IN_SHORT_BUF) 
    {
      peakbuf->snum[bufpos] = 0;
      peakbuf->freelists[bufpos] = TRUE; // mark as free
    } else   
    {
      peakbuf->lnum[bufpos] = 0;
      peakbuf->freelistl[bufpos] = TRUE; // mark as free
    }
    peakbuf->PeaksInBuf--;
  } 
   
  if (i >= 0 && i < psbuf->length) {  
    int n = psbuf->length - 1 - i;
    if (n>0) {
      memmove (psbuf->mz + i, psbuf->mz + i+1, n*sizeof(RAMPREAL));
      memmove (psbuf->bufnr + i, psbuf->bufnr + i+1, n*sizeof(bufnrType));
      memmove (psbuf->slbuf + i, psbuf->slbuf + i+1, n*sizeof(int));
    }
  } 
  psbuf->length--;
}
  
int lower_bound(RAMPREAL val,struct mzvalStruct *mzval,int first, int length){
int half,mid;  // mzval->mz[first]
  while (length > 0) {
    half = length >> 1;
    mid = first;
    mid += half;
    if ( mzval->mz[mid] < val){
      first = mid;
      first ++;
      length = length - half -1;
    }
    else length = half;
  }
  return(first);
}

int upper_bound(RAMPREAL val,struct mzvalStruct *mzval,int first, int length){
int half,mid;  
  while (length > 0) {
    half = length >> 1;
    mid = first;
    mid += half;
    if (val < mzval->mz[mid]){
      length = half;
    }
    else {
      first = mid;
      first ++;
      length = length - half -1;
    }
  }
  return(first);
}

int lowerBound(double val,double *mzval,int first, int length){
int half,mid; 
  while (length > 0) {
    half = length >> 1;
    mid = first;
    mid += half;
    if ( mzval[mid] < val){
      first = mid;
      first ++;
      length = length - half -1;
    }
    else length = half;
  }
  return(first);
}

int upperBound(double val,double *mzval,int first, int length){
int half,mid;  
  while (length > 0) {
    half = length >> 1;
    mid = first;
    mid += half;
    if (val < mzval[mid]){
      length = half;
    }
    else {
      first = mid;
      first ++;
      length = length - half -1;
    }
  }
  return(first);
}

bufnrType getFreeBufPos(const int which_buf, struct peakbufStruct *peakbuf){
  bufnrType pos=0;
  if (which_buf == IN_SHORT_BUF) {
    while((peakbuf->freelists[pos] == FALSE) && pos < SIZE_PEAKBUFS) pos++;
    if (pos >= SIZE_PEAKBUFS-1) error("SIZE_PEAKBUFS too small ! \n"); 
  } else
  {
    while((peakbuf->freelistl[pos] == FALSE) && pos < SIZE_PEAKBUFL) pos++;
    if (pos >= SIZE_PEAKBUFL-1) error("SIZE_PEAKBUFL too small ! \n"); 
  }
  
  return(pos);
}

void insertpeak(RAMPREAL fMass,RAMPREAL fInten,const int scan, struct peakbufStruct *peakbuf,struct mzvalStruct *mzval,struct pickOptionsStruct pickOptions)
{
  int i,j,wasfound=FALSE;
  bufnrType num;
  RAMPREAL mzmean=0;
  
  for (i=0; i < mzval->length; i++) // TODO: Ausnutzen der Sortierungseigenschaft -> binary search, direkte Suche nach (fMass +- (pickOptions.dev * (float) fMass))
   { 
    double ddiff = fabs( mzval->mz[i] -  fMass);
    double ddev = (pickOptions.dev *  fMass);
    
    if (ddiff <= ddev) 
    { // match -> extend this box
      wasfound = TRUE;
      if (mzval->slbuf[i] == UNDEF_BUF) error("UNDEF_BUF error! \n"); 
      if (mzval->slbuf[i] == IN_SHORT_BUF) 
      { // APPEND TO SHORT_BUF
        bufnrType bnr = mzval->bufnr[i];
        num = peakbuf->snum[bnr];
        if (num >= DIM_PEAKBUFS-1) // its getting too long for small buf, we got to move it, move it !
          { 
             bufnrType newbufpos = getFreeBufPos(IN_LONG_BUF,peakbuf);
             mzval->bufnr[i]  = newbufpos;
             mzval->slbuf[i]  = IN_LONG_BUF;
             peakbuf->freelistl[newbufpos] = FALSE; 
             // now copy to long buf
             peakbuf->lnum[newbufpos] = num;
             memmove(peakbuf->lint[newbufpos] , peakbuf->sint[bnr], num*sizeof(RAMPREAL));
             memmove(peakbuf->lmz[newbufpos] , peakbuf->smz[bnr], num*sizeof(RAMPREAL));
             memmove(peakbuf->lscan[newbufpos] , peakbuf->sscan[bnr], num*sizeof(int));
             // free short buf pos
             peakbuf->snum[bnr] = 0;
             peakbuf->freelists[bnr] = TRUE; 
          } else 
          {  // still in short buf, so extend it
            peakbuf->sint[bnr][num] = fInten;
            peakbuf->smz[bnr][num] =  fMass;
            peakbuf->sscan[bnr][num] = scan;
            peakbuf->snum[bnr]++;
            for (j=0,mzmean=0.0;j<=num;j++) mzmean += peakbuf->smz[bnr][j];
            mzmean /= num+1;
          } 
      }
      
      if (mzval->slbuf[i] == IN_LONG_BUF) 
      { // APPEND TO LONG_BUF
        bufnrType bnr = mzval->bufnr[i];
        num = peakbuf->lnum[bnr];
        if (num >= DIM_PEAKBUFL-1)  error("PEAKBUFL BUFFER OVERFLOW ERROR! \n"); 
        peakbuf->lint[bnr][num] = fInten;
        peakbuf->lmz[bnr][num] =  fMass;
        peakbuf->lscan[bnr][num] = scan;
        peakbuf->lnum[bnr]++;
        for (j=0,mzmean=0.0;j<=num;j++) mzmean += peakbuf->lmz[bnr][j];
        mzmean /= num+1; 
      } 
        mzval->mz[i] = mzmean; // update mzval
    } 
   } // for
   
   // if not found
   if (wasfound == FALSE) {  // no, create new entry for mz
      int ipos = lower_bound(fMass,mzval,0,mzval->length);
      // insert into mzval->mz at pos ipos
      // and add values at free buf pos
      // printf("ipos: %d\n",ipos);
      bufnrType bufpos = getFreeBufPos(IN_SHORT_BUF,peakbuf);
      insertMZ(fMass,ipos,bufpos,IN_SHORT_BUF,mzval);
      
      peakbuf->sscan[bufpos][0] = scan;
      peakbuf->smz[bufpos][0]   = fMass;
      peakbuf->sint[bufpos][0]  = fInten;
      peakbuf->snum[bufpos] = 1;
      peakbuf->freelists[bufpos] = FALSE; // mark as occupied
      peakbuf->PeaksInBuf++;
    }
   
}

void insertscan(struct scanStruct *scanbuf ,const int scan, struct peakbufStruct *peakbuf,struct mzvalStruct *mzval, struct pickOptionsStruct pickOptions)
{
  int p;
  RAMPREAL fMass,lastMass=-1;
  RAMPREAL fInten;   
  //unsigned short size;
  
  if ((mzval->length == 0) && (peakbuf->PeaksInBuf == 0) ){ // copy first scan directly
    for (p=0;p < scanbuf->length;p++)
    { 
       fMass  = scanbuf->mz[p];
       fInten = scanbuf->intensity[p];
     
       mzval->mz[p] = fMass; 
       mzval->slbuf[p] = IN_SHORT_BUF;
       //size = (unsigned short) peakbuf->snum[p];
       peakbuf->sscan[p][0] = scan;
       peakbuf->smz[p][0]   = fMass;
       peakbuf->sint[p][0]  = fInten;
       peakbuf->snum[p] = 1;
       peakbuf->freelists[p] = FALSE; // mark as occupied
        
       mzval->bufnr[p] = p;
       mzval->length++;
    }
       peakbuf->PeaksInBuf = p;    
  } else 
  {
    for (p=0;p < scanbuf->length;p++)
    { 
       fMass  = scanbuf->mz[p];
       fInten = scanbuf->intensity[p];   
       // Test auf Sortiertheit
       if (fMass < lastMass)  error("m/z sort assumption violated ! \n");
       lastMass = fMass; //mzval->length
       insertpeak(fMass,fInten,scan,peakbuf,mzval,pickOptions);
    }
  }
}     

void cleanup(const int ctScan, struct peakbufStruct *peakbuf,struct mzvalStruct *mzval, int *scerr, struct pickOptionsStruct pickOptions, const int idebug){
int i;
bufnrType bufnr,entries;
 // check all peaks in mzval
 for (i=0; i < (mzval->length); i++) {
    bufnr = mzval->bufnr[i];
    int lastscan=0;
    if (mzval->slbuf[i] == IN_SHORT_BUF) 
      {
        entries = peakbuf->snum[bufnr];
        if (entries > 0) lastscan= peakbuf->sscan[bufnr][entries-1]; else error("(entries <= 0 ?)  err ! \n");
      } else
      {
        entries = peakbuf->lnum[bufnr];
        if (entries > 0) lastscan= peakbuf->lscan[bufnr][entries-1]; else error("(entries <= 0 ?)  err ! \n"); 
      }
      
  
  // finished (entries >= minEntries)  or just extended
  if ((entries >= pickOptions.minEntries) || (lastscan == ctScan)) // good feature
   { // is it finished ?
     if ((entries >= pickOptions.minEntries) && (lastscan < ctScan) )
       deleteMZ(i,mzval,peakbuf,FALSE); // delete index, keep values in buffer
       if (entries > ctScan) {
          if (idebug == TRUE) printf("Warning : entries > ctScan (is this centroid data ?) i: %d m: %3.4f  %d entries, lastscan %d   (ctScan=%d)\n",i,mzval->mz[i],entries,lastscan,ctScan); 
          (*scerr)++;
       }   
   }
   else 
   {
      deleteMZ(i,mzval,peakbuf,TRUE); // delete index & values 
   }
 
 } // for i
 
}

void getscan(int *pscan, double *pmz, double *pintensity, int *pscanindex,int *nmz, int *lastScan, double *pscanbufmz,double *pscanbufint,int *length) {
  int idx1,idx2;
  if (*pscan == *lastScan) {
    idx1 =  pscanindex[*pscan -1] +1;
    idx2 =  *nmz-1;
  } else
  {
    idx1 =  pscanindex[*pscan -1] +1;
    idx2 =  pscanindex[*pscan];
  }
  #ifdef DEBUG 
      printf("idx1 %d   idx2   %d  \n", idx1,idx2);  // Rprintf("Mx: %f  My: %f  \n",Mx,My);
  #endif
  int idx,i=0;
  if ((idx2 -idx1) > DIM_SCANBUF -1) error("SCANBUF too small ! \n"); 
  for (idx=idx1;idx <= idx2; idx++) 
  {
    pscanbufmz[i]   = pmz[idx-1];
    pscanbufint[i]  = pintensity[idx-1];
    i++;
  }
  *length=i;
}

void getScan(int scan, double *pmz, double *pintensity, int *pscanindex,int nmz, int lastScan, struct scanStruct *pscanbuf) {
  int idx,idx1,idx2,i=0;
  idx1 =  pscanindex[scan -1] +1;
  
  if (scan == lastScan)  idx2 =  nmz-1;   else   idx2 =  pscanindex[scan];
  if ((idx2 -idx1) > DIM_SCANBUF -1) error("SCANBUF too small ! \n");
  
  for (idx=idx1;idx <= idx2; idx++) 
  {
    pscanbuf->mz [i]   = pmz[idx-1];
    pscanbuf->intensity[i]  = pintensity[idx-1];
    i++;
  }
  pscanbuf->length=i;
}

double getScanEIC(int scan, double from, double to, double *pmz, double *pintensity, int *pscanindex,int nmz, int lastScan) {
  int idx,idx1,idx2;
  double sum=0.0;
  
  idx1 =  pscanindex[scan -1] +1;
  if (scan == lastScan)  idx2 =  nmz-1;  
    else idx2 =  pscanindex[scan];
    
    
  int idx1b = lowerBound(from, pmz, idx1-1, idx2-idx1); 
  int idx2b = upperBound(to, pmz, idx1b, idx2-idx1b+1);  
    
  //printf("idx1 %d idx2 %d idx1b %d idx2b %d \n",idx1,idx2,idx1b,idx2b);  
    
  for (idx=idx1b;idx <= idx2b; idx++) 
  {
    double mzval = pmz[idx-1];
    if ((mzval <= to) && (mzval >= from)) sum += pintensity[idx-1];
  }
  return(sum);
}

SEXP getEIC(SEXP mz, SEXP intensity, SEXP scanindex, SEXP massrange, SEXP scanrange, SEXP lastscan) {
  double *pmz, *pintensity,*p_vint, massrangeFrom,massrangeTo;
  int i,*pscanindex, *p_scan,scanrangeFrom, scanrangeTo,ilastScan,nmz,ctScan,buflength;
  SEXP list_names,reslist,vscan,vint;
  pmz = REAL(mz);  
  nmz = GET_LENGTH(mz);
  pintensity = REAL(intensity); 
  pscanindex = INTEGER(scanindex);
  int firstScan = 1;   // is always 1 
  ilastScan = INTEGER(lastscan)[0];
  massrangeFrom = REAL(massrange)[0];  
  massrangeTo =  REAL(massrange)[1];  
  scanrangeFrom = INTEGER(scanrange)[0];
  scanrangeTo = INTEGER(scanrange)[1];
  if ((scanrangeFrom <  firstScan) || (scanrangeFrom > ilastScan) || (scanrangeTo < firstScan) || (scanrangeTo > ilastScan))
     error("Error in scanrange \n"); 
  char *names[2] = {"scan", "intensity"}; 
  PROTECT(list_names = allocVector(STRSXP, 2));
  for(i = 0; i < 2; i++)  
    SET_STRING_ELT(list_names, i,  mkChar(names[i])); 
    
  buflength = scanrangeTo - scanrangeFrom +1;  
  PROTECT(reslist = allocVector(VECSXP, 2));
  PROTECT(vscan = NEW_INTEGER(buflength));   
  p_scan = INTEGER_POINTER(vscan);
  PROTECT(vint = NEW_NUMERIC(buflength));   
  p_vint = NUMERIC_POINTER(vint);
  
  i=0;
  for (ctScan=scanrangeFrom;ctScan<=scanrangeTo;ctScan++)
  {
    p_scan[i] = ctScan;
    p_vint[i] = getScanEIC(ctScan,massrangeFrom,massrangeTo,pmz, pintensity, pscanindex,nmz,ilastScan);
    i++;
  }
  
  SET_VECTOR_ELT(reslist, 0, vscan);// attaching integer vector scan to list
  SET_VECTOR_ELT(reslist, 1, vint); // attaching double vector m/z to list
  setAttrib(reslist, R_NamesSymbol, list_names); //and attaching the vector names
  
  UNPROTECT(4);
  return(reslist);
}

SEXP findmzboxes(SEXP mz, SEXP intensity, SEXP scanindex, SEXP massrange, SEXP scanrange, SEXP lastscan, SEXP dev, SEXP minEntries, SEXP debug) {
  double *pmz, *pintensity,*p_vmz,*p_vint, massrangeFrom,massrangeTo;
  int i,*pscanindex, *p_scan,scanrangeFrom, scanrangeTo,ctScan,nmz,lastScan, idebug;
  int scerr = 0;  // count of peak insertion errors, due to missing/bad centroidisation
  SEXP peaklist,entrylist,vscan,vmz,vint,list_names;
     
  pmz = REAL(mz);  
  nmz = GET_LENGTH(mz);
  pintensity = REAL(intensity); 
  pscanindex = INTEGER(scanindex);
  lastScan = INTEGER(lastscan)[0];
  idebug = INTEGER(debug)[0];
  
  pickOptions.dev = REAL(dev)[0];                         // 140e-6;
  pickOptions.minEntries = INTEGER(minEntries)[0];        //4;   
  massrangeFrom = REAL(massrange)[0];  
  massrangeTo =  REAL(massrange)[1];  
  scanrangeFrom = INTEGER(scanrange)[0];
  scanrangeTo = INTEGER(scanrange)[1];
  
  memset(peakbuf.snum, 0, sizeof(peakbuf.snum));
  memset(peakbuf.lnum, 0, sizeof(peakbuf.lnum));
  memset(peakbuf.freelists, TRUE, sizeof(peakbuf.freelists)); // TRUE means FREE
  memset(peakbuf.freelistl, TRUE, sizeof(peakbuf.freelistl));
  peakbuf.PeaksInBuf = 0;
  
  memset(mzval.mz, 0, sizeof(mzval.mz));
  memset(mzval.slbuf, UNDEF_BUF, sizeof(mzval.slbuf));
  mzval.length = 0;   
  
  memset(scanbuf.mz, 0, sizeof(scanbuf.mz));
  memset(scanbuf.intensity, 0, sizeof(scanbuf.intensity));
  scanbuf.length = 0;
  
  // a character string vector of the "names" attribute of the objects in our list
  char *names[3] = {"scan","mz", "intensity"}; 
  PROTECT(list_names = allocVector(STRSXP, 3));
  for(i = 0; i < 3; i++)
    SET_STRING_ELT(list_names, i,  mkChar(names[i])); 
   
  for (ctScan=scanrangeFrom;ctScan<=scanrangeTo;ctScan++)
  {
    getScan(ctScan, pmz, pintensity, pscanindex,nmz,lastScan, &scanbuf);
    if (scanbuf.length > 0) 
    {
      if (idebug == TRUE) 
          printf("Scan Nr. %d of %d (%d %%) %d peaks -- working at %d m/z boxes, %d boxes found.\n", ctScan, scanrangeTo,  (int)100.0*ctScan/scanrangeTo,scanbuf.length,mzval.length,peakbuf.PeaksInBuf);
      insertscan(&scanbuf,ctScan,&peakbuf,&mzval,pickOptions);
      cleanup(ctScan,&peakbuf,&mzval,&scerr,pickOptions,idebug);
    }
  }
  cleanup(ctScan+1,&peakbuf,&mzval,&scerr,pickOptions,idebug); 
  
  PROTECT(peaklist = allocVector(VECSXP, peakbuf.PeaksInBuf));
  int total = 0;
  for (i=0;i<SIZE_PEAKBUFS;i++) {
    if (peakbuf.freelists[i] == FALSE) {
      bufnrType buflength = peakbuf.snum[i];
      PROTECT(entrylist = allocVector(VECSXP, 3));
      PROTECT(vscan = NEW_INTEGER(buflength));   
      p_scan = INTEGER_POINTER(vscan);
      PROTECT(vmz = NEW_NUMERIC(buflength));   
      p_vmz = NUMERIC_POINTER(vmz);
      PROTECT(vint = NEW_NUMERIC(buflength));   
      p_vint = NUMERIC_POINTER(vint);
      memmove(p_scan,peakbuf.sscan[i],buflength*sizeof(int));
      memmove(p_vmz,peakbuf.smz[i],buflength*sizeof(RAMPREAL));
      memmove(p_vint,peakbuf.sint[i],buflength*sizeof(RAMPREAL));
      SET_VECTOR_ELT(entrylist, 0, vscan);// attaching integer vector scan to list
      SET_VECTOR_ELT(entrylist, 1, vmz); // attaching double vector m/z to list
      SET_VECTOR_ELT(entrylist, 2, vint);// attaching double vector intensity to list
      setAttrib(entrylist, R_NamesSymbol, list_names); //and attaching the vector names
      SET_VECTOR_ELT(peaklist, total,entrylist);
      UNPROTECT(4); //entrylist,vmz,vint,vscan
     total++;
    }
  }
  for (i=0;i<SIZE_PEAKBUFL;i++) {
    if (peakbuf.freelistl[i] == FALSE) {
      bufnrType buflength = peakbuf.lnum[i];
      PROTECT(entrylist = allocVector(VECSXP, 3));
      PROTECT(vscan = NEW_INTEGER(buflength));   
      p_scan = INTEGER_POINTER(vscan);
      PROTECT(vmz = NEW_NUMERIC(buflength));   
      p_vmz = NUMERIC_POINTER(vmz);
      PROTECT(vint = NEW_NUMERIC(buflength));   
      p_vint = NUMERIC_POINTER(vint);
      memmove(p_scan,peakbuf.lscan[i],buflength*sizeof(int));
      memmove(p_vmz,peakbuf.lmz[i],buflength*sizeof(RAMPREAL));
      memmove(p_vint,peakbuf.lint[i],buflength*sizeof(RAMPREAL));
      SET_VECTOR_ELT(entrylist, 0, vscan);// attaching integer vector scan to list
      SET_VECTOR_ELT(entrylist, 1, vmz); // attaching double vector m/z to list
      SET_VECTOR_ELT(entrylist, 2, vint);// attaching double vector intensity to list
      setAttrib(entrylist, R_NamesSymbol, list_names); //and attaching the vector names
      SET_VECTOR_ELT(peaklist, total,entrylist);
      UNPROTECT(4); //entrylist,vmz,vint,vscan
     total++;
    }
  }
  if (scerr > 0) printf("Warning: There were %d peak data insertion problems. \n The reason might be missing/bad centroidisation or wrong parameter settings (\"dev\" too high?).\n", scerr);
  printf(" %d m/z boxes.\n", total);
  
  UNPROTECT(2); // peaklist,list_names
  return(peaklist);
}




