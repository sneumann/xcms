// STDLIB:
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include "string.h"

// Obiwarp
#include "vec.h"
#include "mat.h"
#include "lmat.h"
#include "xcms_dynprog.h"
#include "cmdparser.h"

// R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>


/********************************************/
//char * VERSION = "0.9.2";
/********************************************/

#define DEBUG (0)
/**************************************
bool format_is_labelless(const char *format);
************************************/

//XCMS****************
extern "C" SEXP R_set_from_xcms(SEXP valscantime, SEXP scantime, SEXP mzrange, SEXP mz, SEXP intensity, SEXP valscantime2, SEXP scantime2, SEXP mzrange2, SEXP mz2, SEXP intensity2, SEXP argv){
// einlesen von zwei Datensätzen im lmata-Format
//END XCMS****************


    // NOTE: use outfile as indicator if option passed in as opts.outfile!
    // because we can set opts.outfile to NULL and other routines will
    // automatically write to stdout!
    bool outfile = 0;
    bool outfile_is_stdout = 0;
//XCMS********************************************
    int len,i;
    len=length(argv);
    char*  pargv[len];
    PROTECT(argv = AS_CHARACTER(argv));
    for(i=0;i<len;i++){
      pargv[i] = R_alloc(strlen(CHAR(STRING_ELT(argv, i))), sizeof(char)); 
      strcpy(pargv[i], CHAR(STRING_ELT(argv, i))); 
      }
//END XCMS****************************************
    CmdParser opts(len, pargv, "0.9.2");

   
            outfile_is_stdout = 1;
            opts.outfile = NULL;

       // opts.print_argv();
    
/**************************************************
    char file1[1024];
    char file2[1024];
    strcpy(file1, opts.infiles[0]);
    strcpy(file2, opts.infiles[1]); 
***************************************************/

//XCMS*************
  int pvalscantime, pmzrange;
  int pvalscantime2, pmzrange2;
  double *pscantime, *pmz, *pintensity;
  double *pscantime2, *pmz2, *pintensity2;    
  SEXP corrected;
   
  valscantime = coerceVector(valscantime, INTSXP);
  mzrange = coerceVector(mzrange, INTSXP);
  
  valscantime2 = coerceVector(valscantime2, INTSXP);
  mzrange2 = coerceVector(mzrange2, INTSXP);
  //intensity = coerceVector(intensity, REALSXP);
  
  
   pvalscantime = INTEGER(valscantime)[0];
   pmzrange = INTEGER(mzrange)[0];
   pscantime = REAL(scantime);
   pmz = REAL(mz);
   pintensity = REAL(intensity);

   pvalscantime2 = INTEGER(valscantime2)[0];
   pmzrange2 = INTEGER(mzrange2)[0];
   pscantime2 = REAL(scantime2);
   pmz2 = REAL(mz2);
   pintensity2 = REAL(intensity2);
             
// END XCMS************

    // ************************************************************
    // * READ IN FILES TO GET MAT 
    // ************************************************************
    LMat lmat1;
    LMat lmat2;
    MatF smat;
    DynProg dyn;
//for(int i=0; i<length(t_scantime);i++)
//    printf(" %f ",t_pscantime[i]);

//XCMS*************
        lmat1.set_from_xcms(pvalscantime, pscantime, pmzrange, pmz, pintensity);
        lmat2.set_from_xcms(pvalscantime2, pscantime2, pmzrange2, pmz2, pintensity2);
        //set_from_xcms in lmat.cpp lmat.h
//lmat2.print_xcms();
    
    ////puts("LMAT1 AND LMAT2"); 
   //lmat1.print_xcms(); lmat2.print_xcms();
// END XCMS************


    // ************************************************************
    // * SCORE THE MATRICES
    // ************************************************************
    if (DEBUG) {
        std::cerr << "Scoring the mats!\n";
    }

  if (opts.smat_in != NULL) {
        smat.set_from_binary(opts.smat_in);
        dyn._smat = &smat;
    }
//XCMS****************
      // smat.set(pvalscantime,pmzrange,(float*)pintensity);
      //  dyn._smat = &smat;
//END XCMS *****************    

        dyn.score(*(lmat1.mat()), *(lmat2.mat()), smat, opts.score);
        // SETTING THE SMAT TO BE std normal
        if (!opts.nostdnrm) {
            if (!smat.all_equal()) { 
                smat.std_normal();
            }
        }
        if (!strcmp(opts.score,"euc")) {
            smat *= -1; // inverting euclidean
        }


 //smat.print();


/*****************************************    
    if (opts.smat_out != NULL) {
        std::cerr << "Writing binary smat to '" << opts.smat_out << "'\n";
        smat.write(opts.smat_out);
        //smat.print(smat_out_files[0]);
        exit(0);
    }
*********************************************/

    // ************************************************************
    // * PREPARE GAP PENALTY ARRAY
    // ************************************************************
   
    MatF time_tester;
    MatF time_tester_trans;
    VecF mpt;
    VecF npt;
    VecF mOut_tm;
    VecF nOut_tm;

    int gp_length = smat.rows() + smat.cols();
 //printf("smat.rows:%d\n",smat.rows());
    VecF gp_array;
    dyn.linear_less_before(opts.gap_extend,opts.gap_init,gp_length,gp_array);

    // ************************************************************
    // * DYNAMIC PROGRAM
    // ************************************************************ 
    int minimize = 0;
    if (DEBUG) {
        std::cerr << "Dynamic Time Warping Score Matrix!\n";
    }
    dyn.find_path(smat, gp_array, minimize, opts.factor_diag, opts.factor_gap, opts.local, opts.init_penalty);

    VecI mOut;
    VecI nOut;
    dyn.warp_map(mOut, nOut, opts.response, minimize);
    //puts("mOUT"); mOut.print(); nOut.print();

    // Major output unless its the only case where we don't need warped time
    // values
/*********************************************
    if (!(outfile_is_stdout && format_is_labelless(opts.format))) {
        // MAJOR OUTPUT:
        VecF nOutF;
        VecF mOutF;
        lmat1.tm_axis_vals(mOut, mOutF);
        lmat2.tm_axis_vals(nOut, nOutF); //
        lmat2.warp_tm(nOutF, mOutF); 
        lmat2.tm()->print(1);
    }
**************************************************/
    // No labels on matrix and we have an outfile to produce
    // Needs to be after MAJOR OUTPUT since it warps the data!
/**********************************************************
    if (format_is_labelless(opts.format) && outfile) {
        // @TODO: implement data warping here
    }
*********************************************************/

    // All subroutines below should write to the specified file
    // if the file == NULL then they should write to stdout!
    // opts.outfile is set to NULL if "STDOUT" is specified!
/*************************************  
  if (outfile) {
        if (!strcmp(opts.format, "mat")) {
            lmat2.mat()->write(opts.outfile);
        }
        else if (!strcmp(opts.format, "mata")) {
            lmat2.mat()->print(opts.outfile);
        }
        else if (!strcmp(opts.format, "lmat")) {
            lmat2.write(opts.outfile);
        }
        else if (!strcmp(opts.format, "lmata")) {
            lmat2.print(opts.outfile);
        }
        else {
            std::cerr << "Can't output to" << opts.format << "format (yet)\n";
            exit(0);
        }
    }
**************************************/
    // After all other output to stdout
/*********************************************
    if (opts.timefile != NULL) {
        time_tester.set_from_ascii(opts.timefile, 1);  // no headers on the files
        time_tester.transpose(time_tester_trans);
        mpt.set(time_tester_trans.cols(), time_tester_trans.pointer(0));
        npt.set(time_tester_trans.cols(), time_tester_trans.pointer(1));
        float ssr, asr, sad, aad;
        dyn.path_accuracy((*lmat1._tm), (*lmat2._tm), mOut, nOut, mpt, npt, ssr, asr, sad, aad);
        printf("%f %f %f %f\n", ssr, asr, sad, aad);
    }
****************************************************/


/*******************************************
    if (opts.images) {
        PngIO wrt(1);
        char base_fn[1024];
        strcpy(base_fn, "obi-warp_");
        char tb_fn[1024];
        strcpy(tb_fn, base_fn);
        strcat(tb_fn, "tb.png");
        //char *tb_fn = "tb.png";
        wrt.write(tb_fn, dyn._tb);
        char tbpath_fn[1024];
        strcpy(tbpath_fn, base_fn);
        strcat(tbpath_fn, "tbpath.png");
        wrt.write(tbpath_fn, dyn._tbpath);

        char asmat_fn[1024];
        strcpy(asmat_fn, base_fn);
        strcat(asmat_fn, "asmat.png");
        //wrt.write(asmat_fn, dyn._asmat);

        //strcpy(base_fn, "tb.png");
        //char *tbpath_fn = "tbpath.png";
        //char *tbscores_fn = "tbscores.png";
        //wrt.write(tbscores_fn, dyn._tbscores);
        //char *asmat_fn = "asmat.png";
        //wrt.write(asmat_fn, dyn._asmat);
        char *smat_fn = "smat.png";
        //wrt.write(smat_fn, *dyn._smat);
    }
**************************************************/

/*
   char silly[100];
   strcpy(silly, "png_");
   char tmpp[5];
   sprintf(tmpp, "%d", i);
   strcat(silly, tmpp); 
   strcat(silly, ".png");

   PngIO wrt(0);
//wrt.write(silly, dyn._tbpath);
wrt.write(silly, _scorepath);
*/
        VecF nOutF;
        VecF mOutF;
        lmat1.tm_axis_vals(mOut, mOutF);
        lmat2.tm_axis_vals(nOut, nOutF); 
        lmat2.warp_tm(nOutF, mOutF); 

        //s_lmat.tm()->print(1);

//printf("lmat2.tm()[0]:%f\n",lmat2.tm());
//printf("lmat2.tm()first:%f\n",lmat2.tm()->first());
//printf("lmat2.tm()last:%f\n",lmat2.tm()->last());
//for(int i; i < lmat2.tm()->len();i++){
   //printf("lmat2.tm()back[0]:%f\n",lmat2.tm()->back()[i]);}
//XCMS****************************************
PROTECT(corrected = allocVector(REALSXP, length(scantime2)));
for(int i; i < length(scantime2);i++){
   REAL(corrected)[i] = lmat2.tm()->back()[i];}
UNPROTECT(2);
//END XCMS************************************

//printf("STOP\n");
//lmat2.print_xcms();
//lmat2.mat()->print();

//printf("lo_tm():%f\n",lmat2.lo_tm());

//printf("output:%f\n",lmat2.mat());

//XCMS *********************
return corrected;
//END XCMS *****************

//return R_NilValue;
}
/***********************************************
bool format_is_labelless(const char *format) {
    if (!strcmp(format,"mat") || !strcmp(format,"mata")) {
        return 1;
    }
    else {
        return 0;
    }
}
**************************************************/
