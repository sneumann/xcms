
#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include "cmdparser.h"
#include "string.h"

float cor_gap_init = 0.3f;
float cor_gap_extend = 2.4f;
float cov_gap_init = 0.0f;
float cov_gap_extend = 11.7f;
float euc_gap_init = 0.9f;
float euc_gap_extend = 1.8f;
float prd_gap_init = 0.0f;
float prd_gap_extend = 7.8f;

CmdParser::CmdParser(int argc, char ** argv, char * version) {

    int diagnostics = 0;

    char * larg; // long argument
    char * sarg; // long argument
    _argv = argv;

    // IO:
    format = NULL;
    outfile = NULL;
    images = 0;
    timefile = NULL;

    // WARPING:
    score = "cor";
    local = 0;
    factor_diag = 2.f;
    factor_gap = 1.f;
    gap_init = 0.f;
    gap_extend = 0.f;
    init_penalty = 0.f;
    response = 1.f;
    nostdnrm = 0;

    //warp_data = 0;

    int infile_cnt = 0;
    infiles[0] = NULL;
    infiles[1] = NULL;

    //infiles;
    smat_out = NULL;
    smat_in = NULL;
    
    progname = "obiwarp";

    for (_i = 1; _i < argc; ++_i) {
        //printf("Considering argument: %s\n", argv[_i]);
        if (argv[_i][0] == '-') {
            if (argv[_i][1] == '-') {
                larg = &(argv[_i][2]);
                //printf("LONGARG: %s\n", larg);
                if (eq(larg,"format")) {set_string(&format);}
                else if (eq(larg,"outfile")) {set_string(&outfile);}
                else if (eq(larg,"images")) {set_flag(images);}
                else if (eq(larg,"timefile")) {set_string(&timefile);}
                else if (eq(larg,"score")) {set_string(&score);set_defaults_by_score(score);}
                else if (eq(larg,"local")) {set_flag(local);}
                else if (eq(larg,"nostdnrm")) {set_flag(nostdnrm);}
                else if (eq(larg,"factor")) {set_comma_list_float(factor_diag, factor_gap);}
                else if (eq(larg,"gap")) {set_comma_list_float(gap_init, gap_extend);}
                else if (eq(larg,"init")) {set_float(init_penalty);}
                else if (eq(larg,"response")) {set_float(response);}
                else if (eq(larg,"help")) {help_func(argc); exit(1);}
                else if (eq(larg,"version")) {print_version(version); exit(1);}
                else if (eq(larg,"smat_in")) {set_string(&smat_in);}
                else if (eq(larg,"smat_out")) {set_string(&smat_out);}
                else if (eq(larg, "diagnostics")) {diagnostics = 1;}
                else {
                    std::cerr << "********************************************\n";
                    std::cerr << "WARNING: invalid argument: " << argv[_i] << "\n";
                    std::cerr << "********************************************\n";
                }
            }
            else {  // short arg
                sarg = &(argv[_i][1]);
                if (*sarg == 'l') { set_flag(local); }
                else if (*sarg == 'o') { set_string(&outfile); }
                else if (*sarg == 't') { set_string(&timefile); }
                else if (*sarg == 's') { set_string(&score); set_defaults_by_score(score); }
                else if (*sarg == 'f') { set_comma_list_float(factor_diag, factor_gap); }
                else if (*sarg == 'g') { set_comma_list_float(gap_init, gap_extend); }
                else if (*sarg == 'i') { set_float(init_penalty); }
                else if (*sarg == 'r') { set_float(response);  }
                else if (*sarg == 'w') { set_string(&smat_in);  }
                else if (*sarg == 'x') { set_string(&smat_out);  }
                else if (*sarg == 'h') { help_func(argc); exit(1); }

                else if (*sarg == 'v') { print_version(version); exit(1); }
                else { 
                    std::cerr << "********************************************\n";
                    std::cerr << "WARNING: invalid argument: " << argv[_i] << "\n";
                    std::cerr << "********************************************\n";
                }
            }
        }
        else {  // must be an infile!
            infiles[infile_cnt] = argv[_i];
            ++infile_cnt;
        }

    }


/**********************************************************
    if (format == NULL && infiles[0] != NULL ) { get_format_from_file1(); }
    //printf("diagnostics: %d\n", diagnostics);
    if (diagnostics) {
        //printf("FORMATT!!! %s\n", format);
        std::cout <<  "**********************************************\n";
        if (format == NULL) { std::cout << "format: \n"; }
        else { std::cout << "format: "<< format << "\n"; }
        if (outfile == NULL) { std::cout << "outfile: \n"; }
        else { std::cout << "outfile: "<< outfile << "\n"; }
        std::cout << "images: " << images << "\n";
        if (timefile == NULL) { std::cout << "timefile: \n"; }
        else { std::cout << "timefile: "<< timefile << "\n"; }
        std::cout << "score: " << score << "\n";
        std::cout << "local: " << local << "\n";
        std::cout << "nostdnrm: " << nostdnrm << "\n";
        std::cout << "response: " << response << "\n";
        std::cout << "factor_diag: " << factor_diag << "\n";
        std::cout << "factor_gap: " << factor_gap << "\n";
        std::cout << "gap_init: " << gap_init << "\n";
        std::cout << "gap_extend: " << gap_extend<< "\n";
        std::cout << "init_penalty: " << init_penalty << "\n";
        if (infiles[0] == NULL) { std::cout << "file1: \n"; }
        else { std::cout << "file1: "<< infiles[0] << "\n"; }
        if (infiles[1] == NULL) { std::cout << "file2: \n"; }
        else { std::cout << "file2: "<< infiles[1] << "\n"; }
        std::cout << "**********************************************\n";
        exit(0);
    }
    verify_infiles();
****************************************/
}

// If format not supplied, we need to grab it from the file
void CmdParser::get_format_from_file1() {
    format = strrchr(infiles[0], '.');
    format++;
}

void CmdParser::set_flag(bool &val) {
    val = 1;
}

void CmdParser::set_string(char ** string_ptr) {
    ++_i;
    *string_ptr = _argv[_i];
}

void CmdParser::set_float(float &val) {
    ++_i;
    val = atof(_argv[_i]);
}

bool CmdParser::eq(char * first, char * sec) {
    return !strcmp(first, sec);
}

void CmdParser::set_comma_list_float(float &val1, float &val2) {
    ++_i;
    int sz = strlen(_argv[_i]);
    char *copied = new char[sz+1];
    strcpy(copied, _argv[_i]);
    char *comma_ptr;
    comma_ptr = strchr(copied, ',');
    *comma_ptr = '\0';
    val1 = atof(copied);
    val2 = atof(++comma_ptr);
    delete[] copied;
}


void CmdParser::verify_infiles() {
    int i;
    // If either file is null, give them the help!
    for (i = 0; i < 2; ++i) {
        if (infiles[i] == NULL) {
            help_func(1); exit(1);
        }
    }
    for (i = 0; i < 2; ++i) {
        if (file_is_readable(infiles[i])) {
        }
        else {
            printf("Cannot open %s for reading! aborting...", infiles[i]);
            exit(1);
        }
    }
}
bool CmdParser::file_is_readable(char * filename) {
    if (FILE * file = fopen(filename, "r")) {
        fclose(file);
        return true;
    }
    return false;
}


void CmdParser::set_defaults_by_score(char * score_arg) {
    if (!strcmp(score,"cor")) {
        gap_extend = cor_gap_extend;
        gap_init = cor_gap_init;
    }
    else if (!strcmp(score,"cov")) {
        gap_extend = cov_gap_extend;
        gap_init = cov_gap_init;
    }
    else if (!strcmp(score,"euc")) {
        gap_extend = euc_gap_extend;
        gap_init = euc_gap_init;
    }
    else if (!strcmp(score,"prd")) {
        gap_extend = prd_gap_extend;
        gap_init = prd_gap_init;
    }
}

void CmdParser::print_version(char * vers) {
    printf("%s\n", vers);
}
// needs to also print version!!!

//xcms****************
void CmdParser::print_argv() {
int i, j;
for(i=0; i< 4 ;i++)    
 //for(j=0; j< 2 ;j++) 
  printf("argv[%d]:%s\n",i, _argv[i]);
   
  printf("response:%f\n",response);
  printf("gab:%f,%f\n",gap_init,gap_extend); 
  }
//END XCMS******************

void CmdParser::help_func(int arg_cnt) {
    ++_i;
    //printf("_i: %d arg_cnt: %d\n", _i, arg_cnt);
    if (_i < arg_cnt && !strcmp(_argv[_i], "formats")) {
            std::cout << 
                // @TODO: NEED TO FINISH THIS GUY!
                //"       .mzXML  ISB mass spectrometry data format\n" <<
                //"       .mzData PSI mass spectrometry data format\n" <<
                //"       .ma     Matrix Ascii: ascii matrix\n" <<
                "FILE FORMATS:\n" <<
                "   .mat    Matrix Tagged: binary matrix with dim tags\n" <<
                "   .mata   Matrix Tagged Ascii: ascii matrix with dim tags\n" <<
                "   .lmat   Labeled Matrix Tagged: labeled binary matrix with dim tags\n" <<
                "   .lmata  Labeled Matrix Tagged Ascii: labeled ascii matrix with dim tags\n" <<
                "\n" <<
                "   Dim tags indicate the size of a vector or matrix.  These are essential for\n" <<
                "   parsing binary file formats and speed up the parsing of ascii files.\n" <<
                "\n" <<
                "   Labeled matrices contain vectors that give specific locations to each\n" <<
                "   row and column of matrix data.\n" <<
                "\n" <<
                "PASSING AS ARGUMENT:\n" <<
                "   The format of the files may be passed as an argument (with no leading '.'):\n" << 
                "       e.g. " << progname << " --format lmat file1 file2\n" << 
                "\n" <<
                "FILE FORMAT EXAMPLES:\n" <<
                "   Newlines are important.  Comments follow '#'. Spaces separate ascii values.\n" <<
                "   All binary formats use 4 byte ints and floats (unless stated otherwise).\n" <<
                "   (f) = previous number is a float (single precision)\n" << 
                "   (i) = previous number is an integer (signed)\n" <<
                "\n" <<
                "   .ma (a 3x2 matrix in ASCII):\n" <<
                "       20.2 25.4       # <- row1\n" <<
                "       18.1 23.2       # <- row2\n" <<
                "       29.9 22.6       # <- row3\n" <<
                "\n" <<
                "   .mata (a 3x2 matrix in ASCII with dim tags):\n" <<
                "       3 2             # <- 3x2 dim tag (one line in ascii format)\n" <<
                "       20.2 25.4\n" <<
                "       18.1 23.2\n" <<
                "       29.9 22.6\n" <<
                "\n" <<
                "   .mat (a 3x2 matrix in BINARY with dim tags):\n" <<
                "       3(i)2(i)20.2(f)25.4(f)18.1(f)23.2(f)29.9(f)22.6(f)\n" <<
                "\n" <<
                "   .lmata (a 3x2 matrix in ASCII with dim tags and axis labels):\n" <<
                "       3               # <- dim tag for the next line of m axis labels\n" <<
                "       0.0 1.0 2.0     # <- label for the location of the rows\n" <<
                "       2               # <- dim tag for the next line of n axis labels\n" <<
                "       1.5 2.5         # <- label for the location of the cols\n" <<
                "       20.2 25.4       # <- (no dim tags since they were made for axis labels)\n" <<
                "       18.1 23.2\n" <<
                "       29.9 22.6\n" <<
                "\n" <<
                "   .lmat (a 3x2 matrix in BINARY with dim tags and axis labels):\n" <<
                "       3(f)0.0(f)1.0(f)2.0(f)2(i)1.5(f)2.5(f)20.2(f)25.4(f)18.1(f)23.2(f)29.9(f)22.6(f)\n" <<
                ""; // TERMINATE
            //else if (!strcmp(_argv[_i], "long")) {
            //}
    }
    else if (_i < arg_cnt && !strcmp(_argv[_i], "long")) {
        std::cout <<
            "USAGE: " << progname << " [options] <file1> <file2>\n" <<
            "   'Ordered bijective interpolated time warping' aligns matrices along\n" <<
            "   the m axis and by default returns a space separated array of new\n" <<
            "   m axis labels applicable to <file2>.\n" <<
            "\n" <<
            "   Because the default output are m axis labels, formats without axis labels\n" <<
            "   are assumed to be zero indexed by successive integers (i.e., 0, 1, 2, ...).\n" <<
            "\n" <<
            "OPTIONS LEGEND:\n" <<
            "'(*)' [asterik] = the default option;  '<>' = expected parameter\n" <<
            "'[]' = optional parameter;  s = string, f = float, p = percent (0-100)\n" <<
            "\n" <<
            "-h, --help [s]     Generates concise help message\n" <<
            "                   '--help long'    this message\n" <<
            "                   '--help formats' information about file formats\n" <<
            "\n" <<
            "INPUT/OUTPUT:\n" <<
            "--format <s>       The input file format (overriding file extensions)\n" <<
            "                   Otherwise, format is inferred from the file extension:\n" <<
            "                   [dim tag = dimension of the vector or matrix given]\n" <<
            "                   [labeled = vectors specifying position of each row & column]\n" <<
            //"                       .ma     ascii matrix (coming!)\n" <<
            "                       .mat    binary matrix with dim tags\n" <<
            "                       .mata   ascii matrix with dim tags\n" <<
            "                       .lmat   labeled binary matrix with dim tags\n" <<
            "                       .lmata  labeled ascii matrix with dim tags\n" <<
            //"                       .mzXML  ISB mass spectrometry data format (coming!)\n" <<
            //"                       .mzData PSI mass spectrometry data format (coming!)\n" <<
            "                   As argument it is given without the '.' (e.g., --format mata)\n" <<
            "                   \"" << progname << " --help formats\" gives details on the file formats\n" <<
            "                   (Currently) both files should be in the same format\n" <<
            "\n" <<
            "-o, --outfile <s>  Outputs a complete outfile (in same format as input)\n" <<
            "                   The values of the new m axis are still output to STDOUT,\n" <<
            "                   unless the outfile is given as 'STDOUT' (e.g. -o STDOUT),\n" <<
            "                   in which case the file data is piped to STDOUT instead.\n" <<
            "                   Formats WITH axis labels will output a file identical to\n" <<
            "                   <file2> with altered m axis labels.  Formats WITHOUT axis\n" <<
            "                   labels will output a file with data matrix warped.\n" <<
            "\n" <<
            "--images           Creates png images of the alignment process\n" <<
            "                   (if png2matrix installed)\n" <<
            "\n" <<
            "-t, --timefile <s> Given a file of two columns of time points (col1 for file1,\n" <<
            "                   col2 for file2)\n" <<
            "                   emits to STDOUT (after any other output) 4 values on a line:\n" <<
            "                   'ssr asr sad aad'\n" <<
            "                       ssr = sum of the squared residuals\n" <<
            "                       asr = avg of the squared residuals\n" <<
            "                       sad = sum of the absolute values of deltas (time diffs)\n" <<
            "                       aad = avg of the absolute values of deltas\n" <<
            "\n" <<
            "WARPING:\n" <<
            "-s, --score <s>    Score function: (*)cor (Pearson's R), cov (covariance),\n" 
            "                   prd (product), euc (Euclidean distance)\n" <<
            "\n" <<
            "    --nostdnrm     No standard normal (By default the score matrix values are\n" <<
            "                   mean centered and divided by the standard deviation\n" <<
            "                   [unless all values are equal])\n" <<
            "\n" <<
            //"--warp_data        Warps the underlying matrix rather than *the m axis labels.\n" <<
            //"                   NOTE: formats without axis labels default to --warp_data.\n" <<
            //"\n" <<
            "-l, --local        Local rather than *global alignment\n" <<
            "\n" <<
            "-f, --factor <f,f> Local weighting applied to diagonal/gap moves in alignment.\n" <<
            "                   These are given in a comma separated string:\n" <<
            "                       diagonal,gap\n" <<
            "                   (*)Default: (diagonal,gap) = '" << factor_diag << "," << factor_gap << "'\n" <<
            "                   '2,1' gives an unbiased alignment.\n" <<
            "\n" <<
            "-g, --gap <f,f>    Gap penalty given in comma separated string: initiate,extend\n" <<
            "                   (*)Defaults: (gap_init,gap_extend) [by score type]:\n" <<
            "                       'cor' = '" << cor_gap_init << "," << cor_gap_extend << "'\n" <<
            "                       'cov' = '" << cov_gap_init << "," << cov_gap_extend << "'\n" <<
            "                       'prd' = '" << prd_gap_init << "," << prd_gap_extend << "'\n" <<
            "                       'euc' = '" << euc_gap_init << "," << euc_gap_extend << "'\n" <<
            "\n" <<
            "-i, --init <f>     Penalty for initiating alignment (for local alignment only)\n" <<
            "                   (*)Default: " << init_penalty << "\n" <<
            "\n" <<
            "-r, --response <p> Responsiveness of warping.  0 will give a linear warp based\n" << 
            "                   on start and end points.  100 will use all bijective anchors.\n" <<
            "\n" <<
            "MISCELLANEOUS:\n" <<
            "--version          prints version\n" << 
            ""; // TERMINATE

    }
    else {
        std::cout <<
            "USAGE: " << progname << " [options] <file1> <file2>\n" <<
            "   'Ordered bijective interpolated time warping' aligns matrices along\n" <<
            "   the m axis and by default returns a space separated array of new\n" <<
            "   m axis labels applicable to <file2>.\n" <<
            "\n" <<
            "-h, --help [s]     Generates concise help message\n" <<
            "                   '" << progname << " --help long'    thorough help message\n" <<
            "                   '" << progname << " --help formats' info on file formats\n" <<
            "\n" <<
            "INPUT/OUTPUT:\n" <<
            "    --format <s>   The input format (overriding file extensions)\n" <<
            "-o, --outfile <s>  Outputs a complete outfile (in same format as input)\n" <<
            "    --images       Creates png images of the alignment process (if png2matrix)\n" <<
            "-t, --timefile <s> Given a file of two columns of time points,\n" <<
            "                   emits last to STDOUT: 'ssr asr sad aad'\n" <<
            "\n" <<
            "WARPING:\n" <<
            "-s, --score <s>    Score function: (*)cor, cov, prd, euc\n" <<
            "    --nostdnrm     No standard normalization of the score matrix\n" <<
            //"--warp_data        Warps the underlying matrix rather than *the m axis labels.\n" <<
            //"                   NOTE: formats without axis labels default to --warp_data.\n" <<
            //"\n" <<
            "-l, --local        Local rather than (*)global alignment\n" <<
            "-f, --factor <f,f> Local weighting applied to diagonal/gap moves in alignment.\n" <<
            "-g, --gap <f,f>    Gap penalty given in comma separated string: initiate,extend\n" <<
            "-i, --init <f>     Penalty for initiating alignment (for local alignment only)\n" <<
            "-r, --response <p> Responsiveness of warping.  0-100\n" <<
            "\n" <<
            "MISCELLANEOUS:\n" <<
            "--version          prints version\n" << 
            ""; // TERMINATE

    }

}

