
#ifndef _CMDPARSER_H
#define _CMDPARSER_H


class CmdParser {
    public:
        // internal
        int _i;
        char ** _argv;

        // options:
        char * format;
        char * outfile;
        bool images;
        char * timefile;
        char * score;
        bool nostdnrm;
        bool local;
        float factor_diag;
        float factor_gap;
        float gap_init;
        float gap_extend;
        float init_penalty;
        float response;
        void help_func(int arg_cnt);
        void print_version(char * version);
        //bool warp_data;
        char * smat_out;
        char * smat_in;

        char * infiles[2];

        char * progname;
        
        // initialize
        CmdParser(int argc, char ** argv, char * version);

        // internal:
        void set_flag(bool &val);
        void set_string(char ** string_ptr);
        void set_float(float &val);
        void set_comma_list_float(float &val1, float &val2);
        void verify_infiles();
        void get_format_from_file1();
        bool file_is_readable(char * filename);
        void set_defaults_by_score(char * score_arg);
        bool eq(char * first, char * sec);
//XCMS***********
        void print_argv();
//END XCMS***********
};

#endif
