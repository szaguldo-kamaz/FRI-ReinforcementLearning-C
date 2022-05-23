/*
** Fuzzy Rule Interpolation-based Reinforcement Learning (ANSI C / AVX version)
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** Author: David Vincze <david.vincze@webcode.hu>
**
** Various contributions by Daniel Palko <palko.daniel@uni-miskolc.hu>
**
*/


#include "frirl_utils.h"
#include "frirl.h"
//#include "FIVE.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <getopt.h>

#include "../inl/fast_abs.inl"
#include "../inl/arr.inl"
#include "../inl/min.inl"


///show rulebase
///\param frirl initialized FRIRL struct
///\return void
void frirl_show_rb(struct frirl_desc *frirl) {

    struct FIVERB *frb = frirl->fiverb;

#ifdef BUILD_MPI
    //if this is a worker process, this would be rubbish
    if (frirl->runmode == FRIRL_MPI && frirl->agent_id != 0) {
        DEBUG_MSG("show_rb: this is a worker process, exit\n");
        return;
    }
#endif

    for (unsigned int i = 0; i < frb->numofrules; i++) {
        printf("%d. ", i);
        for (unsigned int j = 0; j < (frb->rulelength - 1); j++) {
            printf("%.18f ", frb->rant[i * (frb->rulelength - 1) + j]);
//            printf("%.18f ", frb->rseqant[j][i]);

// DEBUG
//            printf("%d ", frb->rant_uindex[i * (frb->rulelength - 1) + j]);
//            printf("%.18f ", frb->rant_veval[i * (frb->rulelength - 1) + j]);
//            printf("%d ", frb->rseqant_uindex[j][i]);
//            printf("%.18f ", frb->rseqant_veval[j][i]);
        }
        printf("Q: %.18f\n", frb->rconc[i]);
    }

}


///show rulebase hexadecimal values
///\param frirl initialized FRIRL struct
///\return 0 if success
void frirl_show_hex_rb(struct frirl_desc *frirl) {

    struct FIVERB *frb = frirl->fiverb;

#ifdef BUILD_MPI
    //if this is a worker process, this would be rubbish
    if (frirl->runmode == FRIRL_MPI && frirl->agent_id != 0) {
        DEBUG_MSG("show_hex_rb: this is a worker process, exit\n");
        return;
    }
#endif

    for (unsigned int i = 0; i < frb->numofrules; i++) {
        for (unsigned int j = 0; j < (frb->rulelength - 1); j++) {
            printf("%.13a ", frb->rant[i * (frb->rulelength - 1) + j]);
            /*        printf("%x",(char )frb->rb[i*(frirl->frb->rulelength)+j]);
                    printf("%x",(char )frb->rb[i*(frirl->frb->rulelength)+j+1]);
                    printf("%x",(char )frb->rb[i*(frirl->frb->rulelength)+j+2]);
                    printf("%x",(char )frb->rb[i*(frirl->frb->rulelength)+j+3]);

                    printf("%x",(char )frb->rb[i*(frirl->frb->rulelength)+j+4]);
                    printf("%x",(char )frb->rb[i*(frirl->frb->rulelength)+j+5]);
                    printf("%x",(char )frb->rb[i*(frirl->frb->rulelength)+j+6]);
                    printf("%x ",(char )frb->rb[i*(frirl->frb->rulelength)+j+7]);*/
        }
        printf("%.13a\n",frb->rconc[i]);
    }

}


///save rulebase to test file
///\param frirl initialized FRIRL struct
///\param file_name
///\return 0 if success
int frirl_save_rb_to_text_file(struct frirl_desc *frirl, const char *file_name) {

    struct FIVERB *frb = frirl->fiverb;

#ifdef BUILD_MPI
    //if this is a worker process, this would be rubbish
    if (frirl->runmode == FRIRL_MPI && frirl->agent_id != 0) {
        DEBUG_MSG("save_rb_to_text_file: this is a worker process, exit\n");
        return 0;
    }
#endif

    int rbfd, towrite, written;
    unsigned char strbuf[32];

    rbfd = open(file_name, O_WRONLY | O_CREAT | O_TRUNC, 00644);
    if (rbfd == -1) {
        perror("FRIRL_save_rb_to_file_text: open(): ");
        return -1;
    }

    for (unsigned int i = 0; i < frb->numofrules; i++) {

        for (unsigned int j = 0; j < (frb->rulelength - 1); j++) {
            towrite = sprintf(strbuf, "%.18f ", frb->rant[i * (frb->rulelength - 1) + j]);
            written = write(rbfd, strbuf, towrite);
            if (written == -1) {
                perror("FRIRL_save_rb_to_file_text: write(): ");
                return -1;
            }
        }

        towrite = sprintf(strbuf, "%.18f \n", frb->rconc[i]);
        written = write(rbfd, strbuf, towrite);
        if (written == -1) {
            perror("FRIRL_save_rb_to_file_text: write(rconc): ");
            return -1;
        }

    }

    close(rbfd);

    return 0;
}


///save rulebase to binary file
///\param frirl initialized FRIRL struct
///\param file_name
///\return 0 if success
int frirl_save_rb_to_bin_file(struct frirl_desc *frirl, const char *file_name) {

    struct FIVERB *frb = frirl->fiverb;

#ifdef BUILD_MPI
    //if this is a worker process, this would be rubbish
    if (frirl->runmode == FRIRL_MPI && frirl->agent_id != 0) {
        DEBUG_MSG("save_rb_to_bin_file: this is a worker process, exit\n");
        return 0;
    }
#endif

    int rbfd;
    int written;
    int totalbytes;

    unsigned char strbuf[32];

    rbfd = open(file_name, O_WRONLY | O_CREAT | O_TRUNC, 00644);
    if (rbfd == -1) {
        perror("frirl_dump_rb_to_file: open(): ");
        return -1;
    }

    written = write(rbfd, &frb->numofrules, sizeof(int)); // number of rules
    if (written == -1) {
        perror("frirl_dump_rb_to_file: write(): ");
        return -1;
    }

    //TODO rant+rconc
//    written = write(rbfd, frb->rb, sizeof(fri_float)*frb->numofrules * frb->rulelength);
//    if (written == -1) {
//        perror("frirl_dump_rb_to_file: write(): ");
//        return -1;
//    }
    for (unsigned int r = 0; r < frb->numofrules; r++) {
        written = write(rbfd, frb->rant + r * frirl->numofantecedents, sizeof(double) * frirl->numofantecedents);
        if (written == -1) {
            perror("frirl_dump_rb_to_file: write(): ");
            return -1;
        }

        totalbytes += written;

        written = write(rbfd, frb->rconc + r, sizeof(double));
        if (written == -1) {
            perror("frirl_dump_rb_to_file: write(): ");
            return -1;
        }

        totalbytes += written;
    }

    close(rbfd);

//    return written;
//
//    int rbfd, written;
//
//    rbfd = open(file_name, O_WRONLY | O_CREAT | O_TRUNC, 00644);
//    if (rbfd == -1) {
//        perror("frirl_dump_rb_to_file: open(): ");
//        return -1;
//    }
//    written = write(rbfd, &frb->numofrules, sizeof(int)); // number of rules
//    if (written == -1) {
//        perror("frirl_dump_rb_to_file: write(): ");
//        return -1;
//    }
//    //TODO rant+rconc
//    written = write(rbfd, frb->rb, sizeof(fri_float) * frb->numofrules * frb->rulelength);
//    if (written == -1) {
//        perror("frirl_dump_rb_to_file: write(): ");
//        return -1;
//    }
//    close(rbfd);
//
//    return written;
}


///load rulebase from binary file
///\param frirl initialized FRIRL struct
///\param file_name
///\return 0 if success
int frirl_load_rb_from_bin_file(struct frirl_desc *frirl, const char *file_name) {

    struct FIVERB *frb = frirl->fiverb;

    int rbfd;
    int readbytes;
    int totalbytes = 0;
    int numofrules;
    double rule[frb->rulelength];

    rbfd = open(file_name, O_RDONLY);
    if (rbfd == -1) {
        perror("frirl_load_rb_from_file_bin: open(): ");
        return -1;
    }

    readbytes = read(rbfd, &numofrules, sizeof(int));  // number of rules
    if (readbytes == -1) {
        perror("frirl_load_rb_from_file_bin: read(): ");
        return -1;
    }

//    frirl->rules_len = 0;
    frb->numofrules = 0;
    frb->newrant = frb->rant;
    frb->newrconc = frb->rconc;
    frirl->fus_is_rule_inserted = 0;

    for (int r = 0; r < numofrules; r++) {
        readbytes = read(rbfd, rule, frb->rulelength * sizeof(double));
        if (readbytes == -1) {
            perror("frirl_load_rb_from_file_bin: read(): ");
            return -1;
        }
        totalbytes += readbytes;

        five_add_rule(frirl->fiverb, rule);
//    frirl->rules_len++; //todo frirl.rules_len miert is kell?

    }
//    frirl->rules_len = numofrules;

    close(rbfd);

    return readbytes;

//    int rbfd, readbytes;
//
//    rbfd = open(file_name, O_RDONLY);
//    if (rbfd == -1) {
//        perror("frirl_load_rb_from_file_bin: open(): ");
//        return -1;
//    }
//    readbytes = read(rbfd, &frb->numofrules, sizeof(int)); // number of rules
//    if (readbytes == -1) {
//        perror("frirl_load_rb_from_file_bin: read(): ");
//        return -1;
//    }
//    frirl->rules_len = frb->numofrules;
//    readbytes = read(rbfd, frb->rb, sizeof(fri_float)*frb->numofrules * frb->rulelength);
//    if (readbytes == -1) {
//        perror("frirl_load_rb_from_file_bin: read(): ");
//        return -1;
//    }
//    close(rbfd);
//
//    // Generate "rant"
//    int l = 0;
//    for (unsigned int j = 0; j < frb->numofrules; j++) {
//        for (unsigned int k = 0; k < frb->rulelength - 1; k++) {
//
//            frb->rseqant[k][j] = frb->rb[j * frb->rulelength + k];;
//            frb->rseqant_uindex[k][j] = get_vag_abs_min_i_fixres(
//                frb->u + (frb->univlength) * k,
//                frb->univlength, frb->rseqant[k][j],
//                frb->udivs[k]);
//
//            frb->rant[l] = frb->rb[j * frb->rulelength + k];;
//            frb->rant_uindex[l++] = frb->rseqant_uindex[k][j];
//            DEBUG_MSG("frirl->frb->rb[%d]: %f\n", j * frb->rulelength + k, frb->rb[j * frb->rulelength + k]);
//        }
//        frb->rconc[j] = frb->rb[j * frb->rulelength + frb->rulelength - 1];
//    }
////TODO generate rconc
//
//    return readbytes;
}


void frirl_print_usage() {

    printf("FRIRL learning usage:\n\n");
    printf("    -m --runmode <runmode>\n");
    printf("\tMode of operation. Possible values: seq, omp, mpi, test. With 'test', a rulebase file should also be specified.\n\tDefault: seq\n\n");
    printf("    -f --rbfile <rule-base file path>\n");
    printf("\tPath to binary rule-base file.\n\n");
    printf("    -r --reductionstrategy <rule-base reduction strategy>\n");
    printf("\tRule-base reduction strategy. Possible values: noreduce, default. Incompatible with 'test' mode.\n\tDefault: noreduce\n\n");
    printf("    -d --draw\n");
    printf("\tEnable visualization.\n\n");
    printf("    -q --quiet\n");
    printf("\tEnable quiet mode.\n\n");
//    printf("\t-c construct rulebase"); // TODO

}


void frirl_parse_cmdline(struct frirl_desc *frirl, int argc, char **argv) {

    frirl->argc = argc;
    frirl->argv = argv;

    if (argc == 1)
        return;

//    enum frirl_runmode rmode = frirl->runmode;
    char rmode[81];
    *rmode = 0;

    char reduce[81];
    *reduce = 0;

    static struct option cmdline_options[] = {
      {"runmode",           required_argument, 0, 'm'},
      {"rbfile",            required_argument, 0, 'f'},
      {"reductionstrategy", required_argument, 0, 'r'},
      {"draw",                    no_argument, 0, 'd'},
      {"quiet",                   no_argument, 0, 'q'},
//      {0, 0, 0, 0}
    };

    int option_index = 0;
    int opt;
    while ((opt = getopt_long(argc, argv, "m:f:r:dq", cmdline_options, &option_index)) > 0) {
        switch (opt) {
        case '?':
            frirl_print_usage();
            exit(-1);
            break;
        case 'm':
            strcpy(rmode, optarg);
            break;
        case 'f':
            frirl->rbfile = malloc(strlen(optarg) * sizeof(char));
            strcpy(frirl->rbfile, optarg);
            break;
        case 'r':
            strcpy(reduce, optarg);
            break;
        case 'd':
            frirl->visualization = 1;
            break;
        case 'q':
            frirl->verbose = 0;
            break;
        }
    }

    if (*rmode) {
        if (!strcmp(rmode, "seq")) {
            frirl->runmode = FRIRL_SEQ;
        } else if (!strcmp(rmode, "omp")) {
            frirl->runmode = FRIRL_OMP;
        } else if (!strcmp(rmode, "mpi")) {
            frirl->runmode = FRIRL_MPI;
        } else if (!strcmp(rmode, "test")) {
            frirl->runmode = FRIRL_TEST;
        } else {
            printf("Invalid runmode!\n\n");
            frirl_print_usage();
            exit(-1);
        }
    }

    if (*reduce) {
        if (!strcmp(reduce, "noreduce")) {
            frirl->reduce_rb = 0;
            frirl->reduction_strategy = FRIRL_REDUCTION_STRATEGY_NOREDUCE;
        } else if (!strcmp(reduce, "default")) {
            frirl->reduce_rb = 1;
            frirl->reduction_strategy = FRIRL_REDUCTION_STRATEGY_DEFAULT;
        } else {
            printf("Invalid reduction strategy!\n\n");
            frirl_print_usage();
            exit(-1);
        }
    }

}
