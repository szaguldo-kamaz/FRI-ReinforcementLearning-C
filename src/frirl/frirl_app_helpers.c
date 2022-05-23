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


#include "frirl.h"
#include "frirl_app_helpers.h"
#include "frirl_test.h"
#include <string.h>

#include <stdio.h>


#ifdef BUILD_VISUALIZATION
///empty inline function performance trick (visualization runtime disabled)
static inline void draw_func(struct frirl_desc *frirl, fri_float *new_curr_state, double action, unsigned int steps) {}
#endif //BUILD_VISUALIZATION


///generate linear values
///\param arr array of values
///\param len lenght of the array
///\param div step between the values
///\return void
void frirl_gen_fixres_arr(fri_float *arr, int len, fri_float div) {

    int i;
    fri_float from = -((len - 1) * div) / 2;

    for (i = 0; i < len / 2 + 1; i++) {
        arr[i] = from + div * i;
    }
    for (; i < len; i++) {
        arr[i] = arr[len - 1 - i] * -1;
    }

}

///FRIRL main loop
///\param frirl initialized frirl struct
///\param verbose 0:disable, 1:enable; print some messages, like actual episode num, reward, etc.
///\return void
void frirl_run(struct frirl_desc *frirl, int verbose) {

    // if a rule-base file was supplied, then use it
    if (frirl->rbfile != NULL) {
        if (frirl_load_rb_from_bin_file(frirl, frirl->rbfile) < 0) {
            printf("Error while loading the binary rule-base file: %s!\n", frirl->rbfile);
            exit(-1);
        } else {
            printf("Loaded %d rules from binary rule-base file: %s\n", frirl->fiverb->numofrules, frirl->rbfile);
            if (frirl->verbose > 1) {
                frirl_show_rb(frirl);
            }
        }
    } else {
        if (frirl->construct_rb == 0 && frirl->reduce_rb == 1) {
            printf("Incrementally constructed rule-base file is missing, \nplease run the construction process first, \nthen supply the constructed rule-base file! \n(example -f example.frirlrb.bin)\n");
            exit(-1);
        }
    }

#ifdef BUILD_VISUALIZATION
    if (frirl->visualization) {
        frirl_visualization_init(frirl);
    } else {
        frirl->draw_func = draw_func;
    }
#endif //BUILD_VISUALIZATION

    switch (frirl->runmode) {
        case FRIRL_SEQ:  frirl_sequential_run(frirl); break;
        case FRIRL_OMP:  frirl_omp_run(frirl); break;
        case FRIRL_MPI:  frirl_mpi_run(frirl); break;
        case FRIRL_TEST: frirl_test_run(frirl); break;
    }

#ifdef BUILD_VISUALIZATION
    if (frirl->visualization) {
        frirl_visualization_deinit();
    }
#endif //BUILD_VISUALIZATION

//    //run test? // TODO rendes kapcsolokkal megcsinalni az egeszet
///*    if (frirl->argc >= 3) {
//        if (!strcmp(frirl->argv[1], "-t")) {
//            frirl->reduction_state = 1; // this way there will be no updates in the RB (originally for testing truncated RBs in reduction state)
//
//            frirl_load_rb_from_bin_file(frirl, frirl->argv[2]);
//
//            frirl_test_run(frirl, 1);
//        }
//    } else {*/
//#if defined BUILD_MPIoff
//        //todo
//#elif defined BUILD_OPENMP
//        frirl_omp_run(frirl, verbose);
//#else
//        frirl_sequential_run(frirl, verbose);
//#endif
////    }
}


void frirl_visualization_init(struct frirl_desc *frirl) {

#ifdef BUILD_VISUALIZATION
    char title[81];

    sprintf(title, "FRI Reinforcement Learning: %s", frirl->argv[0]);
    initGLUT(&frirl->argc, frirl->argv, frirl->gui_width, frirl->gui_height, title);
#endif //BUILD_VISUALIZATION

}


void frirl_visualization_deinit() {

#ifdef BUILD_VISUALIZATION
    // Deinit visualization
//    deinitGLUT();
#endif

}
