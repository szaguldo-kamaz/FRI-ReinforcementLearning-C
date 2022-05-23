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


#include "frirl_test.h"
#include "frirl_app_helpers.h"
#include "frirl.h"
#include <stdio.h>
#include <string.h>


void frirl_test_run(struct frirl_desc *frirl) {

//    if (frirl->rbfile == 0 || strlen(frirl->rbfile) == 0) {
//        printf("You need to specify the binary rule-base file!\n\n");
//        frirl_print_usage();
//        exit(-1);
//    }

    if (frirl->reduce_rb) {
        printf("The 'test' mode does not perform reduction. Parameter omitted.\n\n");
//TODO: maybe a warning is enough, no need to exit?
        frirl_print_usage();
        exit(-1);
    }

//        if (argc > 2) {
//        frirl.construct_rb = atoi(argv[2]);
//        if (frirl.construct_rb != 1) {
//            frirl.construct_rb = 0;
//            frirl.reduction_state = 0;
//        }
//        }
//        if (argc > 3) {
//        frirl.reduce_rb = 1;
//        frirl.reduction_strategy = atoi(argv[3]);
//        if ( (frirl.reduction_strategy < 1) || (frirl.reduction_strategy > 2) ) {
//            frirl.reduction_strategy = 1;
//        }
//        }
//    }
//    printf("Construct RB: %d\n", frirl.construct_rb);
//    printf("Reduce RB: %d Strategy: %d\n", frirl.reduce_rb, frirl.reduction_strategy);

//    frirl->reduction_state = 0;

//    if (frirl_load_rb_from_bin_file(frirl, frirl->rbfile) < 0) {
//        printf("Error while loading the binary rule-base file: %s!\n", frirl->rbfile);
//        exit(-1);
//    }
//
//    printf("Loaded %d rules from binary rule-base file: %s\n", frirl->fiverb->numofrules, frirl->rbfile);
//    if (frirl->verbose > 1) {
//        frirl_show_rb(frirl);
//    }

    frirl->construct_rb = 0;
    frirl->reduction_state = 1;

    // run the FRIRL iteration
    frirl_episode(frirl);

    if (frirl->verbose != 0) {
        printf("Steps:\t%d\nRules:\t%d\nReward:\t%f\n",
            frirl->reward.ep_total_steps, frirl->fiverb->numofrules,
            frirl->reward.ep_total_value);
    }

    // did we find the final rule-base?
    if (frirl->reward.ep_total_value > frirl->reward_good_above) {
        printf("\nSuccess!\n");
        exit(0);  // TODO exit otherwise RB will be dumped
    } else {
        printf("\nInvalid!\n");
        exit(-1);  // TODO exit otherwise RB will be dumped
    }

}
