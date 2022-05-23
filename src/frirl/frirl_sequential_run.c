/*
** Fuzzy Rule Interpolation-based Reinforcement Learning (ANSI C / AVX version)
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** Author: David Vincze <david.vincze@webcode.hu>
**
** Various contributions by Daniel Palko <palko.daniel@uni-miskolc.hu>
**
** Reduction strategy 1+2 ported from MATLAB by Alex Martossy
*/


#include "frirl.h"
#include <string.h>
#include <stdio.h>
#include <math.h>


///FRIRL main loop (sequential version, without MPI and openmp)
///\param frirl initialized frirl struct
///\param verbose 0:disable, 1:enable; print some messages, like actual episode num, reward, etc.
///\return void
void frirl_sequential_run(struct frirl_desc *frirl) {

    struct FIVERB *frb = frirl->fiverb;
    // do not modify
    double prev_reward;
//    double *prev_rb;
    double *prev_rconc;
    double *tmp_rconc;
    int prev_numru;
    int prev_steps;
    int epend = 0;
    int redend = 0;
    frirl->epended = 0;

//disabled temporary
#ifdef BUILD_AVX2off
    double  qdiff_final_tolerance_arr[8] __attribute__ ((aligned (32))) = { frirl->qdiff_final_tolerance, frirl->qdiff_final_tolerance, frirl->qdiff_final_tolerance, frirl->qdiff_final_tolerance };
    double *qdiff_final_tolerance_p = qdiff_final_tolerance_arr;
    double  absmask[4] __attribute__ ((aligned (32))) = { -0.0, -0.0, -0.0, -0.0 };
    double *absmask_p = absmask;
//    double *prev_rconc_p = prev_rconc;
//    double *curr_rconc_p = frb->rconc;
#endif

    // local
//    prev_rb = (double *) MALLOC(sizeof(double) * frb->rulelength * FRIRLRULES);
    prev_rconc = (double *) MALLOC(sizeof(double) * frirl->five_maxnumofrules);
    int epchunk = 1;
    int maxep = (frirl->runmode == FRIRL_MPI || frirl->runmode == FRIRL_OMP) ? FRIRL_AGENT_EPCHUNK : frirl->max_episodes;

    // simulation starts here
    if (frirl->construct_rb == 1) {

        for (;;) {

            if (!(epchunk < maxep)) {
                if (!(frirl->episode_num < frirl->max_episodes)) {
                    frirl->is_running = 0;
                }
//                printf("%d exit2\n", frirl->agent_id);
                break;
            }

            // save values from the previous iteration
            prev_numru  = frb->numofrules;
            prev_reward = frirl->reward.ep_total_value;
            prev_steps  = frirl->reward.ep_total_steps;
//            memcpy(prev_rb, frb->rb, sizeof(double) * frb->rulelength * FRIRLRULES);
            memcpy(prev_rconc, frb->rconc, sizeof(double) * frirl->five_maxnumofrules);

            // run the FRIRL iteration
            frirl_episode(frirl);

            printf("#%d Episode: %d\tSteps: %d\tReward: %s%f%s\tRules: %d\n",
                frirl->agent_id, frirl->episode_num, frirl->reward.ep_total_steps,
                frirl->reward.ep_total_value > frirl->reward_good_above ? TERM_GREEN : TERM_RED,
                frirl->reward.ep_total_value, TERM_NC, frb->numofrules);

            // did we find the final rule-base?
            if ( ( (prev_numru == frb->numofrules) &&
                   (prev_steps == frirl->reward.ep_total_steps) &&
                   (frirl->reward.ep_total_value > frirl->reward_good_above) &&
                   (prev_reward == frirl->reward.ep_total_value)
                 ) || (frirl->user_exited == 1) ) {

                epend = 1;
                frirl->epended = 1;
                if (frirl->verbose > 0) {
                    printf("Rule-base size and reward are the same as in the previous iteration.\n");
                }
                // cartpole comes here 2x times
                // acrobot  comes here 8x times

// disabled temporary, not much effect here
#ifdef BUILD_AVX2off
                // it's OK to go beyond numofrules when comparing, because prev_rconc is a full copy (FRIRLRULES) of rconc, so trash is also the same at the end
                asm volatile (
                    "xor    %%rcx, %%rcx       \n" // i = 0
                    "mov    %0, %%rax          \n" // &prev_rconc[0]
                    "mov    %1, %%rbx          \n" // &frb->rconc[0]
                    "mov    %2, %%rdx          \n" // &qdiff_final_tolerance_arr[0]
                    "vmovapd (%%rdx), %%ymm2   \n" // ymm2 = qdiff_final_tolerance_arr[]
                    "mov    %5, %%rdx          \n" // &absmask[0]
                    "vmovapd (%%rdx), %%ymm4   \n" // ymm4 = -0.0 -0.0 -0.0 -0.0 (maszk miatt)
                    "mov    %4, %%rdx          \n" // rdx = frb->avx2_rbsize | arrsize = numofrules*8/32 -> numofrules/4
                "l:  \n"
                    "vmovapd (%%rax), %%ymm0             \n" // ymm0 = prev_rconc[]
                    "vmovapd (%%rbx), %%ymm1             \n" // ymm1 = frb->rconc[]
                    "vsubpd %%ymm0, %%ymm1, %%ymm3       \n" // ymm3 = ymm1 - ymm0 (frb->rconc - prev_rconc)
//                    "vsubpd %%ymm1, %%ymm0, %%ymm4       \n" // ymm4 = ymm0 - ymm1 (prev_rconc - frb->rconc)
//                    "vmaxpd %%ymm3, %%ymm4, %%ymm1       \n" // ymm1 = max(ymm3,ymm4) - we got abs()
                    "vandnpd %%ymm3, %%ymm4, %%ymm1      \n" // ymm1 = ymm3 & absmask - we got abs()
                    "vcmppd $13, %%ymm2, %%ymm1, %%ymm0  \n" // ymm0 = _CMP_GE_OS ymm1 > ymm2
                    "vptest %%ymm0, %%ymm0               \n" // csak flageket valtoztat (hogy az egesz ymm0 nulla-e)
                    "jnz    ski                          \n" // ugrik ha valahol igaz, szoval ymm1 > ymm2 (szoval valamelyik nagyobb mint a tolerancia)
                    "inc    %%rcx              \n"
                    "add    $32, %%rax         \n"
                    "add    $32, %%rbx         \n"
                    "cmp    %%rdx,%%rcx        \n"
                    "jne    l                  \n"
                    "jmp    out\n"
                "ski: \n"
                    "movb   $0, %3             \n" // epend = 0
                "out: \n"

                : : "m" (prev_rconc), "m" (frb->rconc), "m" (qdiff_final_tolerance_p), "m" (epend), "m" (frb->avx2_rbsize), "m" (absmask_p) : "rax", "rbx", "rcx", "rdx", "ymm0", "ymm1", "ymm2", "ymm3", "ymm4", "memory" );
                //clobber list output value m ?

#else // no avx

                for (int i = 0; i < frb->numofrules; i++) {

                    fri_float r1 = prev_rconc[i];
                    fri_float r2 = frb->rconc[i];

                    if (fabs(r2 - r1) >= frirl->qdiff_final_tolerance) {
                        epend = 0;
                        if (frirl->verbose > 0) {
                            printf("Greater at rule %d. %.18f - %.18f = %.18f (max: %.18f)\n", i, frb->rconc[i], prev_rconc[i], frb->rconc[i] - prev_rconc[i], frirl->qdiff_final_tolerance);
                            // no break here, because we would like to see all rules with Qs above tolerance
                        } else {
                            break; // one rule above tolerance is enough
                        }
                    }
                }
#endif
            }

            if ( (epend == 1) || (frirl->user_exited == 1) ) {
                frirl->is_running = 0;
                printf("-----------------------------------------------------------------\n");
                printf("No more significant changes in rule-base. RB considered complete.\n");
                printf("-----------------------------------------------------------------\n");
                break;
            }

            frirl->episode_num++;
            epchunk++;

        }

    }  // construct end

    // ***************************************************************************************************************

    // reduction of the constructed rule-base
    if (frirl->reduce_rb == 1) {

        tmp_rconc = (double *) MALLOC(sizeof(double) * frirl->five_maxnumofrules);

        frirl->original_learning = 1;

        //make copy of rconc, this will be used to mark rules as important, by setting a Q to -nan, while leaving original frb->rconc untouched
        memcpy(tmp_rconc, frb->rconc, sizeof(double) * frirl->five_maxnumofrules);

        int steps_frirl_incremental;
        int numofrules = frb->numofrules;
        frirl->reduction_state = 1;
        int iterations;
        unsigned int mindex;
        double mvalue;
        double *removed_rule = (double *) MALLOC(sizeof(double) * frb->rulelength);
        double diff_rewardf;


        if ( (frirl->reduction_strategy == 1) || (frirl->reduction_strategy == 2) ) {
            iterations = (numofrules + 1);
        } else {
            iterations = 10000;  //used for frirl->reduction_strategy == 3 (not implemented yet)
        }

        frirl_episode(frirl);  //run episode to determine steps_frirl_incremental
        steps_frirl_incremental = frirl->reward.ep_total_steps;

        for (frirl->episode_num = 1; frirl->episode_num <= iterations; frirl->episode_num++) {

            // run the FRIRL iteration
            frirl_episode(frirl);

            printf("Reduction Episode: %d\tSteps: %d\tReward: %f\tEpsilon: %f\tRules: %d\n",
                frirl->episode_num, frirl->reward.ep_total_steps,
                frirl->reward.ep_total_value, frirl->epsilon, frb->numofrules);

            if ( (frirl->reduction_strategy == 1) || (frirl->reduction_strategy == 2) ) {

                if (frirl->episode_num > 1) {

                    diff_rewardf = prev_reward - frirl->reward.ep_total_value;

                    if ( (frirl->reward.ep_total_value > frirl->reward_good_above) &&
                         (frirl->reward.ep_total_steps == steps_frirl_incremental) &&
                         (fabs(diff_rewardf) <= frirl->reduction_reward_tolerance) ) {

                        // omission of this rule could be a good idea
                        diff_rewardf = prev_reward - frirl->reward.ep_total_value;  // redundant?
                        if (frirl->verbose > 0) {
                            printf("Reduction Episode: %d\tEliminated rule: no: %d. - ", frirl->episode_num, mindex + 1);
                            for (unsigned int j = 0; j < frb->rulelength; j++) {
                                printf(" %f", removed_rule[j]);
                            }
                            printf(" \tReward diff was: %f\n", diff_rewardf);
                        }
                        prev_reward = frirl->reward.ep_total_value;

                    } else {

                        // omission of the rule was a bad idea
                        // restore rule base
                        memcpy(tmp_rconc, prev_rconc, sizeof(double) * frirl->five_maxnumofrules);
                        tmp_rconc[mindex] = 0.0 / 0.0;  // mark rule as improtant (set conclusion to -nan)

                        if (frirl_load_rb_from_bin_file(frirl, "reduction_tmp.frirlrb.bin") < 0) {
                            printf("Error while loading the binary rule-base file!\n");
                            exit(-1);
                        }
                        frb = frirl->fiverb;

                        if (frirl->verbose > 0) {
                            printf("Reduction Episode: %d\tRule stays: no: %d. - ", frirl->episode_num, mindex + 1);
                            for (unsigned int j = 0; j < (frb->rulelength - 1); j++) {
                                printf(" %f", frb->rseqant[j][mindex]);
                            }
                            printf(" %f\n", frb->rconc[mindex]);
                        }
                    }

                } else {
                    prev_reward = frirl->reward.ep_total_value;
                }

            }  // if frirl->reduction_strategy == 1 / 2

            // 1. min q
            if (frirl->reduction_strategy == 1) {

                // find min Q
                mvalue = fabs(tmp_rconc[0]);
                mindex = 0;
                for (int j = 1; j < frb->numofrules; j++) {
                    // if next element is less than MIN, OR is a number (not nan)
                    if(mvalue > fabs(tmp_rconc[j]) || (mvalue != mvalue && tmp_rconc[j] == tmp_rconc[j]) ) {  // TODO why these conds?
                        mvalue = fabs(tmp_rconc[j]);
                        mindex = j;
                    }
                }

                // take out next candidate rule
                if (mvalue!=mvalue) {
                    if (frirl->verbose > 0) {
                        printf("Smallest rulebase found. Exiting.\n");
                    }
                    redend = 1;
                } else {
                    // make copy of frb
                    frirl_save_rb_to_bin_file(frirl, "reduction_tmp.frirlrb.bin");
                    for (unsigned int k = 0; k < (frb->rulelength - 1); k++) {
                        removed_rule[k] = frb->rseqant[k][mindex];
                    }
                    removed_rule[frb->rulelength - 1] = tmp_rconc[mindex];

                    // make copy of marked as important conclusions, to be able to restore if needed
                    memcpy(prev_rconc, tmp_rconc, sizeof(double) * frirl->five_maxnumofrules);

                    // remove conclusion, from marked as important conclusions
                    for(int k = mindex; k < (frb->numofrules - 1); k++) {
                        tmp_rconc[k] = tmp_rconc[k + 1];
                    }

                    // remove rule from frb
                    five_remove_rule(frb, mindex);

                    // now tmp_rconc is the same as frb->rconc; some rules may be marked as important, but both have the same rules
                }
            }

            // 2. max q
            if (frirl->reduction_strategy == 2) {

                // find max Q
                mvalue = fabs(tmp_rconc[0]);
                mindex = 0;
                for (int j = 1; j < frb->numofrules; j++) {
                    // if next element is greater than MAX, OR is a number (not nan)
                    if (mvalue < fabs(tmp_rconc[j]) || (mvalue != mvalue && tmp_rconc[j] == tmp_rconc[j]) ) { // TODO what are these conds?
                        mvalue = fabs(tmp_rconc[j]);
                        mindex = j;
                    }
                }

                // take out next candidate rule
                if (mvalue != mvalue) {
                    if (frirl->verbose > 0) {
                        printf("Smallest rulebase found. Exiting.\n");
                    }
                    redend = 1;
                } else {
                    // make a copy of frb
                    frirl_save_rb_to_bin_file(frirl, "reduction_tmp.frirlrb.bin");

                    for (unsigned int k = 0; k < (frb->rulelength - 1); k++) {
                        removed_rule[k] = frb->rseqant[k][mindex];
                    }
                    removed_rule[frb->rulelength - 1] = tmp_rconc[mindex];

                    // make copy of marked as important conclusions, to be able to restore if needed
                    memcpy(prev_rconc, tmp_rconc, sizeof(double) * frirl->five_maxnumofrules);

                    // remove conclusion, from marked as important conclusions
                    for(int k = mindex; k < (frb->numofrules - 1); k++) {
                        tmp_rconc[k] = tmp_rconc[k +1 ];
                    }

                    // remove rule from frb
                    five_remove_rule(frb, mindex);
                    // now tmp_rconc is equal to frb->rconc; some rules may be marked as important, but both have the same rules
                }
            }

            if (redend == 1) {
                break;
            }

        }

        free(tmp_rconc);

    }  // reduce end

//    free(prev_rb);
    free(prev_rconc);

}
