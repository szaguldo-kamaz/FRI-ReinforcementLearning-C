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
//#include "FIVE.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "../inl/fast_abs.inl"
#include "../inl/arr.inl"
#include "../inl/min.inl"


///initialize the given frirl_desc struct based on the struct values
///\param frirl used as configuration params
///\return 0 if success
int frirl_init(struct frirl_desc *frirl)
{
//    struct FIVERB *tmpfrb;
    struct timeval t1;
    int i, j, k, knu, st;
    int realnumofrules;

    //FIVE compatibility, todo
    int usize = frirl->statedims[0].universe_len;
    fri_float *u = (fri_float *)MALLOC((usize * (frirl->statedims_len + 1)) * sizeof(fri_float));
    if (u == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(u)\n");
#endif
        return -1;
    }
    for (i = 0; i < frirl->statedims_len; i++) {
        memcpy(u + i * usize, frirl->statedims[i].universe, usize * sizeof(fri_float));
    }
    memcpy(u + frirl->statedims_len * usize, frirl->actiondim.universe, usize * sizeof(fri_float));

//    double *u = frirl->statedim_arr[0].universe;
//    int usize = frirl->statedim_arr[0].universe_len;

    frirl->numofantecedents = frirl->statedims_len + 1;

    frirl->possible_states = MALLOC(sizeof(struct frirl_values_desc) * frirl->statedims_len);
    if (frirl->possible_states == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(possible_states)\n");
#endif
        return -1;
    }
    for (i = 0; i < frirl->statedims_len; i++) {
        frirl->possible_states[i].values_len = frirl->statedims[i].values_len;
        frirl->possible_states[i].values = MALLOC(sizeof(double) * frirl->statedims[i].values_len);
        if (frirl->possible_states[i].values == NULL) {
#ifdef DEBUG
            printf("frirl_init: Error - malloc(frirl->possible_states[%d].values)\n",i);
#endif
            return -1;
        }
        for (j = 0; j < frirl->statedims[i].values_len; j++) {
            frirl->possible_states[i].values[j] = frirl->statedims[i].values[j];
        }
        frirl->possible_states[i].epsilon = frirl->statedims[i].values_div;
    }

    frirl->possible_actions = MALLOC(sizeof(struct frirl_values_desc));
    if (frirl->possible_actions == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(frirl->possible_actions)\n");
#endif
        return -1;
    }
    frirl->possible_actions->values_len = frirl->actiondim.values_len;
    frirl->possible_actions->values = MALLOC(sizeof(fri_float)*frirl->actiondim.values_len);
    if (frirl->possible_actions->values == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(frirl->possible_actions->values)\n");
#endif
        return -1;
    }
    frirl->possible_actions->vevalues = MALLOC(sizeof(fri_float)*frirl->actiondim.values_len);
    if (frirl->possible_actions->vevalues == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(frirl->possible_actions->vevalues)\n");
#endif
        return -1;
    }
    for (j = 0; j < frirl->actiondim.values_len; j++) {
        frirl->possible_actions->values[j] = frirl->actiondim.values[j];
    }
    frirl->possible_actions->epsilon = frirl->actiondim.values_div;

//    frirl->curr_state = MALLOC(sizeof(double)*frirl->statedim_arr_len);
//    frirl->quant_obs = MALLOC(sizeof(double)*frirl->statedim_arr_len);

    // FIVE init part
//    frirl->fiverb = MALLOC(sizeof(struct FIVERB));
//    if (frirl->fiverb == NULL) {
//#ifdef DEBUG
//        printf("frirl_init: Error - malloc(frirl->fiverb)\n");
//#endif
//        return -1;
//    }
//    tmpfrb=frirl->fiverb;

//    frb->u = u;
//    frb->univlength = usize;
//    frb->numofunivs = frirl->numofantecedents;
    fri_float *ve = MALLOC(frirl->numofantecedents * usize * sizeof(fri_float));
    fri_float *rant = MALLOC(frirl->five_maxnumofrules * (frirl->numofantecedents) * sizeof(fri_float));
    fri_float *rconc = MALLOC(frirl->five_maxnumofrules * sizeof(fri_float));
    int numofrules;
    frirl->reduction_state = 0;

//    frirl_init_five(frirl);
    frirl_init_ve(frirl, ve, u, usize);
    frirl_init_rb(frirl, rant, rconc, &numofrules);
//    frb->numofrules = numofrules;
    int rulelength = frirl->numofantecedents + 1;


//    realnumofrules = numofrules; //printf("rules: %d\n", realnumofrules);
    //    frirl->frb=FIVEInit(u, frirl->frb->ve, 0, frirl->nantecedents, usize, frirl->frb->numofrules, frirl->frb->rulelength, frirl->frb->rb);
    frirl->fiverb = FIVEInit(u, ve, 0, frirl->numofantecedents, usize,
                             numofrules, frirl->five_maxnumofrules, rulelength, rant, rconc);
    struct FIVERB *frb = frirl->fiverb;

//    frb->numofrules = realnumofrules;//printf("rules: %d\n", realnumofrules);
//    frb->newrule = frb->rb + (frb->numofrules) * frb->rulelength;
//    frb->newrconc = frb->rconc + frb->numofrules;
//    frb->newrant = frb->rant + (frb->numofrules) * (frb->rulelength - 1);
//    frirl->rules_len = realnumofrules;
    //frirl->skip_rules=0;
//#ifdef BUILD_AVX2
//        frb->avx2_rbsize = realnumofrules / 4 + (realnumofrules % 4 > 0);
//#endif

    frirl->reward.ep_total_value = -1;
    frirl->reward.ep_total_steps = -1;
    frirl->fus_is_rule_inserted = 0;

    frirl->fiverb_vea = ve + usize * frirl->statedims_len;
    frirl->fiverb_ua  = u  + usize * frirl->statedims_len;
    // precalculate VE values for the possible actions - need to do this after frirl_init_five()
    for (j = 0; j < frirl->actiondim.values_len; j++) {
        frirl->possible_actions->vevalues[j] = frirl->fiverb_vea[get_vag_abs_min_i_fixres(frirl->fiverb_ua, frb->uksize, frirl->possible_actions->values[j], frb->udivs[frb->numofunivs-1])];
    }

    // check for some obvious errors in given parameters
    // u[0]-u[max] == 0 ?
    for (k = 0; k < frb->numofunivs; k++) {
        knu = k * frb->univlength;

        if ((frb->u[knu + frb->univlength - 1] - frb->u[knu]) == 0.0) {
            printf("frirl_init: FATAL ERROR: U dimension %d max-min=0.0! Did you cast the constant to double explicitly? e.g. 1001.0 instead of 1001.\n", k);
            exit(1);
        }
    }

    // check state dims
    for (k = 0; k < (frb->numofunivs - 1); k++) {
        knu = k * frb->univlength;
        //do the state intervals fit their "u"-s?
        //    printf("k.: %d knu: %d frirl->possiblestates[k].statedims: %d min: %f max: %f\n",k, knu, frirl->possiblestates[k].statedims, frirl->possiblestates[k].statedim[0], frirl->possiblestates[k].statedim[frirl->possiblestates[k].statedims-1]);
        if ( (frirl->possible_states[k].values[0] < frb->u[knu]) ||
             ((frirl->possible_states[k].values[frirl->possible_states[k].values_len - 1] > frb->u[(k + 1) * frb->univlength - 1])) ) {
            printf("frirl_init: FATAL ERROR: The points are out of range for statedim%d. Should be between %lf and %lf, but found %lf - %lf! Check u generation and the statedim boundaries!\n", k, frb->u[knu], frb->u[(k + 1) * frb->univlength - 1], frirl->possible_states[k].values[0], frirl->possible_states[k].values[frirl->possible_states[k].values_len - 1]);
            exit(1);
        }
        //default state inside interval?
        //    printf("( frirl->possiblestates[k].statedim[0] < frirl->states_default[k] ): %d\n",( frirl->possiblestates[k].statedim[0] > frirl->states_default[k] ));
        //    printf("( frirl->possiblestates[k].statedim[frirl->possiblestates[k].statedims-1] > frirl->states_default[k] ): %d\n",( frirl->possiblestates[k].statedim[frirl->possiblestates[k].statedims-1] < frirl->states_default[k] ));
        if ( (frirl->possible_states[k].values[0] > frirl->statedims[k].values_def) ||
             (frirl->possible_states[k].values[frirl->possible_states[k].values_len - 1] < frirl->statedims[k].values_def) ) {
            printf("WARNING: default state value for statedim%d is out of range! Should be true: %lf < %lf < %lf\n", k, frirl->possible_states[k].values[0], frirl->statedims[k].values_def, frirl->possible_states[k].values[frirl->possible_states[k].values_len - 1]);
        }
        //check divider vs. actual steps
        for (st = 1; st < frirl->possible_states[k].values_len; st++) {
            //printf("%d[%d]=%lf ",k,st,frirl->possiblestates[k].statedim[st]);
            //printf("%lf ",frirl->possiblestates[k].statedim[st]-frirl->possiblestates[k].statedim[st-1]);
            if ( !( (frirl->possible_states[k].values[st] - frirl->possible_states[k].values[st - 1]) > (frirl->possible_states[k].epsilon - 0.00001) &&
                    (frirl->possible_states[k].values[st] - frirl->possible_states[k].values[st - 1]) < (frirl->possible_states[k].epsilon + 0.00001) ) ) {
                printf("frirl_init: WARNING: state dim step or divider is possibly wrong (or innaccurate): statedim%d[%d]=%lf statedim%d[%d]=%lf diff is: %lf and divider is: %lf\n", k, st - 1, frirl->possible_states[k].values[st - 1], k, st, frirl->possible_states[k].values[st], (frirl->possible_states[k].values[st] - frirl->possible_states[k].values[st - 1]), frirl->possible_states[k].epsilon);
            }
        }
        //    printf(" %f %f %f\n",frirl->possiblestates[k].epsilon,frirl->possiblestates[k].epsilon-0.00001,frirl->possiblestates[k].epsilon+0.00001);
    }

    // also check it for the action dim
    k = frb->numofunivs - 1;
    knu = k * frb->univlength;
    //    printf("action: %d knu: %d frirl->possibleaction.statedims: %d min: %f max: %f\n",k, knu, frirl->possibleaction->statedims, frirl->possibleaction->statedim[0], frirl->possibleaction->statedim[frirl->possibleaction->statedims-1]);
    if ( (frirl->possible_actions->values[0] < frb->u[knu]) ||
         ((frirl->possible_actions->values[frirl->possible_actions->values_len - 1] > frb->u[(k + 1) * frb->univlength - 1])) ) {
        printf("frirl_init: FATAL ERROR: The points are out of range for the action dimension. Should be between %lf and %lf, but found %lf - %lf! Check u generation and the statedim boundaries for the action dimension!\n", frb->u[knu], frb->u[(k + 1) * frb->univlength - 1], frirl->possible_actions->values[0], frirl->possible_actions->values[frirl->possible_actions->values_len - 1]);
        exit(1);
    }

    // check divider vs. actual steps for action dim
    for (st = 1; st < frirl->possible_actions->values_len; st++) {
        //    printf("actions[%d]=%lf %lf\n",st,frirl->possibleaction->statedim[st],frirl->possibleaction->statedim[st]-frirl->possibleaction->statedim[st-1]);
        if ( !( (frirl->possible_actions->values[st] - frirl->possible_actions->values[st - 1]) > (frirl->possible_actions->epsilon - 0.00001) &&
                (frirl->possible_actions->values[st] - frirl->possible_actions->values[st - 1]) < (frirl->possible_actions->epsilon + 0.00001) ) ) {
            printf("frirl_init: WARNING: action dim step or divider is possibly wrong (or innaccurate): actions[%d]=%lf actions[%d]=%lf diff is: %lf and divider is: %lf\n", st - 1, frirl->possible_actions->values[st - 1], st, frirl->possible_actions->values[st], (frirl->possible_actions->values[st] - frirl->possible_actions->values[st - 1]), frirl->possible_actions->epsilon);
        }
    }


    frirl->fgba_vagdist_states = MALLOC(sizeof(fri_float) * (frirl->five_maxnumofrules + 4) * (frb->numofunivs - 1));
    if (frirl->fgba_vagdist_states == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(frirl->fgba_vagdist_states)\n");
#endif
        return -1;
    }

    frirl->fgba_ruledist = MALLOC(sizeof(fri_float) * frirl->five_maxnumofrules);
    if (frirl->fgba_ruledist == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(frirl->fgba_ruledist)\n");
#endif
        return -1;
    }
    memset(frirl->fgba_ruledist, 0, sizeof(fri_float) * frirl->five_maxnumofrules);

    frirl->fgba_actconc = MALLOC(sizeof(fri_float) * frirl->five_maxnumofrules);
    if (frirl->fgba_actconc == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(frirl->fgba_actconc)\n");
#endif
        return -1;
    }

    frirl->fgba_dists = (fri_float *)MALLOC((frirl->five_maxnumofrules + 4) * sizeof(fri_float) * FIVE_MAX_NUM_OF_UNIVERSES);
    if (frirl->fgba_dists == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(frirl->fgba_dists)\n");
#endif
        return -1;
    }
    memset(frirl->fgba_dists, 0, (frirl->five_maxnumofrules + 4) * sizeof(fri_float) * FIVE_MAX_NUM_OF_UNIVERSES);

    frirl->fgba_statedistsum = MALLOC(sizeof(fri_float) * frirl->five_maxnumofrules);
    if (frirl->fgba_statedistsum == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(frirl->fgba_statedistsum)\n");
#endif
        return -1;
    }

    frirl->fus_values = MALLOC(sizeof(fri_float) * frb->rulelength);
    if (frirl->fus_values == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(frirl->fus_values)\n");
#endif
        return -1;
    }

    frirl->fus_proposed_values = MALLOC(sizeof(fri_float) * frb->rulelength);
    if (frirl->fus_proposed_values == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(frirl->fus_proposed_values)\n");
#endif
        return -1;
    }

    frirl->fus_check_states = MALLOC(sizeof(fri_float) * frb->numofantecedents);
    if (frirl->fus_check_states == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc(frirl->fus_check_states)\n");
#endif
        return -1;
    }

    frirl->fep_ant = MALLOC(sizeof(fri_float) * frirl->numofantecedents);
    if (frirl->fep_ant == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc fep_rant\n");
#endif
        return -1;
    }
    frirl->fep_cur_ant = MALLOC(sizeof(fri_float) * frirl->numofantecedents);
    if (frirl->fep_cur_ant == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc fep_cur_rant\n");
#endif
        return -1;
    }
    frirl->fep_q_ant = MALLOC(sizeof(fri_float) * frirl->numofantecedents);
    if (frirl->fep_q_ant == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc fep_q_rant\n");
#endif
        return -1;
    }
    frirl->fep_cur_q_ant = MALLOC(sizeof(fri_float) * frirl->numofantecedents);
    if (frirl->fep_cur_q_ant == NULL) {
#ifdef DEBUG
        printf("frirl_init: Error - malloc fep_cur_q_rant\n");
#endif
        return -1;
    }

#ifdef BUILD_MPI
    //tons of todo

    MPI_Init(NULL, NULL);

    // get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &frirl->agent_world_size);

    // rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &frirl->agent_id);

    // whoami
//    int mpi_whoami_len; //it's a string, who cares about the length?
//    MPI_Get_processor_name(frirl->mpi_name, &mpi_whoami_len);
#endif

    frirl->is_running = 1;
    frirl->episode_num = 1;

    // random seed for random action selection
    gettimeofday(&t1, NULL);
    srand(t1.tv_usec * t1.tv_sec);

//    free(tmpfrb);

    return 0;
}
