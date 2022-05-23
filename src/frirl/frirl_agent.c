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

#include "frirl_types.h"
#include "frirl.h"
#include <string.h>
#include <stdio.h>
#include <math.h>


static void check_possible_states(struct frirl_desc *frirl, fri_float *ant, fri_float *newant) {

    unsigned int i;

    for (i = 0; i < frirl->statedims_len; i++) {
//        [newrulestates(currps)]=frirl_check_possible_states(s(currps),possiblestates{currps},possiblestates_epsilons{currps}); // MATLAB code
        DEBUG_MSG("frirl_update_SARSA: currps: %d frirl->possiblestates[currps]->statedims: %d\n", i, frirl->possible_states[i].values_len);

        newant[i] = frirl_check_possible_states(frirl, ant[i], &frirl->possible_states[i]);

        DEBUG_MSG("frirl_update_SARSA: newrulestates[%d]: %f\n", i, new_rule_states[i]);

        if (newant[i] == INFINITY) {
            printf("FATAL: frirl_update_SARSA: CHKPOSSSTATEBUG: currps: %d\n", i);
            exit(2);
        }

    }

    newant[i] = frirl_check_possible_states(frirl, ant[i], frirl->possible_actions);
    if (newant[i] == INFINITY) {
        printf("FATAL: frirl_update_SARSA: POSSACTBUG\n");
        exit(2);
    }

}

static void update_rules(int numofrules, fri_float *weights, fri_float *rconc, fri_float q, fri_float delta) {

    for (int w = 0; w < numofrules; w++) {
        if (weights[w] > delta) {
            rconc[w] = q * weights[w];
        }
    }

}


/// openmp thread sync
static void merge_rb(struct frirl_desc *rcvr_frirl, int fromid, fri_float *newrant, fri_float *newrconc, int numofrules) {
#ifdef BUILD_OPENMP

//    if (sndr_frirl->is_running == 0 || rcvr_frirl->is_running == 0)
//        goto exit;
//        return;

    if (rcvr_frirl->verbose) {
        printf("Synchronizing the rulebase from #%d to #%d\n", fromid, rcvr_frirl->agent_id);
    }

//    struct FIVERB *sndr_frb = sndr_frirl->fiverb;
    struct FIVERB *rcvr_frb = rcvr_frirl->fiverb;

    for (int r = 0; r < numofrules; r++) {

        fri_float *sndr_rant = newrant + r * rcvr_frirl->numofantecedents;
        fri_float sndr_rconc = newrconc[r];

        fri_float *rcvr_weights = rcvr_frb->weights;
        int rcvr_rule_i = FIVE_vag_concl_weight(rcvr_frb, sndr_rant, rcvr_weights);

        fri_float rcvr_rconc;
        FIVE_vag_concl(rcvr_frb, sndr_rant, &rcvr_rconc);
//        fri_float qdiff = sndr_frirl->alpha * (rcvr_rconc - sndr_frirl->gamma * sndr_rconc);
//        fri_float qdiff = (rcvr_rconc - sndr_rconc);
//        fri_float conclusion = 0.5 * rcvr_rconc + 0.5 * sndr_rconc;
//        fri_float qdiff = 0.1 * (-rcvr_rconc + 0.75 * rcvr_rconc + 0.25 * sndr_rconc);
        fri_float qdiff = -rcvr_rconc + sndr_rconc;
//        printf("%f %f %f\n", qdiff, rcvr_rconc, sndr_rconc);
//        fri_float qdiff = 0;
        fri_float delta = rcvr_frirl->rule_weight_considered_significant_for_update;
        fri_float qdiff_pos_boundary = rcvr_frirl->qdiff_pos_boundary;
        fri_float qdiff_neg_boundary = rcvr_frirl->qdiff_neg_boundary;

        if ((qdiff > qdiff_pos_boundary) || (qdiff < qdiff_neg_boundary)) {
            fri_float *rant = sndr_rant;
            fri_float rconc = rcvr_rconc;
            if (CHECK_STATES) {
                rant = rcvr_frirl->fus_check_states;
                check_possible_states(rcvr_frirl, sndr_rant, rant);
                rcvr_rule_i = FIVE_vag_concl(rcvr_frb, rant, &rconc);
            }
            if (rcvr_rule_i == ~0) {
                FIVE_add_rule(rcvr_frb, rant, 0.5 * rconc + 0.5 * sndr_rconc);
//                goto exit;
                continue;
            }
            rcvr_frb->rconc[rcvr_rule_i] = 0.9 * rconc + 0.1 * sndr_rconc;
//            goto exit;
            continue;
        }

        update_rules(rcvr_frb->numofrules, rcvr_weights, rcvr_frb->rconc, 0.9 * rcvr_rconc + 0.1 * sndr_rconc, delta);
    }

//exit:
    return;
#endif
}


///generate statedims default values with the highest distance between them
static void gen_def_states(struct frirl_desc *frirl) {

    int id = frirl->agent_id;
    int worldsize = frirl->agent_world_size;

    if (id == 0 || worldsize < 2)
        return;

    fri_float *rant = frirl->fiverb->rant;
    int numofrules = frirl->fiverb->numofrules;
    int numofantecedents = frirl->fiverb->numofantecedents;
    int gap = numofrules / (worldsize - 2);
    int i;

    for (i = 0; i < frirl->statedims_len; i++) {
        frirl->statedims[i].values_def = rant[(id - 1) * gap * numofantecedents + i];
    }

}

///initalize a new frirl instance (for one thread)
///\param frirl initialized frirl struct
///\param new frirl struct
///\return void
static void omp_init(struct frirl_desc *frirl, struct frirl_desc *ompfrirl, int id, int world_size) {

    memcpy(ompfrirl, frirl, sizeof(struct frirl_desc));

#ifdef BUILD_OPENMP
//    //is there anyhing more needed?
//    frirl_init(ompfrirl);

//    ompfrirl->agent_id = omp_get_thread_num();
    ompfrirl->agent_id = id;
//    ompfrirl->agent_world_size = omp_get_num_threads();
    ompfrirl->agent_world_size = world_size;
    ompfrirl->is_running = 1;

//    if (ompfrirl->agent_id != 0)
//        ompfrirl->no_random = 0;

    int rulelength = frirl->numofantecedents + 1;

    //copy frirl stuffs
    ompfrirl->statedims = MALLOC(sizeof(struct frirl_dimension_desc) * frirl->statedims_len);
    for (unsigned int i = 0; i < frirl->statedims_len; i++) {
        memcpy(ompfrirl->statedims + i, frirl->statedims + i, sizeof(struct frirl_dimension_desc));
//        ompfrirl->statedims[i].values = MALLOC(sizeof(double) * frirl->statedims[i].values_len);
//        memcpy(ompfrirl->statedims[i].values, frirl->statedims[i].values, sizeof(double) * frirl->statedims[i].values_len);
        /*if (ompfrirl->agent_id > 0) {
            int rand_i = rand() % frirl->statedims[i].values_len;
                        printf("omp: thread: %d statedim: %d def_value: %d\n", ompfrirl->agent_id, i, rand_i);
            ompfrirl->statedims[i].values_def = frirl->statedims[i].values[rand_i];
        }*/
    }
    gen_def_states(ompfrirl);
    for (int i = 0; i < ompfrirl->statedims_len; i++) {
        printf("#%d sdim%d def:%f\n", ompfrirl->agent_id, i,
            ompfrirl->statedims[i].values_def);
    }
//    ompfrirl->actiondim = MALLOC(sizeof(struct frirl_dimension_desc));
//    memcpy(ompfrirl->actiondim, frirl->actiondim, sizeof(struct frirl_dimension_desc));
//    ompfrirl->actiondim.values_def = rand() % frirl->actiondim.values_len;

    ompfrirl->possible_states = MALLOC(sizeof(struct frirl_dimension_desc)*frirl->statedims_len);
    for (unsigned int i = 0; i < frirl->statedims_len; i++) {
        ompfrirl->possible_states[i].values_len = frirl->statedims[i].values_len;
        ompfrirl->possible_states[i].values = MALLOC(sizeof(double)*frirl->statedims[i].values_len);

        for (unsigned int j = 0; j < frirl->statedims[i].values_len; j++)
            ompfrirl->possible_states[i].values[j] = frirl->statedims[i].values[j];

        ompfrirl->possible_states[i].epsilon = frirl->statedims[i].values_div;
        ompfrirl->possible_states[i].vevalues = frirl->possible_states[i].vevalues;
    }

    ompfrirl->possible_actions = MALLOC(sizeof(struct frirl_values_desc));
    ompfrirl->possible_actions->values_len = frirl->actiondim.values_len;

    ompfrirl->possible_actions->values = MALLOC(sizeof(fri_float)*frirl->actiondim.values_len);
    for (unsigned int i = 0; i < frirl->actiondim.values_len; i++)
        ompfrirl->possible_actions->values[i] = frirl->actiondim.values[i];

    ompfrirl->possible_actions->epsilon  = frirl->actiondim.values_div;
    ompfrirl->possible_actions->vevalues = frirl->possible_actions->vevalues;

    ompfrirl->fgba_vagdist_states = MALLOC(sizeof(fri_float) * (frirl->five_maxnumofrules + 4) * frirl->statedims_len);
    ompfrirl->fgba_ruledist       = MALLOC(sizeof(fri_float) * frirl->five_maxnumofrules);
    memset(ompfrirl->fgba_ruledist, 0, sizeof(fri_float) * frirl->five_maxnumofrules);
    ompfrirl->fgba_actconc = MALLOC(sizeof(fri_float) * frirl->five_maxnumofrules);
    ompfrirl->fgba_dists   = (fri_float *)MALLOC((frirl->five_maxnumofrules + 4) * sizeof(fri_float) * FIVE_MAX_NUM_OF_UNIVERSES);
    memset(frirl->fgba_dists, 0, (frirl->five_maxnumofrules + 4) * sizeof(fri_float) * FIVE_MAX_NUM_OF_UNIVERSES);
    ompfrirl->fgba_statedistsum = MALLOC(sizeof(fri_float) * frirl->five_maxnumofrules);

    ompfrirl->fus_values          = MALLOC(sizeof(fri_float) * rulelength);
    ompfrirl->fus_proposed_values = MALLOC(sizeof(fri_float) * rulelength);
    ompfrirl->fus_check_states    = MALLOC(sizeof(fri_float) * frirl->numofantecedents);

    ompfrirl->fep_ant       = MALLOC(sizeof(fri_float) * frirl->numofantecedents);
    ompfrirl->fep_q_ant     = MALLOC(sizeof(fri_float) * frirl->numofantecedents);
    ompfrirl->fep_cur_ant   = MALLOC(sizeof(fri_float) * frirl->numofantecedents);
    ompfrirl->fep_cur_q_ant = MALLOC(sizeof(fri_float) * frirl->numofantecedents);

    // copy from FIVE
    fri_float *u = frirl->fiverb->u;
    int usize = frirl->statedims[0].universe_len;

    fri_float *ve = frirl->fiverb->ve; //MALLOC(frirl->numofantecedents * usize * sizeof(fri_float));
    fri_float *rant = MALLOC(frirl->five_maxnumofrules * (frirl->numofantecedents) * sizeof(fri_float));
    fri_float *rconc = MALLOC(frirl->five_maxnumofrules * sizeof(fri_float));
    int numofrules = frirl->fiverb->numofrules;
    ompfrirl->reduction_state = 0;

//    frirl_init_ve(ompfrirl, ve, u, usize);
    frirl_init_rb(ompfrirl, rant, rconc, &numofrules);

    ompfrirl->fiverb = FIVEInit(u, ve, 0, frirl->numofantecedents, usize,
                                numofrules, frirl->five_maxnumofrules, rulelength, rant, rconc);

    ompfrirl->reward.ep_total_value = -1;
    ompfrirl->reward.ep_total_steps = -1;
    ompfrirl->fus_is_rule_inserted = 0;

//    ompfrirl->fiverb_vea = ve + usize * frirl->statedims_len;
//    ompfrirl->fiverb_ua  = u  + usize * frirl->statedims_len;

    //frirl_show_rb(ompfrirl);
#endif //BUILD_OPENMP

}


///destroy copied frirl struct
static void omp_deinit(struct frirl_desc *ompfrirl) {
#ifdef BUILD_OPENMP

    free(ompfrirl->statedims);

    for (unsigned int i = 0; i < ompfrirl->statedims_len; i++) {
        free(ompfrirl->possible_states[i].values);
    }
    free(ompfrirl->possible_states);
    free(ompfrirl->possible_actions->values);
    free(ompfrirl->possible_actions);

    free(ompfrirl->fgba_vagdist_states);
    free(ompfrirl->fgba_ruledist);
    free(ompfrirl->fgba_actconc);
    free(ompfrirl->fgba_dists);
    free(ompfrirl->fgba_statedistsum);

    free(ompfrirl->fus_values);
    free(ompfrirl->fus_proposed_values);
    free(ompfrirl->fus_check_states);

    free(ompfrirl->fep_ant);
    free(ompfrirl->fep_q_ant);
    free(ompfrirl->fep_cur_ant);
    free(ompfrirl->fep_cur_q_ant);

    free(ompfrirl->fiverb->rant);
    free(ompfrirl->fiverb->rconc);
    five_deinit(ompfrirl->fiverb);

//    free(ompfrirl);

#endif //BUILD_OPENMP

}

///FRIRL main loop (with openmp)
///\param frirl initialized frirl stuct
///\return void
void frirl_omp_run(struct frirl_desc *frirl) {

#ifndef BUILD_OPENMP
    printf("Built without openmp support!\n");
#else
//todo
    int agent_world_size /*= 4*/;
#pragma omp parallel
    {
        agent_world_size = omp_get_num_threads();
    }

    frirl->agent_id = -1;
    frirl->agent_world_size = agent_world_size;

    struct frirl_desc agent[agent_world_size];
#pragma omp parallel for
    for (int id = 0; id < agent_world_size; id++) {
        printf("initializing thread %d\n", id);
        omp_init(frirl, &agent[id], id, agent_world_size);
    }

    for(;;) {

#pragma omp parallel for
        for (int id = 0; id < agent_world_size; id++) {
            /*if (agent[id].is_running == 0) {
                printf("initializing thread %d\n", id);
                omp_init(frirl, &agent[id], id, agent_world_size);
            }*/
            frirl_sequential_run(&agent[id]);
        }
//#pragma omp barrier

        if (agent[0].is_running == 0) {
//            printf("vege\n");
            break;
        }

        if (agent[0].epended == 0) {
            struct FIVERB *tmpfrb = agent[0].fiverb;
    //        fri_float master_reward = agent[0].reward.ep_total_value;
            for (int id = 1; id < agent_world_size; id++) {
    //            fri_float slave_reward = agent[id].reward.ep_total_value;
    //            if (/*agent[id].is_running != 0 &&*/ slave_reward > master_reward)

                    merge_rb(&agent[id], 0, tmpfrb->rant, tmpfrb->rconc, tmpfrb->numofrules);
            }
        } else {
            printf("#0 pended, skipping sync to master\n");
        }
//#pragma omp parallel for

        for (int id = 1; id < agent_world_size; id++) {
            if (agent[id].epended == 0) {
                struct FIVERB *tmpfrb = agent[id].fiverb;
    //            fri_float slave_reward = agent[id].reward.ep_total_value;
//                if (agent[id].is_running != 0 /*&& slave_reward < master_reward*/);;
                merge_rb(&agent[0], id, tmpfrb->rant, tmpfrb->rconc, tmpfrb->numofrules);
            } else {
                printf("#%d pended, skipping sync\n", id);
            }
        }

////#pragma omp barrier

        /*for (int id = 1; id < agent_world_size; id++) {
            if (agent[id].is_running == 0) {
                printf("deinitializing finished thread %d\n", id);
                omp_deinit(&agent[id]);
            }
        }*/
//        if (agent[0].is_running == 0) printf("kaputt\n");
    }

    printf("moving rb from thread #0 to main process\n");
    frirl->fiverb->newrant = frirl->fiverb->rant;
    frirl->fiverb->newrconc = frirl->fiverb->rconc;
    frirl->fiverb->numofrules = 0;
    for (int i = 0; i < agent[0].fiverb->numofrules; i++) {
        fri_float *rant = agent[0].fiverb->rant + i * frirl->numofantecedents;
        fri_float rconc = agent[0].fiverb->rconc[i];
        FIVE_add_rule(frirl->fiverb, rant, rconc);
    }

    for (int id = 0; id < agent_world_size; id++) {
        printf("deinitializing thread %d\n", id);
        omp_deinit(&agent[id]);
    }
#endif

}

void frirl_mpi_run(struct frirl_desc *frirl) {

#ifndef BUILD_MPI
    printf("Built without MPI support!\n");
#else
//todo

    /*if (frirl->agent_rnd_init != 0 && frirl->agent_id != 0) {
        for (int i = 0; i < frirl->statedims_len; i++) {
            frirl->statedims[i].values_def =
                frirl->statedims[i].values[rand() % frirl->statedims[i].values_len];
        }
    }*/
    gen_def_states(frirl);
    for (int i = 0; i < frirl->statedims_len; i++) {
        printf("#%d sdim%d def:%f\n", frirl->agent_id, i,
            frirl->statedims[i].values_def);
    }

    int agent_world_size = frirl->agent_world_size;
    int agent_id = frirl->agent_id;
    struct FIVERB *frb = frirl->fiverb;
    fri_float *rant = MALLOC(frb->numofantecedents * frb->maxnumofrules * sizeof(fri_float));
    fri_float *rconc = MALLOC(frb->maxnumofrules * sizeof(fri_float));
    int rules_num = 0;


    for (int id = 0; id < frirl->agent_world_size; id++) {
        for (int s = 0; s < frirl->statedims_len; s++) {
            printf("#%d statedim: %d def: %f\n", id, s, frirl->statedims[s].values_def);
        }
    }

#pragma omp parallel
#pragma omp single
    for(;;) {

        frirl_sequential_run(frirl);

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&frirl->is_running, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (frirl->is_running == 0) {
            printf("End\n");
            break;
        }

        if (agent_id == 0) {
            //var a tobbiekre, sync onnan
            for (int id = 1; id < agent_world_size; id++){
                MPI_Recv(&rules_num, 1, MPI_INT, id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(rant, rules_num * frb->numofantecedents, MPI_DOUBLE, id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(rconc, rules_num, MPI_DOUBLE, id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                merge_rb(frirl, id, rant, rconc, rules_num);
            }
        } else {
            // sync - send to main controller
            MPI_Send(&frb->numofrules, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(frb->rant, frb->numofantecedents * frb->numofrules, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(frb->rconc, frb->numofrules, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (agent_id == 0) {
            for (int id = 1; id < agent_world_size; id++) {
                // sync - send to other agents
                MPI_Send(&frb->numofrules, 1, MPI_INT, id, 0, MPI_COMM_WORLD);
                MPI_Send(frb->rant, frb->numofantecedents * frb->numofrules, MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
                MPI_Send(frb->rconc, frb->numofrules, MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
            }
        } else {
            MPI_Recv(&rules_num, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(rant, rules_num * frb->numofantecedents, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(rconc, rules_num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            merge_rb(frirl, 0, rant, rconc, rules_num);
        }
    }
#endif

}
