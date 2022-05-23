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

#ifndef _FRIRL_H
#define _FRIRL_H

#include "config.h"
#include "frirl_types.h"
#include "frirl_utils.h"

#include "FIVE.h"

//#define TERM_GRAY  "\033[0;30m"
#define TERM_RED   "\033[0;31m"
#define TERM_GREEN "\033[0;32m"
#define TERM_NC    "\033[0m"

/*
//disable
#define TERM_RED   ""
#define TERM_GREEN ""
#define TERM_NC    ""
*/

///RB size is fixed - more faster and simpler than reallocating every time a new rule is added
//#define FRIRLRULES 16384

//API rc5 (obsoleted, kept for backward compatibility)
//todo restore rc5 API
//todo resolve naming inconsistency
//void frirl_update_sarsa_rc5(struct frirl *frirl, double *s, double a, double r, double *sp, double ap);

//API by daniel (obsoleted, to be removed)
int frirl_init(struct frirl_desc *frirl);
void frirl_deinit(struct frirl_desc *frirl);

int frirl_init_ve(struct frirl_desc *frirl, fri_float *ve, fri_float *u, int univlength);
int frirl_init_rb(struct frirl_desc *frirl, fri_float *rant, fri_float *rconc, int *numofrules);

// TODO to be removed
//void frirl_update_sarsa_old(struct frirl_desc *frirl, fri_float *states, fri_float action,
//                            fri_float rule, fri_float *proposed_states, fri_float proposed_action);

void frirl_episode(struct frirl_desc *frirl);

unsigned int frirl_e_greedy_selection(struct frirl_desc *frirl, fri_float *states);

unsigned int frirl_get_best_action(struct frirl_desc *frirl, fri_float *states);

fri_float frirl_check_possible_states(struct frirl_desc *frirl, fri_float observation, struct frirl_values_desc *possible_states);

void frirl_sequential_run(struct frirl_desc *frirl);
void frirl_omp_run(struct frirl_desc *frirl);
void frirl_mpi_run(struct frirl_desc *frirl);

void frirl_update_sarsa(struct frirl_desc *frirl, fri_float *q_ant, fri_float reward, fri_float *cur_q_ant);

#endif //FRIRL_H
