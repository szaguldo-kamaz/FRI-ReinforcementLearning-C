/*
** Fuzzy Rule Interpolation-based Reinforcement Learning (ANSI C / AVX version)
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** Author: David Vincze <david.vincze@webcode.hu>
**
** Various contributions by Daniel Palko <palko.daniel@uni-miskolc.hu>
**
** Imitation related contributions by Alex Martossy
**
*/


//#define DEBUG

#include "frirl.h"
#include "frirl_imitation.h"
#ifdef DEBUG
#include <stdio.h>
#endif


///evaluate an episode
///\param frirl initialized FRIRL struct
///\param max_steps maximum steps
///\return void
void frirl_episode(struct frirl_desc *frirl) {

    fri_float *q_ant     = frirl->fep_q_ant;  // quantized antecedents (ant. values before doaction())
                                              // quantized means: fitted to the nearest grid value
                                              // antecedents means: <state0,state1,...,stateN,action>
    fri_float *cur_q_ant = frirl->fep_cur_q_ant;  // current quantized antecedents (ant. values after doaction())

    fri_float *states     = frirl->fep_ant;
    fri_float *cur_states = frirl->fep_cur_ant;
    fri_float *action     = q_ant + frirl->statedims_len;

    fri_float *q_states     = q_ant;
    fri_float *cur_q_states = cur_q_ant;
    fri_float *prop_action  = cur_q_ant + frirl->statedims_len;  // proposed action pointer

    unsigned int action_i;


    for (int i = 0; i < frirl->statedims_len; i++) {
        q_states[i] = states[i] = frirl->statedims[i].values_def;
    }

    frirl->reward.ep_total_value = 0;
    frirl->reward.ep_total_steps = 0;

#ifdef BUILD_VISUALIZATION
    frirl->draw_func(frirl, cur_states, frirl->actiondim.values_def, frirl->reward.ep_total_steps);
    // usleep(100);
#endif

    // imitation
    if(frirl->original_learning == 0) {

        while(1) {
            getActionFromTerminal(frirl);

            if (frirl->valid_simulation == 1) {
                if (frirl->keyaction == 32) {
                    action_i = frirl_e_greedy_selection(frirl, states);
                } else {
                    action_i = frirl->keyaction;
                }
                frirl->valid_simulation = 0;
                break;
            }

        }

    } else {
        // select an action using the epsilon greedy selection strategy
        action_i = frirl_e_greedy_selection(frirl, states);
    }

    // convert the index of the action into an action value
    *action = frirl->actiondim.values[action_i];

    DEBUG_MSG("frirl_episode: %d, initial a: %d\n", frirl->episode_num, action_i);

    for (int step_num = 1; step_num <= frirl->max_steps; step_num++) {
#ifdef DEBUG
        printf("STEPNO: %d\n", step_num);
        printf("FRIRL:_episode: acno: %d action: %f\n", action_i, action);
        printf("FRIRL_episode: state before doaction: ");
        for (unsigned int i = 0; i < frirl->statedims_len; i++) {
            printf("%f ", ant[i]);  // was x[i]
        }
        printf("\n");
#endif
        // do action
        frirl->do_action_func(frirl, *action, states, frirl->statedims_len, cur_states);
#ifdef DEBUG
        printf("FRIRL_episode: state after doaction: ");
        for (unsigned int i = 0; i < frirl->statedims_len; i++) {
            printf("%f ", cur_ant[i]);
        }
        printf("\n");
#endif
        // observe the reward at state 'cur_state' and set the final state flag
        frirl->get_reward_func(frirl, cur_states, frirl->statedims_len, &frirl->reward);
        frirl->reward.ep_total_value += frirl->reward.value;
        DEBUG_MSG("FRIRL_episode: step_num: %d reward: %.15f\n", step_num, frirl->reward.value);
        //DEBUG_MSG("FRIRL_episode: step_num: %d reward: %.13a\n",step_num,frirl->reward.value);

        // quantize observations
        frirl->quant_obs_func(frirl, cur_states, frirl->statedims_len, cur_q_states);
#ifdef DEBUG
        printf("FRIRL_episode: state after quant (sp): ");
        for (unsigned int i = 0; i < frirl->statedims_len; i++) {
            printf("%f ", cur_q_ant[i]);
        }
        printf("\n");

#endif

#ifdef BUILD_VISUALIZATION
        frirl->draw_func(frirl, cur_states, frirl->actiondim.values_def, frirl->reward.ep_total_steps + 1);
        // usleep(100);
#endif

        // imitation
        unsigned int proposed_action_i;

        if (frirl->original_learning == 0) {

            while(1) {
                getActionFromTerminal(frirl);

                if (frirl->valid_simulation == 1) {
                    if (frirl->keyaction == 32) {
                        proposed_action_i = frirl_e_greedy_selection(frirl, cur_q_states);
                    } else {
                        proposed_action_i = frirl->keyaction;
                        frirl->valid_simulation = 0;
                        break;
                    }
                }
            }

        } else {
            // propose an action for the current situation
            proposed_action_i = frirl_e_greedy_selection(frirl, cur_q_states);
        }

        *prop_action = frirl->actiondim.values[proposed_action_i];
        DEBUG_MSG("frirl_episode: greedy action: index: %d, %f\n", proposed_action_i, *prop_action);

        // update the rule-base conclusions
        if (frirl->reduction_state == 0) {
            //printf("step_num: %d\n", step_num); // debug
//            frirl_update_sarsa(frirl, qobs, action, frirl->reward.value, cur_qobs, proposed_action);
//            frirl_update_sarsa(frirl, qobs, *action, frirl->reward.value, cur_qobs, *pr_action);
            frirl_update_sarsa(frirl, q_ant, frirl->reward.value, cur_q_ant);
        }

        // update the current variables
        for (int i = 0; i < frirl->statedims_len; i++) {
            states[i] = cur_states[i];
        }
        for (int i = 0; i < frirl->numofantecedents; i++) {
            q_ant[i] = cur_q_ant[i];
        }
//        for (int i = 0; i < frirl->statedims_len; i++) {
//            q_states[i] = cur_q_states[i];
//        }
//        *action = *prop_action;

        frirl->reward.ep_total_steps++;
//        frirl_show_rb_hex(frirl);

//#ifdef BUILD_VISUALIZATION
//        frirl->draw_func(frirl, cur_states, *action, frirl->reward.ep_total_steps);
//        // usleep(100);
//#endif

        // if goal is reached then break the episode
        if (frirl->reward.success == 1) {
            break;
        }

    }

#ifdef BUILD_VISUALIZATION
    frirl->draw_func(frirl, cur_states, *action, frirl->reward.ep_total_steps);
    // usleep(100);
#endif

}
