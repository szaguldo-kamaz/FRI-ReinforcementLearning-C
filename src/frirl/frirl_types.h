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


#ifndef FRIRL_TYPES_H
#define FRIRL_TYPES_H

#include "config.h"


enum frirl_runmode {
    FRIRL_SEQ,
    FRIRL_OMP,
    FRIRL_MPI,
    FRIRL_TEST
};

enum frirl_reduction_strategy {
    FRIRL_REDUCTION_STRATEGY_NOREDUCE,
    FRIRL_REDUCTION_STRATEGY_DEFAULT
};

///universe and state / action (dimension) desc
struct frirl_dimension_desc {
    int values_len;         ///<values array length
    fri_float *values;      ///<values array
    fri_float values_div;   ///<step between values
    fri_float values_steep; ///<steepness of the values
    fri_float values_def;   ///<default value
    int universe_len;       ///<universe array length
    fri_float *universe;    ///<universe array
    fri_float universe_div; ///<step between universe values
};

///possible states / actions desc
struct frirl_values_desc {
    int values_len;      ///<values array length (number of actions)
    fri_float *values;   ///<stores the states/action values
    fri_float *vevalues; ///<to store the precalculated VE values for the possible actions
    fri_float epsilon;   ///<epsilon value
};

///reward desc
struct frirl_reward_desc {
    fri_float value;          ///<reward value
    fri_float ep_total_value; ///<episode total reward
    int ep_total_steps;       ///<episode total steps
    int success;              ///<success
};

///common frirl struct
struct frirl_desc {

    int argc;
    char **argv;
    int runmode;
    char *rbfile;
    int visualization;
    int gui_width;
    int gui_height;
    int verbose;
    int agent_rnd_init; ///<randomize the agent's state default value

    struct frirl_dimension_desc actiondim;  ///<action dimension
    int statedims_len;                     ///<state dimensions array length
    struct frirl_dimension_desc *statedims; ///<state dimensions array

    fri_float alpha;   ///<alpha value
    fri_float gamma;   ///<gamma value
    fri_float epsilon; ///<epsilon value
    fri_float qdiff_pos_boundary;    ///<upper boundary
    fri_float qdiff_neg_boundary;    ///<lower boundary
    fri_float qdiff_final_tolerance; ///<final tolerance
    fri_float reward_good_above;     ///<success if above this value
    fri_float rule_weight_considered_significant_for_update; ///<only update rules when their weights are above this value
    fri_float reduction_reward_tolerance;   ///<reward tolerance for reduction (e.g. 0.0 for exact match, INF if any reward difference is ok)
    unsigned char skip_rules;        ///<skip rules parameter - for backward compatibility (do not set it to 1 for new problems)
    unsigned char no_random;        ///<disable random action selection
    unsigned char construct_rb;        ///<construct rb (yes/no) - default 1
    unsigned char reduce_rb;        ///<reduce rb (yes/no) - default 0
    unsigned char reduction_strategy;    ///<reduction strategy (1, 2, !3, !4) - default 1
    int max_episodes; ///<maximum episodes
    int max_steps;    ///<maximum steps
    unsigned int episode_num;    ///<current episode no
    int five_maxnumofrules;

    ///reward callback function
    ///\param frirl initialized FRIRL struct
    ///\param[in] states current states array
    ///\param[in] states_len lengt of the states array
    ///\param[out] reward
    ///\return void
    void (* get_reward_func)(struct frirl_desc *frirl, fri_float *states, int states_len, struct frirl_reward_desc *reward);

    ///action callback function
    ///\param frirl initialized FRIRL struct
    ///\param[in] action
    ///\param[in] states current states array
    ///\param[in] states_len length of the states array
    ///\param[out] new_states
    ///\return void
    void (* do_action_func) (struct frirl_desc *frirl, fri_float action, fri_float *states, int states_len, fri_float *new_states);

    ///quantize observations callback function
    ///\param frirl initialized FRIRL struct
    ///\param[in] states current states array
    ///\param[in] states_len length of the states array
    ///\param[out] new_states
    ///\return void
    void (* quant_obs_func) (struct frirl_desc *frirl, fri_float *states, int states_len, fri_float *new_states);

//TODO check why is this commented
//#ifdef BUILD_VISUALIZATION
    ///draw callback function
    ///\param frirl initialized FRIRL struct
    ///\param[in] new_curr_states current states array
    ///\param[in] action - value
    ///\param[in] steps - current step no in episode
    ///\return void
    void (* draw_func) (struct frirl_desc *frirl, fri_float *new_curr_state, double action, unsigned int steps);

//#endif
    // "private" members, do not modify!
    struct FIVERB *fiverb; ///<instance of FIVE struct
    double *fiverb_ua;     ///<shortcut pointer to action U (in fiverb->u)
    double *fiverb_vea;    ///<shortcut pointer to action VE (in fiverb->ve)
    unsigned int reduction_state; ///<reduction state - if set to 1 no RB updates will be performed
//    int rules_len;        ///<rules num
    int numofantecedents; ///<antecedents length
    struct frirl_values_desc *possible_states;  ///<possible states array
    struct frirl_values_desc *possible_actions; ///<possible actions array
    struct frirl_reward_desc reward;            ///<reward

    fri_float *fgba_vagdist_states;   ///<to be used in get_best_action()
    fri_float *fgba_ruledist;         ///<to be used in get_best_action()
    fri_float *fgba_actconc;          ///<to be used in get_best_action()
    fri_float *fgba_dists;            ///<to be used in get_best_action()
    fri_float *fgba_statedistsum;     ///<to be used in get_best_action()
    fri_float *fus_proposed_values;   ///<to be used in update_sarsa() //obsoleted, rc5 and old API only
    fri_float *fus_values;            ///<to be used in update_sarsa() //obsoleted, rc5 and old API only
    fri_float *fus_check_states;      ///<to be used in update_sarsa() when CHECK_STATES build options is enabled
    fri_float fus_is_rule_inserted;   ///<to be used in update_sarsa() for skip_rules
    fri_float *fep_ant;       ///<to be used in episode()
    fri_float *fep_cur_ant;   ///<to be used in episode()
    fri_float *fep_q_ant;     ///<to be used in episode()
    fri_float *fep_cur_q_ant; ///<to be used in episode()

    int is_running;
    int agent_id;
    int agent_world_size;
    int epended;

    // imitation
    unsigned int keyaction;
    int valid_simulation;
    int original_learning;
    int user_exited;

};

//rc5 API
#define FRIRL frirl_desc

#endif //FRIRL_TYPES_H
