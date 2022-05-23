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


#ifndef FRIRL_TYPES_DEF_H
#define FRIRL_TYPES_DEF_H

#include "frirl_types.h"


///default values for frirl_desc
static const struct frirl_desc frirl_desc_default = {
    .argc = 1,
    .argv = 0,
    .runmode = FRIRL_SEQ,
    .rbfile = 0,
    .visualization = 0,
    .gui_width = 640,
    .gui_height = 480,
    .verbose = 1,
    .agent_rnd_init = 1,

    .alpha   = 0.5,
    .gamma   = 1.0,
    .epsilon = 0.001,

    .qdiff_pos_boundary    = 1.0,
    .qdiff_neg_boundary    = -250.0,
    .qdiff_final_tolerance = 250.0,

    .reward_good_above = 0.0,
    .rule_weight_considered_significant_for_update = 0.05,
    .reduction_reward_tolerance = 0.0,
    .skip_rules = 0,    // only for backward compatibility, should be zero with new problems!
    .no_random  = 1,
    .construct_rb = 1,
    .reduce_rb = 0,
    .reduction_strategy = FRIRL_REDUCTION_STRATEGY_DEFAULT,

    .max_episodes = 1000,
    .max_steps    = 1000,

    .five_maxnumofrules = 16384,

    .get_reward_func = NULL,
    .do_action_func  = NULL,
    .quant_obs_func  = NULL,
//#ifdef BUILD_VISUALIZATION
    .draw_func       = NULL,
//#endif
    .numofantecedents = 0,
    .reduction_state = 0,
    .statedims_len = 0,
//    .rules_len = 0,
    .fus_is_rule_inserted = 0,

    .is_running = 0,
    .epended = 0,
    .agent_id = 0,

    // imitation
    .keyaction = -1,
    .valid_simulation = 0,
    .original_learning = 1,
    .user_exited = 0,

};

#endif //FRIRL_TYPES_DEF_H
