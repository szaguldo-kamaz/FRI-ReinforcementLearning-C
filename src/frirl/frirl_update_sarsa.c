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


//#define DEBUG

#include "frirl.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>


static inline void update_rules(struct frirl_desc *frirl, fri_float *values, fri_float qnow, fri_float qdiff) {

    struct FIVERB *frb = frirl->fiverb;

#ifdef BUILD_AVX2
    fri_float skiprulesave;
#endif

    int rules = frb->numofrules;
    if (frirl->fus_is_rule_inserted) {
        rules--;
    }

    /*if (frirl->rules_len < frb->numofrules)*/
//    printf("fus update_rules: %c numof: %d len: %d rules: %d\n", rule_found == ~0 ? '~' : '0', frb->numofrules, frirl->rules_len, rules);
    //%- FIVEVagConclWeight_fixres - only weight without the conclusion value!
//    unsigned int rule_hit_where=FIVEVagConclWeight(frirl->fiverb, values);
    fri_float *rules_to_update = frb->weights;
    unsigned int rule_hit_where = FIVE_vag_concl_weight(frirl->fiverb, values, rules_to_update);

    // As there are no NaN-s in the rule-base, when there's a direct rule hit, only one rule hits with 1.0
    // which means, that the current situation is on the grid of possible rule places
    // otherwise if it's not with 1.0 weight, the surrounding rules on the grid should be updated
    if ( (rule_hit_where != ~0) && // exact hit at rule no. "rule_hit_where"
        (
         ( frirl->skip_rules == 0 ) ||
         ( ( frirl->skip_rules == 1 ) && ( rule_hit_where < rules ) )
        )
       )
    {
        // cartpole comes here 25673x times
        // acrobot  comes here  2247x times

        frb->rconc[rule_hit_where] = qnow + qdiff; // * rules_to_update[rule_hit_where]; // as update is always 1.0 in this case, there's no need for the multiplication

        return;

    } else { // no exact hit

        if ( (frirl->skip_rules == 1 ) && (rule_hit_where == rules) ) {
            return;
        }

    }

//    fri_float *rules_to_update = frb->weights;

    // backward compatibility with MATLAB version
    if (frirl->skip_rules == 0) {
//        frirl->rules_len = frb->numofrules;
        rules = frb->numofrules;
        frirl->fus_is_rule_inserted = 0;
#ifdef BUILD_AVX2
    } else {
        skiprulesave = frb->rconc[frb->numofrules - 1];
#endif
    }

    // update Q for every rule

#ifdef BUILD_AVX2

    // cartpole   750x jar itt
    // acrobot  11026x jar itt

    // a frirl->skip_rules opcio miatt az utolso szabaly feluliras gond, mert mogotte olyankor mar van egy uj szabaly - skiprulesave valtozo megoldja

    asm volatile (
            "xor    %%ecx, %%ecx        \n" // i = 0
            "mov    %0, %%rax           \n" // &rules_to_update[0]
            "mov    %1, %%rbx           \n" // &frb->rconc[0]
            "vbroadcastsd %2, %%ymm1    \n" // ymm1 = < qnow qnow qnow qnow >
            "vbroadcastsd %3, %%ymm2    \n" // ymm2 = < qdiff qdiff qdiff qdiff >
            "vbroadcastsd %4, %%ymm3    \n" // ymm3 = < 4x frirl->rule_weight_considered_significant_for_update >
            "mov    %5, %%edx           \n" // edx = arrsize / frb->avx2_rbsize

    "fusl%=:   vmovdqa (%%rax), %%ymm0           \n" // ymm0 = rules_to_update[]
            "vcmppd $14, %%ymm3, %%ymm0, %%ymm4  \n" // ymm4 = _CMP_GT_OS ymm3 < ymm0 (ymm3 = 0.05 0.05 0.05 0.05)

            // depending on the problem, sometimes it could be possible that without test+jz the loop will run faster
            // but it's not worth... this whole loop loop is called relatively few times (cartpole: 750x ...)
            "vptest %%ymm4, %%ymm4               \n" // csak flageket valtoztat
            "jz fusski%=                         \n" // ugrik ha sehol sem volt igaz, szoval ymm1 > ymm0
                                                      // egyebkent meg kell csinalni a szamolast, folyt:
            "vmulpd  %%ymm0, %%ymm2, %%ymm0      \n" // ymm0 = rules_to_update[] * qdiff
            "vaddpd  %%ymm0, %%ymm1, %%ymm0      \n" // ymm0 = (rules_to_update[i] * qdiff) + qnow

            "vmovdqa (%%rbx), %%ymm6             \n" // ymm6 = frb->rconc[i]

            "vblendvpd %%ymm4, %%ymm0, %%ymm6, %%ymm6 \n" // ymm6 = ymm4 alapjan ymm6/ymm0 vegyitese - ahol alacsonyabb volt ott atkopizzuk az eredetit, a tobbihez az uj szamitas - ahol kisebb ott y
            "vmovdqa %%ymm6, (%%rbx)             \n" // frb->rconc[i] = ymm6

    "fusski%=: inc    %%ecx             \n"
            "add    $32, %%rax          \n"
            "add    $32, %%rbx          \n"
            "cmp    %%edx,%%ecx         \n"
            "jne    fusl%=              \n"

    : : "m" (rules_to_update), "m" (frb->rconc), "m" (qnow), "m" (qdiff), "m" (frirl->rule_weight_considered_significant_for_update), "m" (frb->avx2_rbsize) : "rax", "rbx", "rcx", "rdx", "ymm0", "ymm1", "ymm2", "ymm3", "ymm4", "ymm5", "ymm6", "memory" );

    // because of frirl->skip_rules "feature" - in case the new 1 rule was overwritten, then restore
//    if (frirl->rules_len < frb->numofrules) {
    if (frirl->fus_is_rule_inserted) {
        frb->rconc[frb->numofrules - 1] = skiprulesave;
    }

#else // no avx2

    for (unsigned int i = 0; i < rules; i++) {
        // cartpole 4966028x jar itt  - a foron belul - ebbol amikor van is frissites (0.05<x) ennyiszer: 34514
        DEBUG_MSG("frirl_update_sarsa: rule %d weight: %.18f\n", i, rules_to_update[i]);
        //DEBUG_MSG("frirl_update_sarsa: rule %d weight: %.13a\n",ruli,rules_to_update[i]);
        if (rules_to_update[i] > frirl->rule_weight_considered_significant_for_update) { //0.05
//            DEBUG_MSG("frirl_update_sarsa: updating rule %d from %.18f to %.18f\n", i, frb->rconc[i], qnow + qdiff * rulestoupdate[i]);
            DEBUG_MSG("frirl_update_sarsa: updating rule %d from %.13a to %.13a\n", i, frb->rconc[i], qnow + qdiff * rules_to_update[i]);

            frb->rconc[i] = qnow + qdiff * rules_to_update[i];
        }
    }
#endif // BUILD_AVX2

}


static inline void check_possible_states(struct frirl_desc *frirl, fri_float *ant, fri_float *newant) {

    unsigned int i;

    for (i = 0; i < frirl->statedims_len; i++) {
        //        [newrulestates(currps)]=frirl_check_possible_states(s(currps),possiblestates{currps},possiblestates_epsilons{currps});
        DEBUG_MSG("frirl_update_sarsa: currps: %d frirl->possiblestates[currps]->statedims: %d\n", i, frirl->possible_states[i].values_len);

        newant[i] = frirl_check_possible_states(frirl, ant[i], &frirl->possible_states[i]);

        DEBUG_MSG("frirl_update_sarsa: newrulestates[%d]: %f\n", i, new_rule_states[i]);

        if (newant[i] == INFINITY) {
            printf("FATAL: frirl_update_sarsa: CHKPOSSSTATEBUG: currps: %d\n", i);
            exit(2);
        }
    }

    newant[i] = frirl_check_possible_states(frirl, ant[i], frirl->possible_actions);
    if (newant[i] == INFINITY) {
        printf("FATAL: frirl_update_sarsa: POSSACTBUG\n");
        exit(2);
    }

}


////static inline unsigned char is_existing_rule(struct frirl_desc *frirl, fri_float *new_rule)
//static inline unsigned int is_existing_rule(struct frirl_desc *frirl, fri_float *new_rule)
//{
//    struct FIVERB *frb = frirl->fiverb;
//// cartpole 25907x jar itt
//
//    // check whether the supposed new rule exists
////    unsigned char rule_found = 0;
//    unsigned int rule_found = ~0;
//
//    for (unsigned int i = 0; i < frb->numofrules; i++) { // cartpole 1139193x jar itt bent a ciklusban
//        unsigned int ci;
//
//        DEBUG_MSG("frirl_update_sarsa: comparing antecedents of rule %d.\n", i);
//
//        // only the antecedens should be checked
//        for (ci = 0; ci < frb->numofantecedents; ci++) {
//            DEBUG_MSG("%f = %f, ", new_rule[ci], (frb->rant + i * frb->numofantecedents)[ci]);
//
//            if (new_rule[ci] != (frb->rant + i * frb->numofantecedents)[ci])
//                break;
//        }
//
//        DEBUG_MSG("\n");
//        DEBUG_MSG("frirl_update_sarsa: ci: %d frirl->frb->rulelength-1: %d\n", ci, frb->rulelength - 1);
//
//        if (ci == frb->numofantecedents) {
////            rule_found = 1;
//            rule_found = i;
//
//            DEBUG_MSG("frirl_update_sarsa: rulefound=%d\n",i);
//
//            break;
//        }
//    }
//
//    return rule_found;
//}


//static inline void insert_rule(struct frirl_desc *frirl, fri_float *values, int values_len, 
//    fri_float *proposed_values, fri_float qnow, fri_float qdiff)
//{
////    fri_float qnew_rule         = FIVEVagConcl(frirl->fiverb, proposed_values); // % Q(s,a)
//    fri_float qnew_rule;
//    unsigned int rule_found = FIVE_vag_concl(frirl->fiverb, proposed_values, &qnew_rule);
////    proposed_values[values_len] = qnew_rule + qdiff; // conclusion -> rconc[last]
////    fri_float *new_rule         = proposed_values; // full rule with q value
//
//    DEBUG_MSG("frirl_update_sarsa: qnewrule: %.28f qnewrule+qdiff: %.28f\n", qnew_rule, qnew_rule + qdiff);
//
//// check whether the new proposed rule exists
////    unsigned char rule_found = is_existing_rule(frirl, new_rule);
////    unsigned int rule_found = is_existing_rule(frirl, new_rule);
////    unsigned int rule_found = is_existing_rule(frirl, proposed_values);
//
//// if not, then add it
////    if (rule_found == 0) {
//    if (rule_found == ~0) {
//        // append new rule to the end of the existing rulebase
//        frirl->fus_is_rule_inserted = 1;
////        five_add_rule(frirl->fiverb, new_rule);
//        FIVE_add_rule(frirl->fiverb, proposed_values, qnew_rule + qdiff);
//
////        DEBUG_MSG("frirl_update_sarsa: newrule q: %.18f\n", new_rule[frb->rulelength - 1]);
//        DEBUG_MSG("frirl_update_sarsa: newrule q: %.18f\n", qnew_rule + qdiff);
//    } else {
//// if exists, then update the Q values of existing rules based on the current situation (not necessarily only the existing rule! - because the current situation is not limited to the grid points, whereas new_rule can be on predefined grid points)
//        // or update rule 'rno'
////        update_rules(frirl, values, qnow, qdiff);
//        frirl->fus_is_rule_inserted = 0;
////        update_rules(frirl, values, qnow, qdiff, rule_found);
//        update_rules(frirl, values, qnow, qdiff);
//    }
//}

///update sarsa (rc5 API, kept for backward compatibility)
///\param frirl initialized frirl struct
///\param s states
///\param a action
///\param r reward
///\param sp proposed states
///\param ap proposed action
///\return void
void frirl_update_sarsa(struct frirl *frirl, double *s, double a, double r, double *sp, double ap) {

    unsigned int i;
    unsigned int values_len = frirl->statedims_len + 1;
    fri_float *q_ant = frirl->fus_values;
    fri_float *cur_q_ant = frirl->fus_proposed_values;

    for (i = 0; i < values_len - 1; i++) {
        q_ant[i] = s[i];
        cur_q_ant[i] = sp[i];

        DEBUG_MSG("frirl_update_sarsa: values[%d]: %.18f proposed_values[%d]: %.18f ", i, q_ant[i], i, cur_q_ant[i]);
    }
    q_ant[values_len - 1] = a;
    cur_q_ant[values_len - 1] = ap;

    DEBUG_MSG("\nfrirl_update_sarsa: values[%d]: %.18f proposed_values[%d]: %.18f\n", i, q_ant[i], i, cur_q_ant[i]);

    frirl_update_sarsa(frirl, q_ant, r, cur_q_ant);

}


///update sarsa (old API, to be removed)
///\param frirl initialized frirl struct
///\param states
///\param action
///\param reward
///\param proposed_states
///\param proposed_action
///\return void
void frirl_update_sarsa(struct frirl_desc *frirl, fri_float *states, fri_float action,
    fri_float reward, fri_float *proposed_states, fri_float proposed_action) {

    frirl_update_sarsa(frirl, states, action, reward, proposed_states, proposed_action);

}


/*void frirl_update_sarsa(struct frirl_desc *frirl, fri_float *states, fri_float action,
    fri_float reward, fri_float *proposed_states, fri_float proposed_action)
{        // frirl_update_sarsa: frirl-learning framework: Update the rule-base

    // cartpole    33002x jar itt
    // acrobot     21207x jar itt
    // mountaincar 15548x jar itt

    unsigned int i;
    unsigned int values_len = frirl->statedims_len + 1;
    fri_float *values = frirl->fus_values;
    fri_float *proposed_values = frirl->fus_proposed_values;

    for (i = 0; i < values_len - 1; i++) {
        values[i] = states[i];
        proposed_values[i] = proposed_states[i];

        DEBUG_MSG("frirl_update_sarsa: values[%d]: %.18f proposed_values[%d]: %.18f ", i, values[i], i, proposed_values[i]);
    }
    values[values_len - 1] = action;
    proposed_values[values_len - 1] = proposed_action;

    DEBUG_MSG("\nfrirl_update_sarsa: values[%d]: %.18f proposed_values[%d]: %.18f\n", i, values[i], i, proposed_values[i]);

    fri_float qp    = FIVEVagConcl(frirl->fiverb, proposed_values); //% Q(sp,ap)
    fri_float qnow  = FIVEVagConcl(frirl->fiverb, values); //% Q(s,a)
    fri_float qdiff = frirl->alpha * (reward + frirl->gamma * qp - qnow);

    DEBUG_MSG("frirl_update_sarsa: qnow: %.18f qp: %.18f qdiff: %.18f\n", qnow, qp, qdiff);
    DEBUG_MSG("frirl_update_sarsa: qnow: %.13a qp: %.13a qdiff: %.13a\n", qnow, qp, qdiff);

    if ((qdiff > frirl->qdiff_pos_boundary) || (qdiff < frirl->qdiff_neg_boundary)) {
//        frirl->rules_len = frb->numofrules;
        // insert new rule if it does not exist
        get_new_rule_states(frirl, states, values_len - 1, proposed_values);
        get_new_rule_action(frirl, action, &proposed_values[values_len - 1]);
        insert_rule(frirl, values, values_len, proposed_values, qnow, qdiff);
    } else {
        // if deltaq is not so big then just update existing rules (q value)
//        update_rules(frirl, values, qnow, qdiff);
        update_rules(frirl, values, qnow, qdiff, ~0);    // possibly many rules have to be updated with different weights (~0 means that no matching rule is in the rule-base and every rule has to be checked one by one - anyway this was the default before)
    }
 }*/


///update sarsa
///\param frirl initialized FRIRL struct
///\param q_ant (before doaction) the quantized antecedents
///antecedents means: <state0,state1,...,stateN,action>
///\param reward
///\param cur_q_ant (after doaction) the current quantized antecedents
///\return void
void frirl_update_sarsa(struct frirl_desc *frirl, fri_float *q_ant, fri_float reward, fri_float *cur_q_ant) {

    struct FIVERB *frb = frirl->fiverb;
    fri_float qp;
    fri_float qnow;
    fri_float qdiff;
    unsigned int rule_i;

    FIVE_vag_concl(frb, cur_q_ant, &qp);  //% Q(sp,ap)
    rule_i = FIVE_vag_concl(frb, q_ant, &qnow);  //% Q(s,a)
    qdiff = frirl->alpha * (reward + frirl->gamma * qp - qnow);

    DEBUG_MSG("frirl_update_sarsa: qnow: %.18f qp: %.18f qdiff: %.18f\n", qnow, qp, qdiff);
    DEBUG_MSG("frirl_update_sarsa: qnow: %.13a qp: %.13a qdiff: %.13a\n", qnow, qp, qdiff);

    if ((qdiff > frirl->qdiff_pos_boundary) || (qdiff < frirl->qdiff_neg_boundary)) {
        fri_float *rant = q_ant;
        fri_float rconc = qnow;

        if (CHECK_STATES) {
            rant = frirl->fus_check_states;
            check_possible_states(frirl, q_ant, rant);
            rule_i = FIVE_vag_concl(frb, rant, &rconc);
        }

        if (rule_i == ~0) {
            frirl->fus_is_rule_inserted = 1;
            FIVE_add_rule(frb, rant, rconc + qdiff);
            goto exit; // danger zone
        }
        frirl->fus_is_rule_inserted = 0;
    }
    update_rules(frirl, q_ant, qnow, qdiff);

exit:
    return;

}
