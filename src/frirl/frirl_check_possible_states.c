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
#include <math.h>
#include <stdlib.h>
#ifdef DEBUG
#include <stdio.h>
#endif


static inline void increment_place(struct frirl_values_desc *possible_states, int i, fri_float p_state_diff) {

    int pi;
    fri_float *new_possible_states;

    //possiblestate->statedims++;
    new_possible_states = MALLOC(sizeof(fri_float) * ++possible_states->values_len);

    for (pi = 0; pi < i; pi++) {
        new_possible_states[pi] = possible_states->values[pi];
    }
    new_possible_states[pi] = possible_states->values[i] + (p_state_diff / 2);

    //check why?
//        for (; pi < i; pi++) {
//            new_possible_statedim[pi + 1] = possible_states->states[pi];
//        }

#ifdef DEBUG
    for (pi = 0; pi < i; pi++) {
        printf("%f \n", new_possible_states[pi]);
    }
    printf("\n");
#endif

    //!!        possiblestate=[possiblestate(1:psno), possiblestate(psno)+pstatediff/2, possiblestate((psno+1):end)];  // MATLAB code
    free(possible_states->values);
    possible_states->values = new_possible_states;

}


static inline fri_float hit_between_possible_places(struct frirl_values_desc *possible_states, fri_float observation, int i) {

    // hit before first rule ?
    if (observation < possible_states->values[0]) {
        return possible_states->values[0]; 
    }

    i--;

    fri_float p_state_diff         = fabs(possible_states->values[i + 1] - possible_states->values[i]);  // distance between the two adjacent points
    fri_float p_state_center_start = p_state_diff / 4 + possible_states->values[i];  // where to start looking for new rule center ( found state + 1/4 of distance = 1/4 from center (left))
    fri_float p_state_center_stop  = p_state_diff / 4 * 3 + possible_states->values[i];  // where to stop looking for new rule center ( found state + 3/4 of distance = 1/4 from center (right))

    DEBUG_MSG("frirl_check_possible_states: pstatediff: %f pstatecenterstart: %f pstatecenterstop: %f\n", p_state_diff, p_state_center_start, p_state_center_stop);

    if ( (possible_states->epsilon <= (p_state_diff / 2)) &&
         (p_state_center_start <= observation) &&
         (observation <= p_state_center_stop) ) {
        // new possible rule place should be added

        increment_place(possible_states, i, p_state_diff);

        return possible_states->values[i + 1];
    }

    // close to a possible rule place, select the closer one
    fri_float relativeobs = observation - possible_states->values[i];
    fri_float relativeobs_next = possible_states->values[i + 1] - observation;

    // previous one is closer
    if (relativeobs < relativeobs_next)
        return possible_states->values[i];

    // next one is closer
    return possible_states->values[i + 1];
}


///check possible states
///\param frirl initialized FRIRL struct
///\param observation
///\param possible_states
///\return rule
fri_float frirl_check_possible_states(struct frirl_desc *frirl, fri_float observation, struct frirl_values_desc *possible_states) {

    // check possible rule places for state
    int is_place_found = 0;

    // search among possible rule places
    int i;
    for (i = 0; i < possible_states->values_len; i++) {
        if (observation < possible_states->values[i]) {
            is_place_found = 1;
            break; // found
        }
    }

    DEBUG_MSG("frirl_check_possible_states: is_place_found: %d psno: %d obs: %f possiblestate->statedim[0]: %f\n", is_place_found, i, observation, possible_states->values[0]);

    if (is_place_found == 1) {
        //hit between possible rule places
        return hit_between_possible_places(possible_states, observation, i);

    }

    // did not hit, and went over the last possible rule place
    DEBUG_MSG("frirl_check_possible_states: is_place_found=0: \"newrule\": %f\n", possible_states->values[possible_states->values_len - 1]);
    return possible_states->values[possible_states->values_len - 1];

}
