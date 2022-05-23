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
#include <stdlib.h>
#include <math.h>

///E Greedy selection
///\param frirl initialized FRIRL struct
///\param states
///\return selected action index
unsigned int frirl_e_greedy_selection(struct frirl_desc *frirl, fri_float *states) {
    // e_greedy_selection selects an action using Epsilon-greedy strategy

    if ((frirl->no_random == 1) || (frirl->epsilon == 0.0)) {
        return frirl_get_best_action(frirl, states);
    }

    double randres = (double) rand() / (double) RAND_MAX;
    if (randres > frirl->epsilon) {
        return frirl_get_best_action(frirl, states);
    }

    randres = round((double) rand() / (double) RAND_MAX * frirl->actiondim.values_len);

    return randres;

}
