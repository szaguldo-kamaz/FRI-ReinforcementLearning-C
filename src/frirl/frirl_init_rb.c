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


static void get_min_max_arr(struct frirl_desc *frirl, fri_float *minvals, fri_float *maxvals) {
    // frirl_gen_five_fri_params: Generate parameters for FRIRL framework

    int st, i;
    double maxtmp, mintmp;

    for (st = 0; st < frirl->statedims_len; st++) {

        double *curSP;
        unsigned int curstatedim;

        DEBUG_MSG("%d. %d db\n", st, frirl->statedims[st].values_len);

        //    steepness=statedims./2;

        maxtmp = frirl->statedims[st].values[0];
        mintmp = maxtmp;
        for (curstatedim = 0; curstatedim < frirl->statedims[st].values_len; curstatedim++) {

            if (frirl->statedims[st].values[curstatedim] > maxtmp) {
                maxtmp = frirl->statedims[st].values[curstatedim];
            }
            if (frirl->statedims[st].values[curstatedim] < mintmp) {
                mintmp = frirl->statedims[st].values[curstatedim];
            }

            DEBUG_MSG("%f %f %f\n", curSP[curstatedim * 3], curSP[curstatedim * 3 + 1], curSP[curstatedim * 3 + 2]);
        }
        minvals[st] = mintmp;
        maxvals[st] = maxtmp;

#ifdef DEBUG
        printf("%d. min: %f max: %f\n", st, minvals[st], maxvals[st]);
#endif

    }

    // generate actions partition
    mintmp = frirl->actiondim.values[0];
    maxtmp = frirl->actiondim.values[0];

    DEBUG_MSG("mintmp: %f maxtmp: %f\n", mintmp, maxtmp);

    for (i = 0; i < frirl->actiondim.values_len; i++) {

        if (frirl->actiondim.values[i] > maxtmp) {
            maxtmp = frirl->actiondim.values[i];
        }
        if (frirl->actiondim.values[i] < mintmp) {
            mintmp = frirl->actiondim.values[i];
        }
    }
    minvals[frirl->statedims_len] = mintmp;
    maxvals[frirl->statedims_len] = maxtmp;

    DEBUG_MSG("mintmp: %f maxtmp: %f\n", mintmp, maxtmp);
    DEBUG_MSG("mina: %f maxa: %f\n", minvals[frirl->statedims_len], maxvals[frirl->statedims_len]);

}

///generate initial rulebase
///\param frirl initialized frirl struct
///\param max_values maximum values
///\param min_values minimum values
///\return 0 if success
int frirl_init_rb(struct frirl_desc *frirl, fri_float *rant, fri_float *rconc, int *numofrules)
{

    fri_float value;
    unsigned int rules, divider;

    fri_float *maxvals = MALLOC(sizeof(double) * frirl->numofantecedents);
    fri_float *minvals = MALLOC(sizeof(double) * frirl->numofantecedents);


    get_min_max_arr(frirl, minvals, maxvals);

    //rules=pow(2,frirl->nantecedents);
    *numofrules = pow(2, frirl->numofantecedents);

    //nem kell mert ugyis most adunk erteket neki
//    for (int i = 0; i < frirl->numofantecedents * FRIRLRULES; i++) {
////        frirl->fiverb->rb[i] = 0.0;
//        frirl->fiverb->rant[i] = 0.0;
//    }
    //eleg annyi konkluziot kinullazni amennyi rule van, uj szabaly hozzadasaval ugyis kap ertek
//    for (int i = 0; i < FRIRLRULES; i++) {
    for (int i = 0; i < *numofrules; i++) {
        rconc[i] = 0.0;
    }

    for (int i = 0; i < frirl->numofantecedents; i++) {

        divider = *numofrules >> (i + 1);

        for (int j = 0; j <= *numofrules - 1; j++) {

            if (((j / divider) % 2) == 0)
                value = minvals[i];
            else
                value = maxvals[i];

//            frirl->fiverb->rb[j * (frirl->numofantecedents + 1) + i] = value;
            rant[j * frirl->numofantecedents + i] = value;
        }
    }

    // not so nice (these (rant, rconc, etc.) were allocated by FIVEInit)
/*FIVEInit handles this later well
    free(frirl->fiverb->rant);
    frirl->fiverb->rant = MALLOC(FRIRLRULES * (frirl->fiverb->rulelength - 1) * sizeof(fri_float));

    free(frirl->fiverb->rconc);
    frirl->fiverb->rconc = MALLOC(FRIRLRULES * sizeof(fri_float));

    free(frirl->fiverb->ruledists);
    frirl->fiverb->ruledists = MALLOC(FRIRLRULES * sizeof(fri_float));

    free(frirl->fiverb->weights);
    frirl->fiverb->weights = MALLOC(FRIRLRULES * sizeof(fri_float));
*/

    free(minvals);
    free(maxvals);

    return 0;
}
