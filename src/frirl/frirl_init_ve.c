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
#include <stdlib.h>
#include <string.h>
#ifdef DEBUG
#include <stdio.h>
#endif

///generate FIVE params
///\param frirl initialized FRIRL struct
///\return 0 if success
int frirl_init_ve(struct frirl_desc *frirl, fri_float *ve, fri_float *u, int univlength)
{
    // frirl_gen_five_fri_params: Generate parameters for FRIRL framework

    int st, i;
    int usize = frirl->statedims_len * univlength;
    double *scf = MALLOC(usize * sizeof(fri_float));
    double *vetmp;
    double *SPact;
    double divratio;

    for (st = 0; st < frirl->statedims_len; st++) {

        double *curSP;
        unsigned int curstatedim;

        DEBUG_MSG("%d. %d db\n", st, frirl->statedims[st].values_len);

        //    steepness = statedims./2; // from MATLAB

        curSP = MALLOC(sizeof(double)*frirl->statedims[st].values_len * 3);

        for (curstatedim = 0; curstatedim < frirl->statedims[st].values_len; curstatedim++) {
            curSP[curstatedim * 3] = frirl->statedims[st].values[curstatedim];
            curSP[curstatedim * 3 + 1] = frirl->statedims[st].values_steep;  // states_steepness[st];
            curSP[curstatedim * 3 + 2] = frirl->statedims[st].values_steep;  // states_steepness[st];

            DEBUG_MSG("%f %f %f\n", curSP[curstatedim * 3], curSP[curstatedim * 3 + 1], curSP[curstatedim * 3 + 2]);
        }

//        scf = FIVEGScFunc(u + st * univlength, 1, univlength, curSP, frirl->statedims[st].values_len, 3, NAN);
        FIVE_GSc_func(u + st * univlength, 1, univlength, curSP, frirl->statedims[st].values_len, 3, NAN, scf);

#ifdef DEBUG
        for (i = 0; i < univlength; i++) {
            printf("%d. %f\n", i, scf[i]);
        }
#endif

        free(curSP);
        vetmp = FIVEGVagEnv(u + st * univlength, 1, univlength, scf);

        DEBUG_MSG("%x %x %x %d\n", ve + st * univlength, vetmp, univlength * sizeof(double), univlength * sizeof(double));

        memcpy(ve + st * univlength, vetmp, univlength * sizeof(double));
        free(vetmp);

#ifdef DEBUG
        for (i = 0; i < univlength; i++) {
            printf("%d. %f\n", i, ve[i]);
        }

        printf("%d. min: %f max: %f\n", st, minvals[st], maxvals[st]);
#endif

    }

    // generate actions partition

    SPact = MALLOC(sizeof(double)*frirl->actiondim.values_len * 3);
    divratio = 1.0 / (frirl->actiondim.values_len - 1) * 2.0;

    DEBUG_MSG("mintmp: %f maxtmp: %f\n", mintmp, maxtmp);

    for (i = 0; i < frirl->actiondim.values_len; i++) {
        SPact[i * 3] = i * divratio - 1.0;
        SPact[i * 3 + 1] = (frirl->actiondim.values_len - 1) / 2;
        SPact[i * 3 + 2] = (frirl->actiondim.values_len - 1) / 2;

        DEBUG_MSG("%d. %f %f %f\n", i, SPact[i * 3], SPact[i * 3 + 1], SPact[i * 3 + 2]);
    }

//    scf = FIVEGScFunc(u + frirl->statedims_len * univlength, 1, univlength, SPact, frirl->actiondim.values_len, 3, NAN);
    FIVE_GSc_func(u + frirl->statedims_len * univlength, 1, univlength, SPact, frirl->actiondim.values_len, 3, NAN, scf);

#ifdef DEBUG
    for (i = 0; i < univlength; i++) {
        printf("%d. %f\n", i, scf[i]);
    }
#endif

    vetmp = FIVEGVagEnv(u + frirl->statedims_len * univlength, 1, univlength, scf);
    memcpy(ve + frirl->statedims_len * univlength, vetmp, univlength * sizeof(double));

#ifdef DEBUG
    for (i = 0; i < univlength; i++) {
        printf("%d. %f\n", i, vetmp[i]);
    }
#endif

    free(scf);
    free(SPact);
    free(vetmp);

    return 0;

}
