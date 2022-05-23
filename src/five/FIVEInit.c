/*
** Fuzzy Inference by Interpolation in Vague Environment toolbox for ANSI C
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** ANSI C / AVX port and FixRes optimization by David Vincze <david.vincze@webcode.hu>
** Various contributions by Daniel Palko <palko.daniel@uni-miskolc.hu>
** mex port by Sandor Zsuga (zsuga@iit.uni-miskolc.hu)
** Original FIVE MATLAB code by Szilveszter Kovacs <szkovacs@iit.uni-miskolc.hu>
*/

///
/// \file FIVEInit.c
///

#include "FIVE.h"
#include <stdlib.h>
#include <stdio.h>
#include "../inl/fast_abs.inl"
#include "../inl/arr.inl"
#include "../inl/min.inl"

///FIVE initialization
///\param u universe arrays
///\param ve vague environment arrays
///\param p parameter for the Shepard interpolation
///\param numofunivs dimension count of the universe array (number of antecedents + conclusion in case it is also fuzzy, usually onbly no., of antec)
///\param univlength length of one universe (resolution of one dimension, number of elements in one "u")
///\param numofrules number of rules (number of vectors in *rb)
///\param rulelength number of rule dimensions (number of elements in each vector (one rule) in *rb)
///\param rulebase the rule-base array
///\return initialized FIVE struct

/*
e.g.:
RB:
Ant1    Ant2    Ant3    Ant4    Conc
1       1       0       1       0.2
0.5     0.2     0.5     1       1
0       0.5     0       0.6     0.8

then:
mr/numofrules = 3 (3 rules)
nr/rulelength = 5 (5 dimensions, 4 ant + 1 conc)

U:
0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0
0.0 0.05  0.1  0.15 0.2  0.25 0.3  0.35 0.4  0.45 0.5
0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0
Then
nu/univlength: 11 (11 elements in one "u" vector)
mu/numofunivs: 3 (3 universes)
 */

struct FIVERB *FIVEInit(double *u, double *ve, int p, int numofunivs, int univlength, int numofrules, int maxnumofrules, int rulelength, double *rant, double *rconc) {

    struct FIVERB *frb;
    int l = 0, j, k;


    if (numofunivs > FIVE_MAX_NUM_OF_UNIVERSES) {
        printf("FIVEInit: fatal error: given number of universes greater than FIVE_MAX_NUM_OF_UNIVERSES (%d) - increase it in the source code\n!", FIVE_MAX_NUM_OF_UNIVERSES);
        exit(1);
    }

#ifdef DEBUG
    printf("FIVEInit: Preparing FRB\n");
#endif

    frb = MALLOC(sizeof(struct FIVERB));
    if (frb == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(FIVERB)\n");
#endif
        return NULL;
    }

    frb->u = u;
    frb->ve = ve;
    frb->numofunivs = numofunivs;
    frb->univlength = univlength;
    frb->numofrules = numofrules;
    frb->maxnumofrules = maxnumofrules;
#ifdef BUILD_AVX2
    frb->avx2_rbsize = maxnumofrules / 4 + (maxnumofrules % 4 > 0);
#endif
    frb->rulelength = rulelength;
    frb->numofantecedents = rulelength - 1;
    if (p == 0) {
        frb->p = frb->rulelength - 1;
    } else {
        frb->p = p;
    }
//    frb->rb = (double *) rulebase;

#ifdef DEBUG
    printf("FIVEInit: Preparing FRB - malloc()\n");
#endif

    frb->vek = MALLOC(frb->numofunivs * sizeof(double *));
    if (frb->vek == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(vek)\n");
#endif
        return NULL;
    }

    frb->uk = MALLOC(frb->numofunivs * sizeof(double *));
    if (frb->uk == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(uk)\n");
#endif
        return NULL;
    }

    for (unsigned char ui = 0; ui < frb->numofunivs; ui++) {
        frb->uk[ui]  = frb->u  + ui*frb->univlength;
        frb->vek[ui] = frb->ve + ui*frb->univlength;
    }

    frb->ruledists = MALLOC(frb->maxnumofrules * sizeof(double));
    if (frb->ruledists == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(rd)\n");
#endif
        return NULL;
    }

    frb->rant = (double *)rant;
//    frb->rant = MALLOC(frb->numofrules * (frb->rulelength - 1) * sizeof(double));
//    if (frb->rant == NULL) {
//#ifdef DEBUG
//        printf("FIVEInit: Error - malloc(rant)\n");
//#endif
//        return NULL;
//    }

    frb->rant_uindex = MALLOC(frb->maxnumofrules * (frb->rulelength - 1) * sizeof(unsigned int));
    if (frb->rant_uindex == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(rant_uindex)\n");
#endif
        return NULL;
    }

    frb->rant_veval = MALLOC(frb->maxnumofrules * (frb->rulelength - 1) * sizeof(double));
    if (frb->rant_veval == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(rant_veval)\n");
#endif
        return NULL;
    }

    frb->rseqant = MALLOC((frb->rulelength - 1) * sizeof(double *));
    if (frb->rseqant == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(rseqant)\n");
#endif
        return NULL;
    }

    frb->rseqant_uindex = MALLOC((frb->rulelength - 1) * sizeof(unsigned int *));
    if (frb->rseqant_uindex == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(rseqant_uindex)\n");
#endif
        return NULL;
    }

    frb->rseqant_veval = MALLOC((frb->rulelength - 1) * sizeof(double *));
    if (frb->rseqant_veval == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(rseqant_veval)\n");
#endif
        return NULL;
    }

    for (unsigned int currant = 0; currant < (frb->rulelength - 1); currant++) {
        frb->rseqant[currant] = MALLOC(frb->maxnumofrules * sizeof(fri_float));
        frb->rseqant_uindex[currant] = MALLOC(frb->maxnumofrules * sizeof(unsigned int));
        frb->rseqant_veval[currant] = MALLOC(frb->maxnumofrules * sizeof(fri_float));
#ifdef DEBUG
        printf("frb->rseqant[%d]: %lx\n",currant,frb->rseqant[currant]);
        printf("frb->rseqant_uindex[%d]: %lx\n",currant,frb->rseqant_uindex[currant]);
        printf("frb->rseqant_veval[%d]: %lx\n",currant,frb->rseqant_veval[currant]);
#endif
        if (frb->rseqant[currant] == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(rseqant[%d])\n",currant);
#endif
        return NULL;
        }
        if (frb->rseqant_uindex[currant] == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(rseqant_uindex[%d])\n",currant);
#endif
        return NULL;
        }
        if (frb->rseqant_veval[currant] == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(rseqant_veval[%d])\n",currant);
#endif
        return NULL;
        }
        memset(frb->rseqant[currant], 0, frb->maxnumofrules * sizeof(fri_float));
        memset(frb->rseqant_veval[currant], 0, frb->maxnumofrules * sizeof(fri_float));
        memset(frb->rseqant_uindex[currant], 0, frb->maxnumofrules * sizeof(unsigned int));
    }

//TODO remove
frb->epno=0;

    // just to make it more comfortable to access the sequential action dimensions
    frb->ract=frb->rseqant[frb->rulelength - 2];
    frb->ract_uindex=frb->rseqant_uindex[frb->rulelength - 2];
    frb->ract_veval=frb->rseqant_veval[frb->rulelength - 2];

    frb->rconc = (double *)rconc;
//    frb->rconc = MALLOC(frb->numofrules * sizeof(fri_float));
//    if (frb->rconc == NULL) {
//#ifdef DEBUG
//        printf("FIVEInit: Error - malloc(rconc)\n");
//#endif
//        return NULL;
//    }
//    memset(frb->rconc, 0, frb->numofrules * sizeof(fri_float));

    frb->ukdomains = MALLOC(frb->numofunivs * sizeof(fri_float));
    if (frb->ukdomains == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(ukdomains)\n");
#endif
        return NULL;
    }

    frb->udivs = MALLOC(frb->numofunivs * sizeof(fri_float));
    if (frb->udivs == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(udivs)\n");
#endif
        return NULL;
    }

    frb->uksize = frb->univlength - 1;
    for (unsigned char k = 0; k < frb->numofunivs; k++) {
        frb->ukdomains[k] = frb->u[k*frb->univlength + frb->uksize] - frb->u[k*frb->univlength];
        frb->udivs[k] = frb->ukdomains[k] / (frb->uksize);
    }

#ifdef DEBUG
    printf("FIVEInit: Preparing rant and rconc\n");
#endif

    // R must be chopped to contain only the antecedents
    // also extract conclusions from R (rule-base) into rconc

    //eleg csak numofrules (maxnumofrules helyett), mert five_add_rule ugyis beallitja ezeket
    for (j = 0; j < frb->maxnumofrules; j++) {
        for (k = 0; k < frb->rulelength - 1; k++) {
            //rb removed
//            frb->rant[l] = frb->rseqant[k][j] = frb->rb[j * frb->rulelength + k];
//            frb->rant[l] = rant[l];
            frb->rseqant[k][j] = frb->rant[l];
            frb->rant_uindex[l] = frb->rseqant_uindex[k][j] = get_vag_abs_min_i_fixres(frb->u + (frb->univlength) * k, frb->univlength, frb->rseqant[k][j], frb->udivs[k]);
            frb->rant_veval[l] = frb->rseqant_veval[k][j] = frb->vek[k][frb->rseqant_uindex[k][j]];
#ifdef DEBUG
            printf("ruleno: %d k: %d u[0]: %lf len: %d, antval[k][j]: %lf udiv[k]: %lf = %d\n", j, k, *(frb->u + (frb->univlength) * k), frb->univlength, frb->rseqant[k][j], frb->udivs[k], frb->rant_uindex[l]);
#endif
            l++; // l = (j * frb->rulelength - 1) + k
#ifdef DEBUG
            printf("frb->rb[%d]: %f\n", j * frb->rulelength + k, frb->rb[j * frb->rulelength + k]);
#endif
        }
        //rb removed
//        frb->rconc[j] = rconc[j];//frb->rb[j * frb->rulelength + k];
#ifdef DEBUG
        printf("frb->rconc[%d]: %a = frb->rb[%d]: %a\n", j, frb->rconc[j], j * frb->rulelength + k, frb->rb[j * frb->rulelength + k]);
#endif
    }

    // TODO rb removed
//    frb->newrule = frb->rb + (frb->numofrules) * frb->rulelength;
    frb->newrconc = frb->rconc + frb->numofrules;
    frb->newrant = frb->rant + frb->numofrules * (frb->rulelength - 1);

    // to be used in VagConclWeight
#ifdef BUILD_AVX2
    frb->wi = MALLOC(sizeof(double)*(frb->maxnumofrules+4));
#else
    frb->wi = MALLOC(sizeof(double)*(frb->maxnumofrules+1));
#endif
    if (frb->wi == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(wi)\n");
#endif
        return NULL;
    }

    // to be used in RuleDist
    frb->frd_dists = (fri_float *)MALLOC((frb->maxnumofrules + 4) * sizeof(fri_float) * FIVE_MAX_NUM_OF_UNIVERSES);
    if (frb->frd_dists == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(frd_dists)\n");
#endif
        return NULL;
    }
    memset(frb->frd_dists, 0, (frb->maxnumofrules + 4) * sizeof(fri_float) * FIVE_MAX_NUM_OF_UNIVERSES);

    // to be used in VagConcl
    frb->fvc_vagdist = MALLOC(frb->numofunivs * sizeof(double));
    if (frb->fvc_vagdist == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(fvc_vagdist)\n");
#endif
        return NULL;
    }

    // to be used in ValVag
    frb->valvagp = MALLOC(frb->numofunivs * sizeof(double)); // ValVag output point value, in fact only 1 element will be used, because it is only called for the consequent dim if it has a VE
    if (frb->valvagp == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(valvagp)\n");
#endif
        return NULL;
    }

// left for reference only
//    frb->valvagu  = frb->u  + (frb->numofunivs - 1) * frb->univlength; // consequent U pointer
//    frb->valvagve = frb->ve + (frb->numofunivs - 1) * frb->univlength; // consequent VE pointer
    frb->valvagu  = frb->uk [frb->numofunivs - 1]; // consequent U pointer
    frb->valvagve = frb->vek[frb->numofunivs - 1]; // consequent VE pointer
    frb->valvagdims = 1; // FIVEValVag is only called when the conclusion has a VE, and only for the conclusion dim, hence "1" (U for 1 dim is a vector)

    frb->weights = MALLOC(sizeof(double)*frb->maxnumofrules);
    if (frb->weights == NULL) {
#ifdef DEBUG
        printf("FIVEInit: Error - malloc(weights)\n");
#endif
        return NULL;
    }

// not required now
//    memset(frb->weights, 0, (frb->numofrules + 1) * sizeof(fri_float));

    return frb;

}
