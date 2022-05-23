/*
** Fuzzy Inference by Interpolation in Vague Environment toolbox for ANSI C
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** ANSI C / AVX port and FixRes/NoInf/NoNan/etc. optimization by David Vincze <david.vincze@webcode.hu>
** Various contributions by Daniel Palko <palko.daniel@uni-miskolc.hu>
** mex port by Sandor Zsuga (zsuga@iit.uni-miskolc.hu)
** Original FIVE MATLAB code by Szilveszter Kovacs <szkovacs@iit.uni-miskolc.hu>
*/

///
///\file five_remove_rule.c
///

#define DEBUG

#include "FIVE.h"

#include "../inl/fast_abs.inl"
#include "../inl/arr.inl"
#include "../inl/min.inl"


///remove rule from the rule-base
///\param frb initialized FIVE struct
///\param rulenotoremove - rule number which should be removed
///return 0 if success
int five_remove_rule(struct FIVERB *frb, unsigned int rulenotoremove) {

#ifdef DEBUG
    printf("FIVERemoveRule: rule to remove: %d/%d. =", rulenotoremove, frb->numofrules-1);
    for (unsigned int i = 0; i < (frb->rulelength-1); i++) {
        printf(" %f", frb->rseqant[i][rulenotoremove]);
    }
    printf(" %f\n",frb->rconc[rulenotoremove]);
#endif

    if (rulenotoremove > (frb->numofrules-1)) {
        printf("FIVERemoveRule - FATAL: Invalid rule: %d max: %d !\n", rulenotoremove, frb->numofrules-1);
        exit(6);
    }

    // if the last rule is to be removed then just decrease counters
    if (rulenotoremove != frb->numofrules) {
//rb removed
//        DEBUG_MSG("FIVERemoveRule - Copying rb: %x <- %x %d\n", &frb->rb[rulenotoremove*frb->rulelength], &frb->rb[(rulenotoremove+1)*frb->rulelength], (frb->numofrules - rulenotoremove)*frb->rulelength*sizeof(fri_float) );
//        memmove(&frb->rb[rulenotoremove*frb->rulelength], &frb->rb[(rulenotoremove+1)*frb->rulelength], (frb->numofrules - rulenotoremove)*frb->rulelength*sizeof(fri_float) );

        DEBUG_MSG("FIVERemoveRule - Copying rconc: %x(%f) <- %x(%f) %d\n", &frb->rconc[rulenotoremove], frb->rconc[rulenotoremove], &frb->rconc[rulenotoremove+1], frb->rconc[rulenotoremove+1], (frb->numofrules - rulenotoremove)*sizeof(fri_float) );
        memmove(&frb->rconc[rulenotoremove], &frb->rconc[rulenotoremove+1], (frb->numofrules - rulenotoremove)*sizeof(fri_float) );

        DEBUG_MSG("FIVERemoveRule - Copying rant: %x <- %x %d\n", &frb->rant[rulenotoremove * (frb->rulelength-1)], &frb->rant[(rulenotoremove+1) * (frb->rulelength-1)], (frb->numofrules - rulenotoremove) * (frb->rulelength - 1) * sizeof(fri_float) );
        memmove(&frb->rant[rulenotoremove * (frb->rulelength-1)], &frb->rant[(rulenotoremove+1) * (frb->rulelength-1)], (frb->numofrules - rulenotoremove) * (frb->rulelength - 1) * sizeof(fri_float) );

        DEBUG_MSG("FIVERemoveRule - Copying rant_uindex: %x <- %x %d\n", &frb->rant_uindex[rulenotoremove * (frb->rulelength-1)], &frb->rant_uindex[(rulenotoremove+1) * (frb->rulelength-1)], (frb->numofrules - rulenotoremove) * (frb->rulelength - 1) * sizeof(fri_float) );
        memmove(&frb->rant_uindex[rulenotoremove * (frb->rulelength-1)], &frb->rant_uindex[(rulenotoremove+1) * (frb->rulelength-1)], (frb->numofrules - rulenotoremove) * (frb->rulelength - 1) * sizeof(unsigned int) );

        DEBUG_MSG("FIVERemoveRule - Copying rant_veval: %x <- %x %d\n", &frb->rant_veval[rulenotoremove * (frb->rulelength-1)], &frb->rant_veval[(rulenotoremove+1) * (frb->rulelength-1)], (frb->numofrules - rulenotoremove) * (frb->rulelength - 1) * sizeof(fri_float) );
        memmove(&frb->rant_veval[rulenotoremove * (frb->rulelength-1)], &frb->rant_veval[(rulenotoremove+1) * (frb->rulelength-1)], (frb->numofrules - rulenotoremove) * (frb->rulelength - 1) * sizeof(fri_float) );

        for (unsigned char currant = 0; currant < frb->rulelength - 1; currant++) {
        DEBUG_MSG("FIVERemoveRule - Copying rseqant[%d]: %x <- %x %d\n", currant, &frb->rseqant[rulenotoremove], &frb->rseqant[rulenotoremove+1], (frb->numofrules - rulenotoremove) * sizeof(fri_float) );
        memmove(&frb->rseqant[currant][rulenotoremove], &frb->rseqant[currant][rulenotoremove+1], (frb->numofrules - rulenotoremove) * sizeof(fri_float) );

        DEBUG_MSG("FIVERemoveRule - Copying rseqant_uindex[%d]: %x <- %x %d\n", currant, &frb->rseqant_uindex[rulenotoremove], &frb->rseqant_uindex[rulenotoremove+1], (frb->numofrules - rulenotoremove) * sizeof(fri_float) );
        memmove(&frb->rseqant_uindex[currant][rulenotoremove], &frb->rseqant_uindex[currant][rulenotoremove+1], (frb->numofrules - rulenotoremove) * sizeof(unsigned int) );

        DEBUG_MSG("FIVERemoveRule - Copying rseqant_uindex[%d]: %x <- %x %d\n", currant, &frb->rseqant_veval[rulenotoremove], &frb->rseqant_veval[rulenotoremove+1], (frb->numofrules - rulenotoremove) * sizeof(fri_float) );
        memmove(&frb->rseqant_veval[currant][rulenotoremove], &frb->rseqant_veval[currant][rulenotoremove+1], (frb->numofrules - rulenotoremove) * sizeof(fri_float) );
        }
    }

    frb->numofrules--;

//    frb->newrule -= frb->rulelength;
    frb->newrconc--;
    frb->newrant -= (frb->rulelength-1);

#ifdef BUILD_AVX2
    frb->avx2_rbsize = frb->numofrules / 4 + (frb->numofrules % 4 > 0);
#endif

    return 0;
}
