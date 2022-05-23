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
///\file five_add_rule.c
///

//#define DEBUG

#include "FIVE.h"

#include "../inl/fast_abs.inl"
#include "../inl/arr.inl"
#include "../inl/min.inl"

///add rule into the rule base. Old API, kept for backward compatibility.
///\param frb initialized FIVE struct
///\param newrule new rule to be added
///return 0 if success
int FIVEAddRule(struct FIVERB *frb, double *newrule) {
    return FIVE_add_rule(frb, newrule, newrule[frb->numofantecedents]);
}


///add rule into the rule base. Obsolated API, to be removed.
///\param frb initialized FIVE struct
///\param ruletoadd new rule to be added
///return 0 if success
int five_add_rule(struct FIVERB *frb, fri_float *ruletoadd) {
    return FIVE_add_rule(frb, ruletoadd, ruletoadd[frb->numofantecedents]);
}


///add rule into the rule base
///\param frb initialized FIVE struct
///\param rant new rule antecedents to be added
///\param rconc new rule conclusion to be added
///return 0 if success
int FIVE_add_rule(struct FIVERB *frb, fri_float *rant, fri_float rconc) {
    // copy new rule into the rule-base (as the last rule)

#ifdef DEBUG
    printf("FIVEAddRule: new rule: %d. -", frb->numofrules);
    for (unsigned int i = 0; i < frb->rulelength; i++) {
        printf(" %f", rule[i]);
    }
    printf("\n");
#endif

    // no prob if data at the end gets overwritten, it's just trash anyway, so we can copy more with using xmm/ymm regs maybe
//    arr_cpy_dirty(ruletoadd, frb->newrule, frb->rulelength);
//    arr_cpy_dirty(rant, frb->newrant, frb->rulelength - 1); //ez miert kell?? par sorral lentebb megint ezt csinalja.
    *frb->newrconc = rconc;

    DEBUG_MSG("FIVEAddRule - Copying to rconc\n");
    DEBUG_MSG("FIVEAddRule - Q: %f -> frb->rconc[%d]\n", frb->newrule[frb->numofantecedents], frb->numofrules);
    //ez miert is kell?? nem ua mint newrconc?
    frb->rconc[frb->numofrules] =  *frb->newrconc;  //frb->newrule[frb->rulelength - 1];

    DEBUG_MSG("FIVEAddRule - Copying to rant: %d\n", (frb->numofrules) * frb->numofantecedents);
    arr_cpy_dirty(rant, frb->newrant, frb->numofantecedents);

    for (unsigned char currant = 0; currant < frb->numofantecedents; currant++) {
        DEBUG_MSG("FIVEAddRule - Copying to rseqant[%d][%d]: %lf\n", currant, frb->numofrules, frb->newrule[currant]);
        frb->rseqant[currant][frb->numofrules] = frb->newrant[currant];

        DEBUG_MSG("FIVEAddRule - Copying to rseqant_uindex[%d][%d]: %lf\n", currant, frb->numofrules, get_vag_abs_min_i_fixres(frb->u + frb->numofantecedents * currant, frb->rulelength - 1, frb->rseqant[currant][frb->numofrules], frb->udivs[currant]));
        frb->rseqant_uindex[currant][frb->numofrules] = get_vag_abs_min_i_fixres(frb->u + frb->univlength * currant, frb->univlength, frb->rseqant[currant][frb->numofrules], frb->udivs[currant]);
        DEBUG_MSG("FIVEAddRule - Copying to rant_uindex[%d]: %lf\n", frb->numofrules*frb->numofantecedents+currant, frb->rseqant_uindex[currant][frb->numofrules]);
        frb->rant_uindex[frb->numofrules*frb->numofantecedents+currant] = frb->rseqant_uindex[currant][frb->numofrules];

        DEBUG_MSG("FIVEAddRule - Copying to rseqant_veval[%d][%d]: %lf\n", currant, frb->numofrules, frb->vek[currant][ frb->rseqant_uindex[currant][frb->numofrules] ]);
        frb->rseqant_veval[currant][frb->numofrules] = frb->vek[currant][ frb->rseqant_uindex[currant][frb->numofrules] ];
        DEBUG_MSG("FIVEAddRule - Copying to rant_veval[%d]: %lf\n", frb->numofrules*frb->numofantecedents+currant, frb->rseqant_veval[currant][frb->numofrules]);
        frb->rant_veval[frb->numofrules*frb->numofantecedents+currant] = frb->rseqant_veval[currant][frb->numofrules];
    }

    frb->numofrules++;
//    frb->newrule += frb->rulelength;
    frb->newrconc++;
    frb->newrant += frb->numofantecedents;
#ifdef BUILD_AVX2
    frb->avx2_rbsize = frb->numofrules / 4 + (frb->numofrules % 4 > 0);
#endif

    return 0;
}
