/*
** Fuzzy Inference by Interpolation in Vague Environment toolbox for ANSI C
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** ANSI C / AVX port and FixRes/NoInf/NoNan/etc. optimization by David Vincze <david.vincze@webcode.hu>
** Various contributions by Daniel Palko <palko.daniel@uni-miskolc.hu>
** mex port by Sandor Zsuga (zsuga@iit.uni-miskolc.hu)
** Original FIVE MATLAB code by Szilveszter Kovacs <szkovacs@iit.uni-miskolc.hu>
**
**
**FIVEVagConclWeight: Calculate the weight of the rule matches from observation
**
**                            [Y]=FIVEVagConclWeight(U,VE,R,X,P)
**
**          Where:
**          U: is the universe (a vector of discrete values in increasing order),
**          VE: is the generated vague environment (the scaled distances on the universe U)
**             If VE(k,1)=0, it is the primitive integral of the scaling function according
**               to the elements of U
**             If VE(k,1)=-1, it is the integral of the scaling function between
**               the neighboring elements of U
**          R: is the Rulebase
**          X: is the observation
**          P: is the power factor in the Shepard interpolation formula: wi=1./dist(i).^p
**             optional, if not given, by default it is equal to the
**             antecedent dimensions of the rulebase R
**          In case of U,VE,X, the rows are the dimensions
**          Y: is the pointer to the calculated weigths array (the weight for every rule)
*/

//#define DEBUG

#include "FIVE.h"

#include <math.h>
#include <stdlib.h>

#ifdef DEBUG
#include <stdio.h>
#endif

#include "../inl/fast_abs.inl"
#include "../inl/fast_pow.inl"


unsigned int FIVEVagConclWeight(struct FIVERB *frb, double *x) {
    return FIVE_vag_concl_weight(frb, x, frb->weights);
}


unsigned int FIVE_vag_concl_weight(struct FIVERB *frb, double *ant, double *weights) {
    unsigned int     i;

    //double  *vagc = frb->weights;
    double  *vagc = weights;
    double  ws = 0.0;

    // Calculate the distances of the observation from the rule antecendents
    int frd_ret = five_rule_distance(frb,ant); // result distances go into frb->ruledists
    //TODO ez nem ujrahasznosithato korabbi "sima" VagConcl hivasbol?

 /* Go through all rules; R(i,m) is the i. conclusion;
 **    rd(i) is the i. distance of the i. antecendent to the observation */

#if defined(FRIRL_FAST) && defined(FIVE_NONAN)
    // if ret is not -1 then it means that there was a direct hit, and returns with the rule no
    if (frd_ret != -1) {
        return frd_ret;
    }

    // If i<numofrules then there was a rule which exactly hit - return with the rule number
    for (i=0; i<frb->numofrules; i++) {
        if (frb->ruledists[i] == 0.0) {
            return i;
        }
    }
#else
    // If i<numofrules then there was a rule which exactly hit
    for (i=0; i<frb->numofrules; i++) {
        if (frb->ruledists[i] == 0.0) {
            break;
        }
    }

    if (i < frb->numofrules) {  // if sum(frb->ruledists==0)>0 - at least one rule antecedent had exactly hit
        // cartpole 32102x
        // acrobot   9846x
#ifdef DEBUG
        printf("VagConclWeight: exact hit: i: %d\n",i);
#endif
        for (unsigned int i = 0; i < frb->numofrules; i++) {  // for all the rules
            if (frb->ruledists[i] == 0.0) {  // rule distance zero -> hit
                vagc[i] = 1.0;
            } else {
                vagc[i] = 0.0;
            }
        }

    } else {  // Shepard interpolation

#endif // defined(FRIRL_FAST) && defined(FIVE_NONAN)

#ifdef DEBUG
        printf("VagConclWeight: no exact hit!\n");
#endif

#ifndef FIVE_NOINF
        // search for not inf hits
        for (i=0; i < frb->numofrules; i++) {
            if (frb->ruledists[i] > 0.0) {
                break;
            }
        }

        if (i < frb->numofrules) {  // At least one rule antecendent is not in inf distance
            for (i=0; i < frb->numofrules; i++) {
                if (frb->ruledists[i] < 0.0) {
                    frb->ruledists[i] = INFINITY;
                }
            }
        }
#endif

        for (i=0; i < frb->numofrules; i++) { // for all the rules
        // acrobot     1715456x
        // cartpole      83976x
        // mountaincar    2731x

            // TODO: haven't this been calculated before???
            frb->wi[i] = 1.0 / fast_pow(fast_abs(frb->ruledists[i]), frb->p);
            ws += frb->wi[i];
#ifdef DEBUG
            printf("VagConclWeight: singleton consequence: i: %d, frb->wi[i]: %10.25f, frb->ruledists[i]: %f p: %d fabs(frb->ruledists[i]): %f pow(fabs(frb->ruledists[i]),p): %f\n", i, frb->wi[i], frb->ruledists[i], frb->p, fabs(frb->ruledists[i]), pow(fabs(frb->ruledists[i]),frb->p) );
#endif
        }

#ifndef FIVE_NOINF
        if (ws==0.0) { // There is no valid conclusion (infinite distances from all the rules)
            for (i=0; i<frb->numofrules; i++) { // for all the rules
                vagc[i]=NAN;
            }
    #ifdef DEBUG
            printf("VagConclWeight: no valid conclusion: vagc: NAN\n");
    #endif
        } else {
#endif

        // cartpole   750x
        // acrobot  11026x

#ifdef BUILD_AVX2
            asm volatile (
                    "xor     %%rcx, %%rcx              \n" // i=0
                    "shl     $0x2, %3                  \n" // arrsize*=4
                    "vbroadcastsd %2, %%ymm1           \n" // ymm1 = < ws ws ws ws >

            "fvcwl:  vmovdqa (%1,%%rcx,8), %%ymm0      \n" // ymm0 = wi[i] wi[i+1] wi[i+2] wi[i+3]
                    "vdivpd  %%ymm1, %%ymm0, %%ymm0    \n" // ymm0 = ymm0 / ymm1
                    "vmovdqa %%ymm0, 0x0(%0,%%rcx,8)   \n" // vagc[i..i+3] = ymm0

                    "add     $0x4, %%ecx               \n"
                    "cmp     %%ecx, %3                 \n"
                    "jnz     fvcwl                     \n"

             : : "r" (vagc), "r" (frb->wi), "m" (ws), "r" (frb->avx2_rbsize) : "rcx", "ymm0", "ymm1", "memory");

#else

            for (i=0; i < frb->numofrules; i++) { // for all the rules
                vagc[i] = frb->wi[i] / ws;  // The valid interpolated conclusion
    #ifdef DEBUG
                printf("VagConclWeight: valid interpolated conclusions weights: vagc[i]: %f ws: %10.25f vagc[i]/ws: %f\n", vagc[i], ws, vagc[i]/ws);
    #endif
            }
#endif // BUILD_AVX2

#ifndef FIVE_NOINF
        }
#endif

#if !(defined(FRIRL_FAST) && defined(FIVE_NONAN))
    }
#endif

    return ~0;

}
