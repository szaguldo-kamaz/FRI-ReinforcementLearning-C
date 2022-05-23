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
**FIVEVagConcl: Calculate the conclusion from the observation and the rulebase
**
**                            [Y]=FIVEVagConcl(U,VE,R,X,P)
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
**          y: is the conclusion
**
**          Rulebase:
**
**            R1: a1 a2 a3 ... am -> b       [a1 a2 a3 ... am b;...
**            R2: a1 a2 a3 ... am -> b        a1 a2 a3 ... am b;...
**            R3: a1 a2 a3 ... am -> b        a1 a2 a3 ... am b;...
**             ...                             ...
**            Rn: a1 a2 a3 ... am -> b        a1 a2 a3 ... am b]
*/

//#define DEBUG
#include <stdio.h>

#include "FIVE.h"

#include <math.h>
#include <stdlib.h>

#ifdef DEBUG
#include <stdio.h>
#endif

#include "../inl/fast_abs.inl"
#include "../inl/fast_pow.inl"


double FIVEVagConcl_FRIRL_BestAct(struct FIVERB *frb, double *ruledists) {

    unsigned int i;

    double vagcz;
#ifdef BUILD_AVX2
    double __attribute__((aligned(32))) vagc = 0;  // Same as vagc in MatLab code
    double __attribute__((aligned(32))) ws = 0;    // Same as ws and wi in MatLab code
#else
    double vagc = 0; // Same as vagc in MatLab code
    double ws = 0;   // Same as ws and wi in MatLab code
#endif
    double wi;
    double wsz;

#ifndef FIVE_NOINF
    int infc = 0; // Count of 'inf's
#endif
    double nc;
    double y, *yp; // the return value


    /* Go through all rules; R(i,m) is the i. conclusion;
     **    ruledists(i) is the i. distance of the i. antecendent to the observation */

    /* Note: I use a different approach for working with inf distances than how
     **       it is done in the original MatLab code. Instead of filling an array
     **       and analyzing it later, i switch to infinite summing when hitting
     **       an 'inf' - Zsuga */

#ifdef FIVE_NONAN
    // with FRIRL - NANs are not allowed in the rule-base, so only one rule can match, hence if the first hit found, we're done
//TODO -1 check? atadni azt is ide frd_ret ? es akkor nem kell forciklus megegyszer
    for (i = 0; i < frb->numofrules; i++) {
        if (ruledists[i] == 0.0) {
            return frb->rconc[i];  // exact hit ?
        }
    }
#else
    for (i = 0; i < frb->numofrules; i++) {
        if (ruledists[i] == 0.0) {
            break; // exact hit ?
        }
    }

    if (i < frb->numofrules) {  // if sum(ruledists==0)>0 - at least one rule antecedent had exactly hit

#ifdef DEBUG
        printf("VagConcl: exact hit: i: %d\n", i);
#endif
        nc = 0; // num of consequents gathered
        //  FRIRL does not have VE in the consequent universe so this is never reached - The rule consequences are singletons; do finite summing only
        for (i = 0; i < frb->numofrules; i++) { // for all the rules

            if (ruledists[i] == 0.0) { // rule distance zero -> hit
                // cartpole ide 90603x jon ha nincs FIVE_NONAN
                // vrconc=[vrconc,R(i,m)]; // conclusion value
                vagc += frb->rconc[i]; // value of rule consequent
                nc++;
#ifdef DEBUG
//                printf("VagConcl: singleton consequence: i: %d, frb->rb[i*frb->rulelength+(frb->rulelength-1)]: %f vagc: %f\n", i, frb->rb[i * frb->rulelength + (frb->rulelength - 1)], vagc);
                printf("VagConcl: singleton consequence: i: %d, frb->rconc[i]: %f vagc: %f\n", i, frb->rconc[i], vagc);
#endif
            } // if ruledists[i]==0.0
        } // for all rules

        vagc /= nc; // Let the conclusion be the average of the conclusions exactly hit

    } else {  // Shepard interpolation
#endif  // FIVE_NONAN

// cartpole 603657x

#ifdef DEBUG
        printf("VagConcl: no exact hit!\n");
#endif

        wsz = 0;
        vagcz = 0; // Sum zero elements seperately to include them later in inf calculation

#ifndef FIVE_NOINF
        // search for not inf hits
        for (i = 0; i < frb->numofrules; i++) {
            if (ruledists[i] > 0.0) {
                break;
            }
        }

        if (i < frb->numofrules) { // At least one rule antecendent is not in inf distance
        // cartpole 603657x jon ide is - tehat mindig
            for (i = 0; i < frb->numofrules; i++) {
                if (ruledists[i] < 0.0) {
                    ruledists[i] = INFINITY;
                }
            }
        }
#endif

#if defined(BUILD_AVX2) && defined(BUILD_POW64)
//TODO: test this to generate new and different rulebase for a given problem, then check whether the RB makes the system work
        double absmask = -0.0;
        double one = 1.0;

// sequence of the additions are not the same as in C -> results can differ! (FP is like that...)
        asm volatile (
                        "mov      %3, %%cx                 \n" // cx = avx2_rbsize | arrsize
                        "mov      %2, %%bx                 \n" // bx = frb->p - 1
                        "dec      %%bx                     \n"
                        "mov      %1, %%r8                 \n" // &ruledists[0]
                        "mov      %0, %%r9                 \n" // &rconc[0]

                        "vbroadcastsd %5, %%ymm1           \n" // ymm1 = < 1.0 1.0 1.0 1.0 >
                        "vpxor %%ymm3, %%ymm3, %%ymm3      \n" // ymm3 = vagc = < 0.0 0.0 0.0 0.0 >
                        "vpxor %%ymm4, %%ymm4, %%ymm4      \n" // ymm4 =   ws = < 0.0 0.0 0.0 0.0 >
                        "vbroadcastsd %4, %%ymm5           \n" // ymm5 = < -0.0 -0.0 -0.0 -0.0 >

              "fvcfba:   mov        %%bx, %%dx             \n" // dx = sumi = bx = frb->p - 1
                        "vmovdqa  (%%r8), %%ymm0           \n" // ymm0 = < ruledist[i] i+1 i+2 i+3 >
                        "vandnpd  %%ymm0, %%ymm5, %%ymm0   \n" // ymm0 = ymm0 & absmask - we got abs()
                        "vmovdqa  (%%r9), %%ymm2           \n" // ymm2 = < rconc[i] i+1 i+2 i+3 >
                        "vmovdqa  %%ymm0, %%ymm6           \n" // ymm6 = ymm0
// 64-bit precision pow()-ing! not the 'default' 80-bit -> results can differ!
              "fvcfba2:  vmulpd   %%ymm0, %%ymm6, %%ymm6   \n" // ymm6 = ymm6 * ymm0
                        "sub          $1, %%dx             \n" // sumi-- (p--)
                        "jne      fvcfba2                  \n" // ymm0 = ymm0 ^ p

                        "vdivpd   %%ymm6, %%ymm1, %%ymm0   \n" // ymm0 = ymm1 / ymm6 (wi=1.0/ymm0)
                        "vaddpd   %%ymm0, %%ymm4, %%ymm4   \n" // ymm4 = ymm4 + ymm0 (ws=ws+wi)
                        "vmulpd   %%ymm0, %%ymm2, %%ymm0   \n" // ymm0 = ymm0 * ymm2 (wi * rconc)
                        "vaddpd   %%ymm0, %%ymm3, %%ymm3   \n" // ymm3 = ymm3 + ymm0 (vagc=vagc+wi*rconc)
// values in rconc[] and ruledist[] were set to 0.0 on init, so over is not a problem

                        "add     $32, %%r8         \n"
                        "add     $32, %%r9         \n"
//                      "dec      %%cx         \n" // sub $1 faster??
                        "sub      $1, %%cx         \n"
                        "jne      fvcfba           \n"

                        "vhaddpd %%ymm3, %%ymm3, %%ymm1       \n" // ymm1 = < ymm3[0] + ymm3[1], ymm3[0] + ymm3[1], ymm3[2] + ymm3[3], ymm3[2] + ymm3[3] >
                        "vextractf128    $0x1, %%ymm3, %%xmm0 \n" // xmm0 = < ymm3[2] + ymm3[3], ymm3[2] + ymm3[3] >
                        "vaddsd  %%xmm0, %%xmm1, %%xmm3       \n" // xmm3 = < vagc 0.0 > (vagc=ymm3[0] + ymm3[1] + ymm3[2] + ymm3[3])
                        "vmovdqa %%xmm3, %6   \n"
//                        "vmovdqu %%xmm3, %6   \n"

                        "vhaddpd %%ymm4, %%ymm4, %%ymm1       \n" // ymm1 = < ymm4[0] + ymm4[1], ymm4[0] + ymm4[1], ymm4[2] + ymm4[3], ymm4[2] + ymm4[3] >
                        "vextractf128    $0x1, %%ymm3, %%xmm0 \n" // xmm0 = < ymm4[2] + ymm4[3], ymm4[2] + ymm4[3] >
                        "vaddsd  %%xmm0, %%xmm1, %%xmm4       \n" // xmm4 = < ws 0.0 > (ws=ymm4[0] + ymm4[1] + ymm4[2] + ymm4[3])
                        "vmovdqa %%xmm4, %7   \n"
//                        "vmovdqu %%xmm4, %7   \n"

        : : "m" (frb->rconc), "m" (ruledists), "m" (frb->p), "m" (frb->avx2_rbsize), "m" (absmask), "m" (one), "m" (vagc), "m" (ws) :
        "r8", "r9", "cx", "dx", "ymm0", "ymm1", "ymm2", "ymm3", "ymm4", "ymm5", "ymm6", "memory" );

#else  // C code

        // FRIRL - The rule consequences are always singletons
        for (i = 0; i < frb->numofrules; i++) {  // for all the rules
            // cartpole 90948003x
            // acrobot   3520876x
            wi = 1.0 / fast_pow(fast_abs(ruledists[i]), frb->p);
            vagc += wi * frb->rconc[i]; // value of rule consequent
            ws += wi;

/* NOTES
// TODO - szekvencialisan osszeadni, skalaris avx
wi[0] = fast_pow(fast_abs(ruledists[i]), frb->p);
wi[1] = fast_pow(fast_abs(ruledists[i+1]), frb->p);
wi[2] = fast_pow(fast_abs(ruledists[i+2]), frb->p);
wi[3] = fast_pow(fast_abs(ruledists[i+3]), frb->p);
ymm0=wi
ymm1=frb->rconc[i].. [i+1]...[i+3]
ymm0=ymm0*ymm1
extracf128 , 1 -> xmm2
punpckhqdq xmm0, xmm0    # broadcast the high half of xmm0 to both halves
movhlps
unpcklpd

    movshdup    xmm1, xmm0
    addps       xmm0, xmm1
    movhlps     xmm1, xmm0
    addss       xmm0, xmm1

    movhlps xmm1, xmm0
    addsd   xmm0, xmm1

*/
#ifdef DEBUG
            //todo rb removed
//            if (frb->rb[i * frb->rulelength + (frb->rulelength - 1)] != frb->rconc[i]) {
//                printf("This should not happen: %a != %a\n", frb->rb[i * frb->rulelength + (frb->rulelength - 1)], frb->rconc[i]);
//                printf("FATAL!: i: %d frb->rconc[i]: %f i * frb->rulelength + (frb->rulelength - 1): %d frb->rb[i * frb->rulelength + (frb->rulelength - 1)]: %f\n",i,frb->rconc[i], i * frb->rulelength + (frb->rulelength - 1), frb->rb[i * frb->rulelength + (frb->rulelength - 1)]);
//                exit(50);
//            }
//
//            printf("VagConcl: singleton consequence: i: %d, wi: %f, ruledists[i]: %f p: %d fabs(ruledists[i]): %f pow(fabs(ruledists[i]),p): %f\n", i, wi, ruledists[i], frb->p, fabs(ruledists[i]), pow(fabs(ruledists[i]), frb->p));
//            printf("VagConcl: singleton consequence2: i: %d, wi: %10.25f, frb->rconc[i]: %f vagc: %f rulelength: %d\n", i, wi, frb->rconc[i], vagc, frb->rulelength);
#endif
        }

#endif // BUILD_AVX2 && BUILD_POW64

#ifndef FIVE_NOINF
        if (ws == 0.0) {  // There is no valid conclusion (infinite distances from all the rules)
            vagc = NAN;
    #ifdef DEBUG
            printf("VagConcl: no valid conclusion: vagc: %f\n", vagc);
    #endif
        } else {
#endif
            vagc /= ws;  // The valid interpolated conclusion
#ifdef DEBUG
            printf("VagConcl: valid interpolated conclusion: vagc/ws: %f ws: %f\n", vagc, ws);
#endif
#ifndef FIVE_NOINF
        }
#endif

#ifndef FIVE_NONAN
    }
#endif

#ifndef FIVE_NONAN
    if (isnan(vagc)) {
        y = NAN;
    #ifdef DEBUG
        printf("VagConcl: y=NAN\n");
    #endif
    } else {
    // with FRIRL, consequent has no VE - The rule consequences are always singletons
        y = vagc;
    #ifdef DEBUG
        printf("VagConcl: consequence singleton: y=vagc: %f\n", y);
    #endif
    }

    return y;

#else
    #ifdef DEBUG
    printf("VagConcl: consequence singleton: vagc: %f\n", vagc);
    #endif

    return vagc;
#endif

}
