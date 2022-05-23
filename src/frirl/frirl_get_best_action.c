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

#include "config.h"
#include "frirl.h"
#include <stdlib.h>
#include <math.h>

#include "../inl/fast_sqrt.inl"
#include "../inl/fast_abs.inl"
#include "../inl/arr.inl"
#include "../inl/min.inl"
#include "../inl/max.inl"
#include "../inl/eucledian.inl"

//#define DEBUG
//#include <stdio.h>

///get best action
///\param frirl initialized FRIRL struct
///\param states
///\return the selected action
unsigned int frirl_get_best_action(struct frirl_desc *frirl, fri_float *states) {
// frirl_get_best_action: return the proposed best action for the current state
    // cartpole calls this 33060 times

    unsigned int maxacti;
    struct FIVERB *frb = frirl->fiverb;
//    int i = 0;
    int j = 0;
    int actno;

    // calculate distances from rules for the 'state' antecedents excluding 'action' (hence -1)
    int numofstateunivs = frb->numofunivs - 1;

    fri_float *vagdist_states = frirl->fgba_vagdist_states;
    fri_float *ruledist = frirl->fgba_ruledist;
    fri_float *actconc  = frirl->fgba_actconc;
    fri_float *dists    = frirl->fgba_dists;
    fri_float *diststmp;
    fri_float *statedistsum = frirl->fgba_statedistsum; // for every rule - ruledist based on states without sqrt()

    fri_float Da = 0.0;

    unsigned int knu = 0;


    DEBUG_MSG("frirl_get_best_action: numofrules: %d rulelength: %d univlength: %d numofstateunivs: %d\n", frb->numofrules, frb->rulelength, frb->univlength, numofstateunivs);

#ifdef BUILD_AVX2

    unsigned int minjs[frb->numofunivs - 1];
    for (unsigned int k = 0; k < (frb->numofunivs - 1); k++) {
        minjs[k] = get_vag_abs_min_i_fixres(frb->uk[k], frb->univlength, states[k], frb->udivs[k]);
    }

    for (unsigned int k = 0; k < (frb->numofunivs-1); k++) {

        diststmp = dists + k*4;

        asm volatile (
                    "mov      %3, %%cx                 \n" // cx = arrsize
                    "mov      %2, %%r9                 \n" // &frb->rseqant_veval[k][0]
                    "mov      %1, %%r10                \n" // &diststmp[0]
                    "vbroadcastsd %0, %%ymm1           \n" // ymm1 = < vek[k][minj] vek[k][minj] vek[k][minj] vek[k][minj] >

             // no need for abs, will be squared anyway
            "frdl:   vmovdqa  (%%r9), %%ymm0           \n" // ymm0 = < frb->rseqant_veval[k][i] i+1 i+2 i+3 >
                    "vsubpd   %%ymm0, %%ymm1, %%ymm0   \n" // ymm0 = ymm1 - ymm0
                    "vmulpd   %%ymm0, %%ymm0, %%ymm0   \n" // ymm0 = ymm0 * ymm0 // these will be summed later
                    "vmovdqa  %%ymm0, (%%r10)          \n" // diststmp[0] = ymm0

                    "add     $32, %%r9         \n"
                    "add    $256, %%r10        \n" // 256 - 4 * FIVE_MAX_NUM_OF_UNIVERSES * sizeof(double)
//                    "dec      %%cx             \n" // sub $1 faster??
                    "sub      $1, %%cx         \n"
                    "jne      frdl\n"

//            : : "m" (frb->vek[k][minjs[k]]), "m" (diststmp), "m" (frb->rseqant_veval[k]), "m" (arrsize) : "r9", "r10", "rcx", "ymm0", "ymm1", "memory" );
            : : "m" (frb->vek[k][minjs[k]]), "m" (diststmp), "m" (frb->rseqant_veval[k]), "m" (frb->avx2_rbsize) : "r9", "r10", "cx", "ymm0", "ymm1", "memory" );

    } // for k

    // FRIRL-specific!
    // precalc ruledist but only for the states, actions will be added later
    // sum of squares for the state dim distances for each rule
    // there should be no INF values

    asm volatile (
                "mov %0,    %%rax;" // &statedistsum[0]
                "mov %1,    %%rbx;" // &dists[0]
                "mov %2,    %%cx;" // arrsize

        "fgba_l: ;"
// TODO: spekulativ vegrehajtas? tobb kulonbozo regiszterbe 1-1 kiteritett iteraciot?
                // FIVE_MAX_NUM_OF_UNIVERSES = 8 - unrolled loop for 4 rules at once - rule n...n+3
                "vmovapd (%%rbx), %%ymm3;"       // ymm3 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//                "vmovapd (%%rbx), %%ymm1;"       // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//                "vmulpd %%ymm1, %%ymm1, %%ymm3;" // ymm3 = ymm1 * ymm1

//                "vmovapd 32(%%rbx), %%ymm1;"     // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
                "vmovapd 32(%%rbx), %%ymm2;"     // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//                "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
                "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2

//                "vmovapd 64(%%rbx), %%ymm1;"     // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
                "vmovapd 64(%%rbx), %%ymm2;"     // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//                "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
                "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2

//                "vmovapd 96(%%rbx), %%ymm1;"     // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
                "vmovapd 96(%%rbx), %%ymm2;"     // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//                "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
                "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2

//                "vmovapd 128(%%rbx), %%ymm1;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
                "vmovapd 128(%%rbx), %%ymm2;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//                "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
                "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2

//                "vmovapd 160(%%rbx), %%ymm1;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
                "vmovapd 160(%%rbx), %%ymm2;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//                "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
                "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2

//                "vmovapd 192(%%rbx), %%ymm1;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
                "vmovapd 192(%%rbx), %%ymm2;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//                "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
                "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2

//                "vmovapd 224(%%rbx), %%ymm1;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
                "vmovapd 224(%%rbx), %%ymm2;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//                "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
                "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2
                // unroll end

                "vmovapd %%ymm3, (%%rax);" // frb->ruledists[i, ..., i+3] = ymm3

                "add  $32,  %%rax;" // statedistsum* += 4
                "add $256,  %%rbx;" // dists* += 32
                "sub   $1,  %%cx;"  // i-- , faster than dec cx?
                "jne fgba_l;"       // if arrsize != i, then goto fgba_l

        : //empty output, frb->ruledists in memory modification
        : "m" (statedistsum), "m" (dists), "m" (frb->avx2_rbsize)
        : "rax", "rbx", "cx", "ymm2", "ymm3", "memory"
    );


#else // C code begins here

    // ignore action universe - only the states now
    frb->numofunivs--;  // because of five_vague_distance

    for (unsigned int i = 0; i < frb->numofrules; i++) {
#ifdef DEBUG
        printf("frirl_get_best_action: vagdist inputs: ruleno(i): %d - rant: ", i);
        for (unsigned int debugi = 0; debugi < frb->rulelength - 2; debugi++) {
            printf(" %f ", frb->rant[i * (frb->rulelength - 1) + debugi]);
        }
        printf("\n");
        printf("frirl_get_best_action: vagdist inputs: ruleno(i): %d - states: ", i);
        for (unsigned int debugi = 0; debugi < frb->rulelength - 2; debugi++) {
            printf(" %f ", states[debugi]);
        }
        printf("\n");
#endif
        five_vague_distance(frirl->fiverb, frb->rant + i * (frb->rulelength - 1), states, vagdist_states + i * (frb->rulelength - 2));
    }

    frb->numofunivs++;
    //reconsider action universe

#ifdef DEBUG
    printf("frirl_get_best_action: vagdist output: vagdist_states: \n");
    for (int debugi = 0; debugi < frb->numofrules; debugi++) {
        for (int debugj = 0; debugj < (frb->rulelength - 2); debugj++) {
            printf("%f ", vagdist_states[debugi*(frb->rulelength - 2) + debugj]);
        }
        printf("\n");
    }
#endif

    // precalc ruledist but only for the states, actions will be added later
    // there should be no INF values

    memset(statedistsum, 0, sizeof(fri_float)*frb->numofrules);
    unsigned int dimstart = 0;

    for (unsigned int curruleno = 0; curruleno < frb->numofrules; curruleno++) {
        for (unsigned int dmi = 0; dmi < (frb->rulelength - 2); dmi++) { // state distances
            statedistsum[curruleno] += vagdist_states[dimstart + dmi] * vagdist_states[dimstart + dmi];
        }
        dimstart += (frb->rulelength - 2);
    }

#endif // BUILD_AVX2 - C code ends here

    // now calculate distance for possible actions:
    // for 'action' values no alignments are needed, because they must exactly hit in the action dimension

    for (actno = 0; actno < frirl->actiondim.values_len; actno++) {
        // This is taken from FIVEVagDist_FixRes !
        // cartpole 694260x hivja ezt a reszt (33060*21...)

//        double P2 = frirl->actiondim.values[actno];
// P2 is not used anymore

#ifndef FRIRL_FAST
        // as FRIRL uses Ruspini partitions NAN values are not possible - but maybe in a later version
        if (isnan(frirl->actiondim.values[actno])) { // Not a valid distance (D=0)
            Da = 0.0; // Indifferent (non existing) rule antecendent
        } else { // Valid scaled distance
#endif

#ifndef FRIRL_FAST
//TODO: was this already checked in frirl_init() ?
            // check boundaries - if ((P2 < ua[0]) || (P2 > ua[uksize])) {
            if ((frirl->actiondim.values[actno] < frirl->fiverb_ua[0]) || (frirl->actiondim.values[actno] > frirl->fiverb_ua[frb->uksize])) {
                DEBUG_MSG("frirl_get_best_action: FATAL: VagDist: frirl->actiondim.values[actno] is out of range!\n");
                DEBUG_MSG("frirl_get_best_action: frirl->actiondim.values[%d]: %lf first U: %lf last U: %lf\n", actno, frirl->actiondim.values[actno], frirl->fiverb_ua[0], frirl->fiverb_ua[frb->uksize]);

                exit(20);
            }
#endif

//            j = get_vag_abs_min_i_fixres(ua, uksize, P2, ukdomsizprod);
//            j = get_vag_abs_min_i_fixres(frirl->fiverb_ua, frb->uksize, P2, frb->udivs[frb->numofunivs-1]);
//debug volt, majd remove
//printf("act: %d %d %f == %f vallen: %d vea[j]: %lf vevalues[actno]: %lf\n", actno, j, P2, frirl->possible_actions->values[actno], frirl->actiondim.values_len, frirl->fiverb_vea[j], frirl->possible_actions->vevalues[actno]);

#ifndef FRIRL_FAST
        } // if isnan(P2)
#endif

        // check actions (hence -2: state1 state2 ... stateX action Q)

#ifdef BUILD_AVX2

        // for reference: these 2 lines were avx-ized
        //    Da = frirl->possible_actions->vevalues[actno] - frb->ract_veval[curruleno];
        //    ruledist[curruleno]=fast_sqrt(statedistsum[curruleno]+Da*Da);

        asm volatile (
                        "mov      %4, %%cx                 \n" // i = cx = frb->avx2_rbsize
                        "mov      %1, %%r8                 \n" // &frb->ract_veval[0]
                        "mov      %2, %%r9                 \n" // &statedistsum[0]
                        "mov      %3, %%r10                \n" // &ruledist[0]
                        "vbroadcastsd %0, %%ymm1           \n" // ymm1 = < vevalues[actno] vevalues[actno] vevalues[actno] vevalues[actno] >

             // no need for abs, will be squared anyway
            "fgba_l2:    vmovdqa  (%%r8), %%ymm0           \n" // ymm0 = < frb->ract_veval[i] i+1 i+2 i+3 >
                        "vsubpd   %%ymm0, %%ymm1, %%ymm0   \n" // ymm0 = ymm1 - ymm0
                        "vmulpd   %%ymm0, %%ymm0, %%ymm0   \n" // ymm0 = ymm0 * ymm0
                        "vmovdqa  (%%r9), %%ymm2           \n" // ymm2 = < statedistsum[i] i+1 i+2 i+3 >
                        "vaddpd   %%ymm2, %%ymm0, %%ymm0   \n" // ymm0 = ymm0 + ymm2
                        "vsqrtpd  %%ymm0, %%ymm0           \n" // ymm0 = sqrt(ymm0)
                        "vmovdqa  %%ymm0, (%%r10)          \n" // ruledist[i] = ymm0

                        "add     $32, %%r8         \n"
                        "add     $32, %%r9         \n"
                        "add     $32, %%r10        \n"
//                      "dec      %%cx         \n" // sub $1 faster??
                        "sub      $1, %%cx         \n"
                        "jne      fgba_l2          \n"

        : : "m" (frirl->possible_actions->vevalues[actno]), "m" (frb->ract_veval), "m" (statedistsum), "m" (ruledist), "m" (frb->avx2_rbsize) : "r8", "r9", "r10", "rcx", "ymm0", "ymm1", "ymm2", "memory" );

#else // orig C code

        for (unsigned int curruleno = 0; curruleno < frb->numofrules; curruleno++) {
            // cartpole comes here 104 779 563x times
            // acrobot  comes here  13 411 548x times

            Da = 0.0;

#ifndef FRIRL_FAST
            // as FRIRL uses Ruspini partitions NAN values are not possible
            if (isnan(frb->ract[curruleno])) {  // Invalid distance (D=0)
                Da = 0.0;  // Indifferent (non existing) rule antecendent
            } else {  // Valid scaled distance
#endif

#ifndef FRIRL_FAST
// TODO, ez kiveheto? Az elerjen az init-ben nem eleg 1x ellenorizni?
                // Safety check - should not happen...
                if ((frb->ract[curruleno] < frirl->fiverb_ua[0]) || (frb->ract[curruleno] > frirl->fiverb_ua[frb->uksize])) {
                    DEBUG_MSG("frirl_get_best_action: FATAL: VagDist: frb->ract[curruleno] is out of range!\n");
                    DEBUG_MSG("frirl_get_best_action: frb->ract[%d]: %lf first U: %lf last U: %lf\n", curruleno, frb->ract[curruleno], frirl->fiverb_ua[0], frirl->fiverb_ua[frb->uksize]);

                    exit(21);
                }
#endif
                // distance in VE between two points - proposed action vs. action value in current rule
                // no abs needed, because it will be squared anyway
                Da = frirl->possible_actions->vevalues[actno] - frb->ract_veval[curruleno];

#ifndef FRIRL_FAST
            } // if isnan()
#endif

            // There should be no INF values
            ruledist[curruleno] = fast_sqrt(statedistsum[curruleno] + Da * Da);

        }

#endif // BUILD_AVX2

#ifdef DEBUG
        printf("frirl_get_best_action: ruledist= actno: %d : ", actno);
        for (unsigned int i = 0; i < frb->numofrules; i++) {
            printf("%f ", ruledist[i]);
        }
        printf("\n");
#endif

        actconc[actno] = FIVEVagConcl_FRIRL_BestAct(frirl->fiverb, ruledist);
    }  // for actions

#ifdef DEBUG
    printf("frirl_get_best_action: actno: %d actconc: ", actno);
    for (unsigned int i = 0; i < frirl->actiondim.values_len; i++) {
        printf("%f ", actconc[i]);
    }
    printf("\n");
#endif

    maxacti = get_arr_max_i(actconc, frirl->actiondim.values_len);

    DEBUG_MSG("frirl_get_best_action: max act i: %d val: %f\n", maxacti, actconc[maxacti]);

    return maxacti;
}
