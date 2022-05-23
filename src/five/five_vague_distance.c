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
**FIVEVagDist:  Calculate the scaled distance of two points in the vague environment
**
**                            [D]=FIVEVagDist(U,VE,P1,P2)
**
**          Where:
**          U: is the universe (a vector of discrete values in increasing order),
**          VE: is the generated vague environment (the scaled distances on the universe U)
**             If VE(k,1)=0, it is the primitive integral of the scaling function according
**               to the elements of U
**             If VE(k,1)=-1, it is the integral of the scaling function between
**               the neighboring elements of U
**          P1,P2 is two points (values), inside U
**               If P1=NaN, or P2=NaN, then the vague distance is null (D=0)
**          D: is the scaled distance of the two points P1,P2.
**               D<0 denotes that D=inf. This case abs(D) is the number of inf
**               values in VE between P1 and P2
**          In case of U,VE,P,D the rows are the dimensions
*/

///
///\file five_vague_distance.c
///

#include "FIVE.h"

#include "../inl/fast_abs.inl"
#include "../inl/arr.inl"
#include "../inl/min.inl"
#include "../inl/sort.inl"
#include <stdlib.h>
#include <stdio.h>


/// Calculate the scaled distance of two points in the vague environment
///\param frb five struct
///\param p1 point A
///\param p2 point B
///\param d output array
///\return void
int five_vague_distance(struct FIVERB *frb, fri_float *p1, fri_float *p2, fri_float *d) {

    for (unsigned int k = 0; k < frb->numofunivs; k++) {
        unsigned int knu = k * frb->univlength;
        unsigned int knu2 = knu + frb->univlength - 1;

#ifdef DEBUG
        printf("k: %d p1[k]: %f knu: %d frb->u[knu]: %f p2[k]: %f knu2: %d frb->u[knu2]: %f\n",k,p1[k],knu,frb->u[knu],p2[k],knu2,frb->u[knu2]);
#endif

#ifndef FIVE_NONAN
        if (unlikely(isnan(p1[k]) || isnan(p2[k]))) {
            d[k] = 0.f; // Indifferent (non existing) rule antecendent

            DEBUG_MSG("VagDist: the points are NAN values\n");

            continue;
        }
#endif

#ifndef FRIRL_FAST
        if (unlikely(((p1[k] < frb->u[knu]) || (p2[k] < frb->u[knu])) ||
            ((p1[k] > frb->u[knu2]) || (p2[k] > frb->u[knu2])))) {

            printf("FATAL: VagDist: The points are out of range!\n");

            exit(-1);
        }
#endif

#ifdef FIVE_FIXRES
// left here for reference
//        fri_float ukdomain = frb->u[knu + frb->univlength - 1] - frb->u[knu]; // u domain length (max(u) - min(u))
//        fri_float udiv = ukdomain / (frb->univlength - 1);

        unsigned int i = get_vag_abs_min_i_fixres(frb->u + knu, frb->univlength, p1[k], frb->udivs[k]);
        unsigned int j = get_vag_abs_min_i_fixres(frb->u + knu, frb->univlength, p2[k], frb->udivs[k]);
#else
        unsigned int i = get_vag_abs_min_i(frb->u + knu, frb->univlength, p1[k]);
        unsigned int j = get_vag_abs_min_i(frb->u + knu, frb->univlength, p2[k]);
#endif

        sort_swap(&i, &j);

        fri_float *veknu = frb->ve + knu;

        if (likely(*veknu >= 0)) {

            d[k] = veknu[j] - veknu[i];

        } else {
            d[k] = get_arr_sum(veknu + i + 1, j - i - 1);

#ifndef FIVE_NOINF
            if (unlikely(d[k] == INFINITY)) {
                d[k] = -get_arr_inf_count(veknu + i + 1, j - i - 1);
            }
#endif
        }
    }

    return 0;
}
