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
*/

#include "FIVE.h"

#include "../inl/fast_abs.inl"
#include "../inl/arr.inl"
#include "../inl/min.inl"
#include "../inl/sort.inl"

#ifdef DEBUG
#include <stdio.h>
#endif


///
///\file five_vague_distance_parallel.c
///
/// Calculate parallel the scaled distance of 4 x two points in the vague environment
///
///    d =    r1,d1 r2,d1 r3,d1 r4,d1 r1,d2 r2,d2 r3,d2 r4,d2 ...
///
///        ... r1,dj-1, r2,dj-1, r3,dj-1, r4,dj-1, r1,dj, r2,dj, r3,dj, r4,dj
///
///    r=(1..4), d=(1..j), j=numofunivs
///
/// This vagdist version leaves out from the calculation the infinite and NaN values!! Fixres only!!
///
///\param frb five struct
///\param p1 point A array
///\param p2 point B
///\param d output array
///\return void
int five_vague_distance_parallel(struct FIVERB *frb, fri_float *p1, int p1_offset, fri_float *p2, fri_float *d) {
#ifdef BUILD_AVX2off
    // TODO, easy :-P
#else
    for (unsigned int k = 0; k < frb->numofunivs; k++) {

        unsigned int knu = k * frb->univlength;

        fri_float ukdomain = frb->u[knu + frb->univlength - 1] - frb->u[knu];  // u domain length (max(u) - min(u))
        fri_float udiv = ukdomain / (frb->univlength - 1);

        unsigned int minj = get_vag_abs_min_i_fixres(frb->u + knu, frb->univlength, *(p2 + k), udiv);

        fri_float *veknu = frb->ve + knu;

        for (unsigned int i = 0; i < 4; i++) {
            int minjj = minj;
            fri_float *p1tmp = p1 + i * p1_offset;
            unsigned int mini = get_vag_abs_min_i_fixres(frb->u + knu, frb->univlength, *(p1tmp + k), udiv);

            sort_swap(&mini, &minjj);

            d[k * 4 + i] = veknu[minjj] - veknu[mini];
        }
    }
#endif

    return 0;
}
