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
**FIVEValVag:  Value of a vague point VP
**
**                            [P]=FIVEValVag(U,VE,VP)
**
**          Where:
**          U: is the universe (a vector of discrete values in increasing order),
**          VE: is the generated vague environment (the scaled distances on the universe U)
**             If VE(1)=0, it is the primitive integral of the scaling function according
**               to the elements of U
**             If VE(1)=-1, it is the integral of the scaling function between
**               the neighboring elements of U
**          VP: is the scaled distance of a point P from the first element 
**               of the vague environment 
**               VP<0 denotes infinite conclusion distance
**          P: is the point which has the VP scaled distance from the first element 
**               of the vague environment
**          In case of U,VE,VP,P the rows are the dimensions
*/

//#define DEBUG

#include "FIVE.h"

#include <math.h>
#include <stdlib.h>

#ifdef DEBUG
#include <stdio.h>
#endif

#include "../inl/fast_abs.inl"


int FIVEValVag(struct FIVERB *frb, double *vp) {

    int    i, k, l;
    double svp, vpc;
    double t;


    for (k=0; k < frb->valvagdims; k++) {  // for all the dimensions - (normally only 1 - because this is called only for the consequent dim)
#ifdef DEBUG
        printf("ValVag: k: %d univlength: %d ve[k*univlength]: %f vp[k]: %f *vp: %f\n", k, frb->univlength, frb->valvagve[k * frb->univlength], vp[k], *vp);
#endif
        if (frb->ve[k*frb->univlength] < 0) {  // check first VE element in the current dim - VE(k,1)=-1, VE is the integral of the scaling function between the neighboring elements of U
            i = 1;
#ifndef FIVE_NOINF
            if (vp[k] < 0) { // Denotes infinite conclusion distance
                vpc = vp[k];
                while ( (vpc < 0) && (i < (frb->univlength-1)) ) {
                    if (frb->valvagve[k*frb->univlength+i] == INFINITY) {  // count inf
                        vpc++;
                    }
                    i++;
                }
            } else {
#endif
                svp = frb->valvagve[k * frb->univlength + i];  // init
                while ( (svp < vp[k]) && (i < (frb->univlength-1)) ) {

        /* Note: determining if the input is out of range is faster this way
        **       than by seperately summing (like how it was in the original MatLab code) */
/*        if (l==univlength) {
            printf("ValVag: The points are out of range!");
            return -1;
        }*/
                    i++;
                    svp += frb->valvagve[k * frb->univlength + i];  // sum up the neighbouring differences in VE
                }
                // the previous was closer
                if ( (svp-vp[k]) > (vp[k] - (svp-frb->valvagve[k*frb->univlength+i])) ) {
                    i--;
                }
#ifndef FIVE_NOINF
            }
#endif
#ifdef DEBUG
            printf("ValVag: svp: %f k*univlength+i: %d ve[k*univlength+i]: %f\n", svp, k*frb->univlength+i, frb->valvagve[k*frb->univlength+i]);
#endif

        } else { // VE(k, 1) = 0, VE is the primitive integral of the scaling function according to the elements of U

#ifdef DEBUG
            printf("ValVag: k: %d univlength: %d k*univlength: %d univlength-1: %d k*univlength+(univlength-1): %d\n",k,frb->univlength,k*frb->univlength,frb->univlength-1,k*frb->univlength+(frb->univlength-1));
            printf("ValVag: vp[k]: %f ve[k*univlength+(univlength-1)]: %f\n",vp[k],frb->valvagve[k*frb->univlength+(frb->univlength-1)]);
#endif
            if ( (vp[k] < 0) || (vp[k] > frb->valvagve[k*frb->univlength+(frb->univlength-1)]) ) {
#ifdef DEBUG
                printf("The points are out of range!");
#endif
                return -1;
            }

//        [dm,i]=min(abs(VE(k,:)-VP(k)));   % i is the index of the min distance from vp (VE(k,1)=0)

            i = 0; // min index
#ifdef DEBUG
            printf("ValVag: minindex: k: %d vp[k]: %f t: %f valvagve[0]: %f\n", k, vp[k], t, frb->valvagve[0]);
#endif
            t = fast_abs(frb->valvagve[0] - vp[k]); // Minimum value
#ifdef DEBUG
            printf("ValVag: minvalue: k: %d vp[k]: %f t: %f valvagve[0]: %f\n", k, vp[k], t, frb->valvagve[0]);
#endif
            for (l=1; l < frb->univlength; l++) {
                if (t > fast_abs(frb->valvagve[k*frb->univlength+l] - vp[k])) {
                    t = fast_abs(frb->valvagve[k*frb->univlength+l] - vp[k]);
#ifdef DEBUG
                    printf("ValVag: l: %d k: %d vp[k]: %f t: %f valvagve[k*univlength+l]: %f\n", l, k, vp[k], t, frb->valvagve[k * frb->univlength + l]);
#endif
                    i = l;
                }
            }

#ifdef DEBUG
            printf("ValVag: i: %d t: %f\n", i, t);
#endif

        }

        // the point which has the VP scaled distance from the first element of the vague environment
        frb->valvagp[k] = frb->valvagu[i];

    }  // for

    return 0;

}
