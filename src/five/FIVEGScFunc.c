/*
** Fuzzy Inference by Interpolation in Vague Environment toolbox for ANSI C
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** ANSI C / AVX port and FixRes optimization by David Vincze <david.vincze@webcode.hu>
** Various contributions by Daniel Palko <palko.daniel@uni-miskolc.hu>
** mex port by Sandor Zsuga (zsuga@iit.uni-miskolc.hu)
** Original FIVE MATLAB code by Szilveszter Kovacs <szkovacs@iit.uni-miskolc.hu>

**
**
** FIVEGScFunc:  Generate the scaling function SCF
**
**                            [SCF]=FIVEGScFunc(U,mu,nu,PSC,mp,np,NLS)
**
**          Where:
**          U: is the universe (a vector of discrete values in increasing order),
**          PSC: contains the points of the scaling function
**               PSC: S
**                   where S is the scaling value of the whole universe, or
**               PSC: P1,S1;P2,S2;...
**                   where P is the point, S is the scaling function value, or
**               PSC: P1,S1l,S1r;P2,S2l,s2r;...
**                   where P is the point, Sl the left, Sr the right scaling function value.
**          NLS: optional,
**               if not given: linear scaling function approximation
**               if given: NLS is the constant factor of sensitivity for
**                   neighbouring scaling factor differences in nonlinear
**                   scaling function approximation
**          SCF: is the generated scaling function 
**               (scaling factors according to the elements of U)
**          In case of U,PSC,SCF the rows are the dimensions
**
**          In nonlinear case, the scaling function is approximated by the
**          following function (y = 1/(x^w)) between the neighbouring scaling factor
**          values (Pn,Pn+1):
**
**               SCF(x)=w.*(((d+1)./(x+1)).^w-1)./((d+1).^w-1), or
**               SCF(x)=w.*(((d+1)./(d-x+1)).^w-1)./((d+1).^w-1),
**                   where
**                       w = NLS*abs(S(n+1)-S(n)), d=P(n+1)-P(n)
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


//rc5 API
double *FIVEGScFunc(double *u, int numofunivs, int univlength, double *psc, int mp, int np, double nls) {

    int usize = numofunivs * univlength; /* size of U */
    double *scf = MALLOC(usize * sizeof(double));
    if (scf == NULL) {
        return NULL;
    }

    int ret = FIVE_GSc_func(u, numofunivs, univlength, psc, mp, np, nls, scf);
    if (ret != 0)
        return NULL;

    return scf;

}


int FIVE_GSc_func(double *u, int numofunivs, int univlength, double *psc, int mp, int np, double nls, double *scf) {

    int su = numofunivs * univlength; // size of U
    double c;  // Approximation type, scaling
    int i, j, p;
    double p1, p2, s1, s2, d, w;

    c = nls;  // Nonlinear approximation, scaling of the power function
    if ((c <= 0) || (c == NAN)) {
#ifdef DEBUG
        printf("NLS must be positive!\n");
#endif
        return -1;
    } else {
        c = NAN;
    }

    /*
    mu/numofunivs=1 - number of rows
    nu/univlength=1001 - number of columns
    l=2 mp - number of rows
    k=3 np - number of columns
    */

    if (mp > 1) { // More than one point

        for (i = 0; i < numofunivs; i++) {

            for (j = 0; j < univlength; j++) {

                if (u[i * univlength + j] < psc[0]) {
                    scf[i * univlength + j] = psc[1];
                } // First

#ifdef DEBUG
                printf("GScFunc: i: %d j: %d\n", i, j);
#endif
                for (p = 0; p < (mp - 1); p++) {
#ifdef DEBUG
                    printf("GScFunc: p  : %d p*np   : %d psc[p*np]    : %f    i*univlength+j: %d u[i*univlength+j]: %f\n", p, p*np, psc[p * np], i * univlength + j, u[i * univlength + j]);
                    printf("GScFunc: p+1: %d (p+1)*n: %d psc[(p+1)*np]: %f\n", p + 1, (p + 1) * np, psc[(p + 1) * np]);
#endif
                    if ((u[i * univlength + j] >= psc[p * np]) && (u[i * univlength + j] < psc[(p + 1) * np])) {
                        p1 = psc[p * np];
                        p2 = psc[(p + 1) * np];
#ifdef DEBUG
                        printf("GScFunc: p1: %f p2: %f\n", p1, p2);
#endif
                        if (np == 2) {
                            s1 = psc[p * np + 1]; // Sl and Sr are the same
                        } else {
                            s1 = psc[p * np + 2]; // Sl and Sr are different
                        }
                        s2 = psc[(p + 1) * np + 1];
#ifdef DEBUG
                        printf("GScFunc: s1: %f s2: %f\n", s1, s2);
#endif

                        /* The scaling function is constant */
                        if (s1 == s2) {
                            scf[i * univlength + j] = s1;
                            continue;
                        }

                        /* The scaling function is not constant */
                        if (isnan(c)) {
#ifdef DEBUG
                            printf("GScFunc: i*univlength+j: %d (s2-s1)/(p2-p1): %f u[i*univlength+j]-p1: %f s1: %f\n", i * univlength + j, (s2 - s1) / (p2 - p1), u[i * univlength + j] - p1, s1);
#endif
                            /* Linear approximation */
                            scf[i * univlength + j] = ((s2 - s1) / (p2 - p1))*(u[i * univlength + j] - p1) + s1;
#ifdef DEBUG
                            printf("GScFunc: scf[i*univlength+j]: %f\n", scf[i * univlength + j]);
#endif
                            continue;
                        }

                        /* Nonlinear approximation */
                        d = p2 - p1;
                        if (s1 > s2) {
                            w = s1 - s2;
#ifdef DEBUG
                            printf("GScFunc: d: %f w: %f u[i*univlength+j]: %f\n", d, w, u[i * univlength + j]);
#endif
                            if (w == INFINITY) {  // s1 == inf
                                if (u[i * univlength + j] > p1) {
                                    scf[i * univlength + j] = s2;
                                } else {
                                    scf[i * univlength + j] = s2 + w;
                                }
                                continue;
                            }
                            scf[i * univlength + j] = s2 + w * (fast_pow((d + 1) / (u[i * univlength + j] - p1 + 1), c * w) - 1) / (fast_pow((d + 1), c * w) - 1);

                        } else {

                            w = s2 - s1;
#ifdef DEBUG
                            printf("GScFunc: d: %f w: %f u[i*univlength+j]: %f\n", d, w, u[i * univlength + j]);
#endif
                            if (w == INFINITY) { /* s2==inf */
                                if (u[i * univlength + j] > p2) {
                                    scf[i * univlength + j] = s1;
                                } else {
                                    scf[i * univlength + j] = s1 + w;
                                }
                                continue;
                            }
                            scf[i * univlength + j] = s1 + w * (fast_pow((d + 1) / (p2 - u[i * univlength + j] + 1), c * w) - 1) / (fast_pow((d + 1), c * w) - 1);
                        }

                    } // end if
                } // end for(p ...

#ifdef DEBUG
                printf("GScFunc: psc[(mp-1)*np]:   %f (mp-1)*np:   %d\n", psc[(mp - 1) * np], (mp - 1) * np);
                printf("GScFunc: psc[(mp-1)*np+1]: %f (mp-1)*np+1: %d\n", psc[(mp - 1) * np], (mp - 1) * np + 1);
#endif
                if (u[i * univlength + j] >= psc[(mp - 1) * np]) { // Last

                    if (np == 2) { /* Sl and Sr are the same */
                        scf[i * univlength + j] = psc[(mp - 1) * np + 1];
                        continue;
                    }

                    /* Sl and Sr are different */
#ifdef DEBUG
                    printf("GScFunc: j: %d univlength: %d i*univlength+j: %d (mp-1)*np]: %d psc[(mp-1)*np]:   A%f\n", j, univlength, i * univlength + j, (mp - 1) * np, psc[(mp - 1) * np]);
#endif
                    if ((j == univlength - 1) && (u[i * univlength + j] == psc[(mp - 1) * np])) { /* The last scaling point is at the end */
                        scf[i * univlength + j] = psc[(mp - 1) * np + 1];
                    }/* The Sl left scaling value is applied */
                    else {
#ifdef DEBUG
                        printf("GScFunc: j: %d univlength: %d i*univlength+j: %d (mp-1)*np+2: %d psc[(mp-1)*np+2]: %f\n", j, univlength, i * univlength + j, (mp - 1) * np + 2, psc[(mp - 1) * np + 2]);
#endif
                        scf[i * univlength + j] = psc[(mp - 1) * np + 2];
                    } /* The Sr right scaling value is applied */
                }
#ifdef DEBUG
                printf("\n");
#endif
            } // for j
        } // for i

        return 0;
    } // if mp > 1


    /* Only one element (single value or one point) */
    if (np == 1) { /* Single value */
        for (i = 0; i < su; i++) {
            scf[i] = psc[0];
        }
        return 0;
    }

    /* Single point */
    for (i = 0; i < numofunivs; i++) {

        for (j = 0; j < univlength; j++) {

            if (u[i * univlength + j] < psc[0]) { /* Before the given point (left side scaling) */
                scf[i * univlength + j] = psc[1];
                continue;
            }

            /* After the given point (right side scaling) */
            if (np == 2) {
                scf[i * univlength + j] = psc[i * np + 1]; /* Sl and Sr are the same */
            } else {
                scf[i * univlength + j] = psc[i * np + 2]; /* Sl and Sr are different */
            }

        }

    }

    return 0;
}
//#undef DEBUG
