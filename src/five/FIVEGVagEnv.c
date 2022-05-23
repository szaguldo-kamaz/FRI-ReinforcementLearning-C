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
** FIVEGVagEnv:  Generate the vague environment VE (the scaled distances of the U universe)
**
**                            [VE]=FIVEGVagEnv(U,SCF)
**
**          Where:
**          U: is the universe (a vector of discrete values in increasing order),
**          SCF: is the scaling function 
**               (scaling factors according to the elements of U)
**          VE: is the generated vague environment (the scaled distances on the universe U)
**             If VE(k,1)=0, it is the primitive integral of the scaling function according
**               to the elements of U
**             If VE(k,1)=-1, it is the integral of the scaling function between
**               the neighboring elements of U
**          In case of U,SCF,VE the rows are the dimensions
**
*/

//#define DEBUG

#include "FIVE.h"

#include <stdlib.h>
#include <math.h>
#ifdef DEBUG
#include <stdio.h>
#endif


double *FIVEGVagEnv(double *u, int numofunivs, int univlength, double *scf) {

    int    su = numofunivs * univlength;
    int    i, k, j;
    double *ve;


    ve = MALLOC(su * sizeof(double));
    if (ve == NULL) {
        return NULL;
    }

    for (k = 0; k<numofunivs; k++) {

        for (i = k*univlength; i < (k + 1)*univlength; i++) {
            if (scf[i] == INFINITY) {
                break;
            }
        }

#ifdef DEBUG
        printf("GVagEnv: k: %d, i: %d, (k + 1)*univlength: %d\n", k, i, (k + 1)*univlength);
#endif
        if (i < ((k + 1)*univlength)) { // An inf scaling factor exists

            ve[k*univlength] = -1; // First

//VE(k,2:n) = diff(U(k,:)).*(SCF(k,1:n-1) + SCF(k,2:n))/2;
            for (j = k*univlength; j < (k + 1)*univlength-1; j++) { /* Trapezoidal area between neighbors */
//                ve[j + 1] = fabs(u[j + 1]-u[j])*(scf[j] + scf[j + 1])*0.5;
#ifdef DEBUG
                printf("GVagEnv: j: %d, ve[j]: %f, u[j + 1]-u[j]: %f scf[j]: %f scf[j + 1]: %f\n", j, ve[j], u[j + 1] - u[j], scf[j], scf[j + 1]);
#endif
                ve[j + 1] = (u[j + 1] - u[j]) * (scf[j] + scf[j + 1]) * 0.5;
            }

//        for (j = k + numofunivs; j<su; j+=numofunivs) /* Trapezoidal area between neighbors */
//            ve[j] = fabs(u[j]-u[j-numofunivs])*(scf[j] + scf[j-numofunivs])*0.5;

        } else {

            ve[k * univlength] = 0;  /* First element */

//e = diff(U(k,:)).*(SCF(k,1:n-1) + SCF(k,2:n))/2;
//for i = 1:n-1
//    VE(k,i + 1) = VE(k,i) + e(i); % primitive integral of the scaling function
//end % for i

            for (j = k*univlength; j<(k + 1)*univlength-1; j++) { // Primitive integral of the scaling function
#ifdef DEBUG
                printf("GVagEnv: j: %d, ve[j]: %f, u[j + 1]-u[j]: %f scf[j]: %f scf[j + 1]: %f\n", j, ve[j], u[j + 1] - u[j], scf[j], scf[j + 1]);
#endif
                ve[j + 1] = ve[j] + (u[j + 1]-u[j])*(scf[j] + scf[j + 1]) * 0.5;
//                 ve[j + 1] = (u[j] - u[j - numofunivs])*(scf[j] + scf[j - numofunivs]) * 0.5 + ve[j - numofunivs];
            }

        }

    }

    return ve;

}
//#undef DEBUG
