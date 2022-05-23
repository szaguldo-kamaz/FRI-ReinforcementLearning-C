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
**
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


double FIVEVagConcl(struct FIVERB *frb, double *x) {

    double conc;
    FIVE_vag_concl(frb, x, &conc);
    return conc;

}


unsigned int FIVE_vag_concl(struct FIVERB *frb, double *ant, double *conc) {

    unsigned int i;

    double vagc = 0; // Same as vagc in MatLab code
    double vagcz;
    double ws, wi;   // Same as ws and wi in MatLab code
    double wsz;

    double *vagdist = frb->fvc_vagdist;
    int infc = 0;  // Count of 'inf's
    double nc, zc;
    double y, *yp; // the return value


    /* Calculate the distances of the observation from the rule antecendents */
    int frd_ret = five_rule_distance(frb, ant); // result distance go into frb->ruledists

    /* Go through all rules; R(i,m) is the i. conclusion;
     **    rd(i) is the i. distance of the i. antecendent to the observation */

    /* Note: I use a different approach for working with inf distances than how
     **       it is done in the original MatLab code. Instead of filling an array
     **       and analyzing it later, i switch to infinite summing when hitting
     **       an 'inf' - Zsuga */

    // cartpole    comes here 91911x times
    // acrobot     comes here 44996x times
    // mountaincar comes here 45210x times

#if defined(FRIRL_FAST) && defined(FIVE_NONAN)
    // with FRIRL - NANs are not allowed in the rule-base, so only one rule can match, hence if the first hit found, we're done
    if (frd_ret != -1) {
        *conc = frb->rconc[frd_ret]; // not -1 means that an exact hit was found - carpotle 89963x - acrobot 21944x - mountaincar 44780x
        return frd_ret;
    }
    for (i = 0; i < frb->numofrules; i++) {
        if (frb->ruledists[i] == 0.0) {
            *conc = frb->rconc[i];
            return i; // If i<numofrules then there was a rule which hit exactly
        }
    }
#else
     // If i<numofrules then there was a rule which hit exactly
    for (i = 0; i < frb->numofrules; i++) {
        if (frb->ruledists[i] == 0.0) {
            break;
        }
    }

    if (i < frb->numofrules) { // if sum(frb->ruledists==0)>0 - at least one rule antecedent had exactly hit
#ifdef DEBUG
/*        if ( i != frd_ret ) {
            printf("bug %d\n",i);
            exit(5);
        }*/
        printf("VagConcl: exact hit: i: %d vs %d\n", i,frd_ret);
#endif
        nc = 0; // num of consequents gathered
        zc = 0; // Counts zeros until reaching first 'inf'

        for (i = 0; i < frb->numofrules; i++) { // for all the rules
#ifdef DEBUG
            printf("frb->ruledists[%d]: %.18f\n", i, frb->ruledists[i]);
#endif
            if (frb->ruledists[i] == 0.0) { // rule distance zero -> hit

                if (frb->rulelength != frb->univlength) { // The rule consequences are singletons; do finite summing only
                    vagc += frb->rconc[i]; // value of rule consequent
                    nc++;
#ifdef DEBUG
                    //todo rb removed
//                    printf("VagConcl: singleton consequence: i: %d, frb->rb[i*frb->rulelength+(frb->rulelength-1)]: %f vagc: %f\n", i, frb->rb[i * frb->rulelength + (frb->rulelength - 1)], vagc);
                    printf("VagConcl: singleton consequence: i: %d, frb->rconc[i]: %f vagc: %f\n", i, frb->rconc[i], vagc);
#endif
                    continue;
                }

                // else is here

                // The consequent universe also has a vague environment
                // Vague distance from the first (smallest) element of the conclusion universe
                // FIVEVagDist(U(m,:),VE(m,:),U(m,1),R(i,m))]

                // params: frb, p1 (first element of consequent U), p2 (rule-base, consequent of i-th rule), retval[numofunivs]
                five_vague_distance(frb, frb->uk[frb->numofunivs - 1], &frb->rconc[i], vagdist);

                // if min(vrconc)<0    % at least one rule consequent is in inf distance
                //   vagc=sum(vrconc(vrconc<=0))./size(vrconc(vrconc<=0),2);

                if (infc) { // Do inf summing
                    if (*vagdist <= 0.0) {
                        vagc += *vagdist;
                        nc++;
                    }
#ifdef DEBUG
                    printf("VagConcl: consequence has VE: inf sum: vagc: %f\n", vagc);
#endif
                    continue;
                }
                // else is here

                if (*vagdist >= 0.0) { // Do finite summing
                    if (*vagdist == 0.0) {
                        zc++;
                    } // why add when it's zero? should put an "else" here
                    vagc += *vagdist;

                    nc++;
#ifdef DEBUG
                    printf("VagConcl: consequence has VE: fin sum: vagc: %f\n", vagc);
#endif
                    continue;
                }

                // Switch to inf summing (vagdist<0 with infc==0) 
                vagc = *vagdist;
                nc = zc + 1; // Add the count of zeros passed so far
                infc = -1;

            } // if frb->ruledists[i]==0.0

        } // for all rules

        vagc /= nc; /* Let the conclusion be the average of the conclusions exactly hit */

    } else { // Shepard interpolation
#endif // defined(FRIRL_FAST) && defined(FIVE_NONAN)

        // cartpole comes here 1948x times
        // acrobot  comes here 23052x times
#ifdef DEBUG
        printf("VagConcl: no exact hit!\n");
#endif

        ws = 0;
        wsz = 0;
        vagcz = 0;  // Sum zero elements seperately to include them later in inf calculation

#ifndef FIVE_NOINF
        // search for not inf hits
        for (i = 0; i < frb->numofrules; i++) {
            if (frb->ruledists[i] > 0.0) {
                break;
            }
        }
    #ifdef DEBUG
        for (i = 0; i < frb->numofrules; i++) {
            printf("frb->ruledists[%d]: %.18f\n", i, frb->ruledists[i]);
        }
    #endif
        if (i < frb->numofrules) {  // At least one rule antecendent is not in inf distance
            for (i = 0; i < frb->numofrules; i++) {
                if (frb->ruledists[i] < 0.0) {
                    frb->ruledists[i] = INFINITY;
                }
            }
        }
#endif

        for (i = 0; i < frb->numofrules; i++) {  // for all the rules

            wi = 1.0 / fast_pow(fast_abs(frb->ruledists[i]), frb->p);
#ifdef DEBUG
            printf("VagConcl: singleton consequence: i: %d, wi: %f, frb->ruledists[i]: %f p: %d fabs(frb->ruledists[i]): %f pow(fabs(frb->ruledists[i]),p): %f\n", i, wi, frb->ruledists[i], frb->p, fabs(frb->ruledists[i]), pow(fabs(frb->ruledists[i]), frb->p));
#endif

#ifndef FRIRL_FAST // always true for FRIRL
            if (frb->rulelength != frb->univlength) { // The rule consequences are singletons
#endif
                vagc += wi * frb->rconc[i]; // value of rule consequent
                ws += wi;
#ifdef DEBUG
                //todo rb removed
//                printf("VagConcl: singleton consequence2: i: %d, wi: %f, frb->rb[i*frb->rulelength+(frb->rulelength-1)]: %f vagc: %f\n", i, wi, frb->rb[i * frb->rulelength + (frb->rulelength - 1)], vagc);
                printf("VagConcl: singleton consequence2: i: %d, wi: %f, frb->rconc[i]: %f vagc: %f\n", i, wi, frb->rconc[i], vagc);
#endif
                continue;
#ifndef FRIRL_FAST
            }
#endif

            // The consequent universe also has a vague environment
            // Vague distance from the first (smallest) element of the conclusion universe

            // FIVEVagDist(U(m,:),VE(m,:),U(m,1),R(i,m))
#ifdef DEBUG
            //todo rb removed
//            printf("VagConcl: consequent has VE: execute: FIVEVagDist frb->u[%d], 1, %d, frb->ve[%d], %f, %f\n", (frb->rulelength - 1) * frb->numofunivs, frb->numofunivs, (frb->rulelength - 1) * frb->numofunivs, frb->u[(frb->rulelength - 1) * frb->numofunivs], frb->rb[(i + 1) * frb->numofrules - 1]);
            printf("VagConcl: consequent has VE: execute: FIVEVagDist frb->u[%d], 1, %d, frb->ve[%d], %f, %f\n", (frb->rulelength - 1) * frb->numofunivs, frb->numofunivs, (frb->rulelength - 1) * frb->numofunivs, frb->u[(frb->rulelength - 1) * frb->numofunivs], frb->rconc[i]);
#endif

            // params: frb, p1 (first element of consequent U), p2 (rule-base, consequent of i-th rule), retval[numofunivs]
            five_vague_distance(frb, frb->uk[frb->numofunivs - 1], &frb->rconc[i], vagdist);

#ifdef DEBUG
            printf("VagConcl: FIVEVagDist result: vagdist: %f\n", *vagdist);
#endif
            if (infc) { // Do inf summing
                if (*vagdist <= 0.0) {
                    vagc += *vagdist * wi;
                    ws += wi;
                }
#ifdef DEBUG
                printf("VagConcl: consequence has VE: inf sum: vagc: %f\n", vagc);
#endif
                continue;
            }

            if (*vagdist >= 0.0) { // Do finite summing
                if (*vagdist == 0.0) {
                    // if vagdist is zero then there's no need to multiply and accumulate
                    // vagcz += *vagdist*wi;
                    wsz += wi;
                } // if vagdist is zero then ther's no need to mul+acc, so an "else" should be here?
                vagc += *vagdist*wi;

                ws += wi;
#ifdef DEBUG
                printf("VagConcl: consequence has VE: fin sum: vagc: %f vagcz: %f\n", vagc, vagcz);
#endif
                continue;
            }

            // Switch to inf summing
            vagc = vagcz + *vagdist * wi;
            ws = wsz + wi;
            infc = -1;
        }

#ifndef FIVE_NOINF
        if (ws == 0.0) { // There is no valid conclusion (infinite distances from all the rules)
            vagc = NAN;
    #ifdef DEBUG
            printf("VagConcl: no valid conclusion: vagc: %f\n", vagc);
    #endif
        } else {
#endif
            vagc /= ws; /* The valid interpolated conclusion */
#ifdef DEBUG
            printf("VagConcl: valid interpolated conclusion: vagc/ws: %f ws: %f\n", vagc, ws);
#endif
#ifndef FIVE_NOINF
        }
#endif

#if !(defined(FRIRL_FAST) && defined(FIVE_NONAN))
    }
#endif

#ifndef FIVE_NONAN
    if (isnan(vagc)) {
        y = NAN;
    #ifdef DEBUG
        printf("VagConcl: y=NAN");
    #endif
    } else {
#endif
        if (frb->rulelength == frb->univlength) { // The consequent universe also has a vague environment
            //        Y=FIVEValVag(U(m,:),VE(m,:),vagc);  % The valid conclusion
#ifdef DEBUG
            printf("VagConcl: execute FIVEValVag(frb,%f)\n", vagc);
#endif
            if (FIVEValVag(frb, &vagc) == -1) {
#ifdef DEBUG
                printf("VagConcl: FIVEValVag error\n");
#endif
                return -1;
            }
            y = *frb->valvagp;
#ifdef DEBUG
            printf("VagConcl: consequence has VE: ValVag: %f\n", y);
#endif
        } else { // The rule consequences are singletons
            y = vagc;
#ifdef DEBUG
            printf("VagConcl: consequence singleton: y=vagc: %f\n", y);
#endif
        }
#ifndef FIVE_NONAN
    }
#endif

    *conc = y;

    return ~0;

}
