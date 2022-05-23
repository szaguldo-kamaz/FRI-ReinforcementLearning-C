/*
** Fuzzy Inference by Interpolation in Vague Environment toolbox for ANSI C
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** ANSI C / AVX port and FixRes optimization by David Vincze <david.vincze@webcode.hu>
** Various contibutions by Daniel Palko <palko.daniel@uni-miskolc.hu>
** mex port by Sandor Zsuga (zsuga@iit.uni-miskolc.hu)
** Original FIVE MATLAB code by Szilveszter Kovacs <szkovacs@iit.uni-miskolc.hu>
*/

///
///\file five_deinit.c
///

#include "FIVE.h"
#include <stdlib.h>

///deinit FIVE
///\param frb FIVE struct
///return void
void five_deinit(struct FIVERB *frb) {

    DEBUG_MSG("five_deinit: freeing FRB\n");

    for (unsigned int currant = 0; currant < (frb->rulelength - 1); currant++) {
        free(frb->rseqant[currant]);
        free(frb->rseqant_uindex[currant]);
        free(frb->rseqant_veval[currant]);
    }

    free(frb->fvc_vagdist);
    free(frb->frd_dists);
    free(frb->wi);
    free(frb->uk);
    free(frb->vek);
    free(frb->rseqant);
    free(frb->rseqant_uindex);
    free(frb->rseqant_veval);
//    free(frb->rant);
    free(frb->rant_uindex);
    free(frb->rant_veval);
//    free(frb->rconc);
    free(frb->ruledists);
    free(frb->weights);
    free(frb->ukdomains);
    free(frb->udivs);
    free(frb->valvagp);
    free(frb);

}
