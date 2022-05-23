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


#include "frirl.h"

///deinit FRIRL
///\param frirl FRIRL struct
///\return void
void frirl_deinit(struct frirl_desc *frirl) {

    DEBUG_MSG("frirl_deinit: freeing FIVERB\n");
    //todo, recheck what should be deleted

    free(frirl->fiverb->u);
    free(frirl->fiverb->ve);
    free(frirl->fiverb->rant);
    free(frirl->fiverb->rconc);
    five_deinit(frirl->fiverb);

    free(frirl->fus_proposed_values);
    free(frirl->fus_values);
    free(frirl->fus_check_states);
    free(frirl->fgba_statedistsum);
    free(frirl->fgba_vagdist_states);
    free(frirl->fgba_ruledist);
    free(frirl->fgba_actconc);
    free(frirl->fgba_dists);
    free(frirl->fep_ant);
    free(frirl->fep_cur_ant);
    free(frirl->fep_q_ant);
    free(frirl->fep_cur_q_ant);

    for (unsigned char i = 0; i < frirl->statedims_len; i++) {
        free(frirl->possible_states[i].values);
    }
    free(frirl->possible_states);

    free(frirl->possible_actions->values);
    free(frirl->possible_actions->vevalues);
    free(frirl->possible_actions);

    if (frirl->rbfile) {
        free(frirl->rbfile);
    }

#ifdef BUILD_MPI
    // Finalize the MPI environment.
    MPI_Finalize();
#endif

}
