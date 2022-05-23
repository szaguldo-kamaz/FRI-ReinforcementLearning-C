/*
** Fuzzy Rule Interpolation-based Reinforcement Learning (ANSI C / AVX version)
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** Author: David Vincze <david.vincze@webcode.hu>
**
** Various contributions by Daniel Palko <palko.daniel@uni-miskolc.hu>
**
** Imitation related contributions by Alex Martossy
**
*/

#ifndef FRIRL_IMITATION_H
#define FRIRL_IMITATION_H

#include "frirl.h"


char getch();
void getActionFromTerminal(struct frirl_desc *frirl);

#endif
