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

#include "frirl_imitation.h"
#define FALSE_ACTION 0

#include <unistd.h>
#include <termios.h>
#include <stdio.h>


char getch() {

    char buf = 0;
    struct termios old = {0};

    if (tcgetattr(0, &old) < 0) {
        perror("tcsetattr()");
    }
    old.c_lflag &= ~ICANON;
    old.c_lflag &= ~ECHO;
    old.c_cc[VMIN] = 1;
    old.c_cc[VTIME] = 0;
    if (tcsetattr(0, TCSANOW, &old) < 0) {
        perror("tcsetattr ICANON");
    }
    if (read(0, &buf, 1) < 0) {
        perror ("read()");
    }
    old.c_lflag |= ICANON;
    old.c_lflag |= ECHO;
    if (tcsetattr(0, TCSADRAIN, &old) < 0) {
        perror ("tcsetattr ~ICANON");
    }

    return (buf);

}


void getActionFromTerminal(struct frirl_desc *frirl) {

    char c = 'v';

    c=getch();

    if (c == 'a') {
        frirl->keyaction = 0;
        frirl->valid_simulation = 1;
        return;
    }
    if (c == 's') {
        frirl->keyaction = 1;
        frirl->valid_simulation = 1;
        return;
    }
    if (c == 'd') {
        frirl->keyaction = 2;
        frirl->valid_simulation = 1;
        return;
    }
    if (c == 'p') {
        if (frirl->user_exited == 0) {
            frirl->user_exited = 1;
            printf("\nSimulation will end at the end of this episode.\n");
        }
        return;
    }
    if (c == 'u') {
        frirl->valid_simulation = 1;
        frirl->original_learning = 1;
        frirl->keyaction = 32;
        return;
    }

    //if (invalid key input || something goes wrong)
    frirl->valid_simulation = 0;

}
