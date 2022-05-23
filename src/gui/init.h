/*
** Fuzzy Rule Interpolation-based Reinforcement Learning (ANSI C / AVX version)
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** Author: David Vincze <david.vincze@webcode.hu>
**
** GUI code by Alex Toth <tothalex95@gmail.com>
**
*/


#ifndef INIT_H
#define INIT_H

#include <GL/freeglut.h>

void initGLUT(int *argc, char **argv, int width, int height, const char *title);
void deinitGLUT();

#endif
