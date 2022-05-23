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

///swap if first greater than second
///\param i first variable
///\param j second variable
static inline void sort_swap(int *i, int *j) {

    if (*j < *i) {
        int tmp = *i;
        *i = *j;
        *j = tmp;
    }

}
