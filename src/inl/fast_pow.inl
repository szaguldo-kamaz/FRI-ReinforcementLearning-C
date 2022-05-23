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


///calc pow faster in case of smaller power
///\param b base
///\param p power value
///return pow
static inline fri_float fast_pow(fri_float b, int p) {

#ifdef FAST_POW
    // TODO HACK
    long double ret = b;
//    double ret = b;

    for (int i = 0; i < p - 1; i++)
        ret *= b;

    return ret;
#else
    return pow(b, p);
#endif

}
