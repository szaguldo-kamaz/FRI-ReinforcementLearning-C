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

#include <math.h>

static inline fri_float fast_sqrt(fri_float a) {

#if defined FAST_SQRT && defined BUILD_AVX2
    __asm__ __volatile__ (
        "vmovsd     %0,     %%xmm2         \n"
        "vsqrtsd %%xmm2, %%xmm1, %%xmm0 \n"
        "vmovsd     %%xmm0, %0             \n"
        : : "m" (a) : "xmm0", "xmm1", "xmm2", "memory"
    );

    return a;
#else
    return sqrt(a);
#endif

}
