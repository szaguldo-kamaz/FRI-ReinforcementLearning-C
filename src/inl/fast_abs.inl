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


#include <stdint.h>
#include <math.h>

///calc absoulute value faster
///\param a value
///\return absolute value
static inline fri_float fast_abs(fri_float a) {

/*#if defined(FAST_ABS) && defined(BUILD_AVX2)
//    const __m256d mm_mask = _mm256_set1_pd(-0.0);
//    __m256d mm_a = _mm256_set1_pd(a);

//    __m256d mm_abs   = _mm256_andnot_pd(mm_mask, mm_a);
    __m256d mm_abs   = _mm256_andnot_pd(_mm256_set1_pd(-0.0), _mm256_set1_pd(a));

    double_t *abs_mm = (double_t *)&mm_abs;

    return abs_mm[0];
#elif defined(FAST_ABS)*/

#ifdef FAST_ABS
    union {
        double ret;
        uint64_t i;
    } u = {a};

    u.i &= 0x7fffffffffffffff;

    return u.ret;
#else
    return fabs(a);
#endif

}

/*#ifdef BUILD_AVX2
static inline __m256d fast_abs_avx2(__m256d mm_a) {

#ifdef FAST_ABS
//    const __m256d mm_mask = _mm256_set1_pd(-0.0);

//    return _mm256_andnot_pd(mm_mask, mm_a);
    return _mm256_andnot_pd(_mm256_set1_pd(-0.0), mm_a);
#else
    double *a_mm = (double_t *)&mm_a;

    return _mm256_set_pd(fabs(a_mm[3]), fabs(a_mm[2]), fabs(a_mm[1]), fabs(a_mm[0]));
#endif
}
#endif
*/
