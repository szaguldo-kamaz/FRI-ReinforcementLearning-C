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


///calc eucledian distance in case of infinite values, denoted by lesser than zero
///\param arr array
///\param len length of the array
///\return eucledian distance
static inline fri_float calc_eucledian_inf_distance(fri_float *arr, int len) {

    fri_float ret = 0.0;

    int i = 0;

/*#ifdef BUILD_AVX2
    // TODO TEST!
    __m256d mm_ret = _mm256_setzero_pd();
    int avxlen = len - len % 4;
    for(; i < avxlen; i += 4) {
        __m256d mm_arr = _mm256_loadu_pd(&(arr[i]));  // unaligned?
        __m256d mm_mask = _mm256_cmp_pd(mm_arr, _mm256_setzero_pd(), _CMP_LT_OQ); // arr < 0.0 ? T : F
        mm_arr = _mm256_and_pd(mm_arr, mm_mask);         //arr = mask ? arr : 0.0
        __m256d mm_tmp = _mm256_mul_pd(mm_arr, mm_arr);  //tmp = arr * arr
        mm_ret = _mm256_add_pd(mm_ret, mm_tmp);          //ret = ret + tmp
    }
    double_t *ret_mm = (double_t *)&mm_ret;
    ret = ret_mm[0] + ret_mm[1] + ret_mm[2] + ret_mm[3];
#endif*/

    for (; i < len; i++) {
        if (arr[i] < 0.0) {
            // there are inf elements
            ret += arr[i] * arr[i];
        }
    }

    return -fast_sqrt(ret);

}


///calc eucledian distance in case of positive values
///\param arr array
///\param len length of the array
///\return eucledian distance
static inline fri_float calc_eucledian_pos_distance(fri_float *arr, int len) {

// cartpole calls this 123 692 233x times

    fri_float ret = 0.0;

    int i = 0;

/*#ifdef BUILD_AVX2
    __m256d mm_ret = _mm256_setzero_pd();
    int avxlen = len - len % 4;
    for (; i < avxlen; i += 4) {
        __m256d mm_arr = _mm256_loadu_pd(&(arr[i]));  // unaligned?
        __m256d mm_tmp = _mm256_mul_pd(mm_arr, mm_arr);  // tmp = arr * arr
        mm_ret = _mm256_add_pd(mm_ret, mm_tmp);          // ret = ret + arr
    }
    double_t *ret_mm = (double_t *)&mm_ret;
    ret = ret_mm[0] + ret_mm[1] + ret_mm[2] + ret_mm[3];
#endif*/

// through the universe (1 rule - states+action)
    for (; i < len; i++) {
        ret += arr[i] * arr[i];
        //printf("rd: %d k: %d: dm[j]: %f %f\n", rni, k, dm[k], rd[rni]);
    }

    return fast_sqrt(ret);

}


///calc eucledian distance
///\param arr array
///\param len length of the array
///\return eucledian distance
static inline fri_float calc_eucledian_distance(fri_float *arr, int len) {

    int i = 0;

/*#ifdef BUILD_AVX2
    int avxlen = len - len % 4;
    for (; i < avxlen; i += 4) {
        __m256d mm_arr = _mm256_loadu_pd(&(arr[i]));  // unaligned?
        __m256d mm_cmp = _mm256_cmp_pd(mm_arr, _mm256_setzero_pd(), _CMP_LT_OQ);

        if (_mm256_movemask_pd(mm_cmp) != 0)
            return calc_eucledian_inf_distance(arr, len);
    }
#endif*/

// ezt valoszinuleg sosem hivja a FRIRL, nincs olyan helyzet amiben inf lenne - vannak sarokpontok, es kvantalas
// TODO: megnezni kiveheto-e, hogy tenyleg nem lesz sosem INF itt
    for (; i < len; i++) {
        if (arr[i] < 0.0) {
            return calc_eucledian_inf_distance(arr, len);
        }
    }

    return calc_eucledian_pos_distance(arr, len);

}
