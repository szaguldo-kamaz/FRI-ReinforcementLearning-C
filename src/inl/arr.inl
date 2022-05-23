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
#include <string.h>


///substract two arrays (out = arr1 - arr2)
///\param arr1 array one
///\param arr2 array two
///\parma len length of the arrays
static inline void arr_sub(fri_float *arr1, fri_float *arr2, fri_float *out, int len) {

    int i = 0;

/*#ifdef BUILD_AVX2
    // TODO TEST!
    int avxlen = len - len % 4;
    for (; i < avxlen; i += 4) {
        __m256d mm_arr1 = _mm256_loadu_pd(&(arr1[i]));  // unaligned?
        __m256d mm_arr2 = _mm256_loadu_pd(&(arr2[i]));  // unaligned?
        __m256d mm_sub = _mm256_sub_pd(mm_arr1, mm_arr2);
        _mm256_storeu_pd(&(out[i]), mm_sub);  // unaligned?
    }
#endif*/

    for (; i < len; i++) {
        out[i] = arr1[i] - arr2[i];
    }

}


///substract two arrays (out = arr - val)
///\param arr array
///\param val value
///\parma len length of the arrays
static inline void arr_sub_val(fri_float *arr, fri_float *out, int len, int val) {

    int i = 0;

/*#ifdef BUILD_AVX2
    // TODO TEST!
    int avxlen = len - len % 4;
    for (; i < avxlen; i += 4) {
        __m256d mm_arr = _mm256_loadu_pd(&(arr[i]));  // unaligned?
        __m256d mm_val = _mm256_set1_pd(val);
        __m256d mm_sub = _mm256_sub_pd(mm_arr, mm_val);
        _mm256_storeu_pd(&(out[i]), mm_sub);  // unaligned?
    }
#endif*/

    for (; i < len; i++) {
        out[i] = arr[i] - val;
    }

}


///calc array absolute values (out = fabs(arr))
///\param arr array
///\param out output array
///\param len length of the arrays
static inline void arr_abs(fri_float *arr, fri_float *out, int len) {

    int i = 0;

/*#ifdef BUILD_AVX2
    // TODO TEST!
    int avxlen = len - len % 4;
    for (; i < avxlen; i += 4) {
        __m256d mm_arr = _mm256_loadu_pd(&(arr[i]));  // unaligned?
        __m256d mm_abs = fast_abs_avx2(mm_arr);
        _mm256_storeu_pd(&(out[i]), mm_abs);
    }
#endif*/

    for (; i < len; i++) {
        out[i] = fast_abs(arr[i]);
    }

}


///sum the array values
///\param arr array
///\param len length of the array
///\return sum
static inline fri_float get_arr_sum(fri_float *arr, int len) {

    fri_float sum = 0.0;
    int i = 0;

// ezt igazabol frirl nem is hivja egyszer sem
/*#ifdef BUILD_AVX2
    // TODO TEST!
    int avxlen = len - len % 4;
    __m256d mm_sum = _mm256_setzero_pd();
    for (; i < avxlen; i += 4) {
//        __m256d mm_arr = _mm256_loadu_pd(&(arr[i])); // unaligned?
        __m256d mm_arr = _mm256_load_pd(&(arr[i]));
        _mm256_add_pd(mm_sum, mm_arr);
    }
    double_t *sum_mm = (double_t *)&mm_sum;
    sum = sum_mm[0] + sum_mm[1] + sum_mm[2] + sum_mm[3];
#endif*/

    for (; i < len; i++) {
        sum += arr[i];
    }

    return sum;

}


///count the inf elements (denoted by < 0!)
///\param arr array
///\param len length of the array
///\return count
static inline unsigned int get_arr_ltz_count(fri_float *arr, int len) {

    unsigned int count = 0;

    for (int i = 0; i < len; i++) {
        if (arr[i] < 0) {
            count++;
        }
    }

    return count;
}


///count the inf elements (denoted INFINITY!)
///\param arr array
///\param len length of the array
///\return count
static inline unsigned int get_arr_inf_count(fri_float *arr, int len) {

    unsigned int count = 0;

    for (int i = 0; i < len; i++) {
        if (arr[i] == INFINITY) {
            count++;
        }
    }

    return count;
}


///copy array
///\param src source array
///\param dest destionation array
///\param len length of the arrays
static inline void arr_cpy(fri_float *src, fri_float *dest, int len) {

    // TODO benchmark
    int i = 0;

/*#ifdef BUILD_AVX2
    int avxlen = len - len % 4;
    for (; i < avxlen; i += 4) {
        _mm256_storeu_pd(&(dest[i]), _mm256_loadu_pd(&(src[i])));
    }
    for (; i < len; i++) {
        dest[i]=src[i];
    }
#endif*/

    memcpy(dest, src, len * sizeof(fri_float));

}


///copy array, and who cares if it is overflowing
///\param src source array
///\param dest destionation array
///\param len length of the arrays
static inline void arr_cpy_dirty(fri_float *src, fri_float *dest, int len) {

    // TODO benchmark
/*#ifdef BUILD_AVX2
    int i = 0;
    for (; i < len; i += 4)
        _mm256_storeu_pd(&(dest[i]), _mm256_loadu_pd(&(src[i])));
#else*/

    memcpy(dest, src, len * sizeof(fri_float));

/*#endif*/

}
