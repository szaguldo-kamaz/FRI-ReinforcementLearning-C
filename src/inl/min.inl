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


///get array minimum value index
///\param arr array
///\param len length of the array
///\return minimum index
static inline unsigned int get_arr_min_i(fri_float *arr, int len) {
    int min_i = 0;

    for (int i = 1; i < len; i++) {
        if (arr[min_i] > arr[i]) {
            min_i = i;
        }
    }

    return min_i;

}


///get array absolute minimum value index
///\param arr array
///\param len length of the array
///\return absolute minimum index
static inline unsigned int get_arr_abs_min_i(fri_float *arr, int len) {

    //fri_float *tmp = MALLOC(len * sizeof(fri_float));
    fri_float tmp[len];

    arr_abs(arr, tmp, len);
    int min_i = get_arr_min_i(tmp, len);

    //free (tmp);

    return min_i;

}


///get vague environment array absolute minimum index
///\param universe subarray
///\param len lengt of the subarray
///\param point point value
///\return absolute minimum index
static inline unsigned int get_vag_abs_min_i(fri_float *universe, int len, fri_float point) {

    fri_float tmp[len];
    arr_sub_val(universe, tmp, len, point);

    return get_arr_abs_min_i(tmp, len);

}


///get vague environment array absolute minimum index (fixres)
///\param universe universe subarray
///\param len length of the subarray
///\param point point value
///\param universe_div universe subarray values division
///\return absolute minimum index
static inline unsigned int get_vag_abs_min_i_fixres(fri_float *universe, int len, fri_float point, fri_float arr_div) {

    int low_i = (point - *universe) / arr_div;

    if (low_i < 0) {
        return 0;
    } else {
        if (low_i >= len) {
            return len - 1;
        }
    }

    double d1 = universe[low_i] - point;
    double d2 = universe[low_i + 1] - point;

    if (fast_abs(d1) <= fast_abs(d2)) {
        return low_i;
    } else {
        return low_i + 1;
    }

}
