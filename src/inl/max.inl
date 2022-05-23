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


///get array maximum value index
///\param arr array
///\len length of the array
static inline unsigned int get_arr_max_i(fri_float *arr, int len) {

    int max_i = 0;

    for (int i = 1; i < len; i++) {
        if (arr[max_i] < arr[i]) {
            max_i = i;
        }
    }

    return max_i;

}

///get array maximum value
///\param arr array
///\len length of the array
static inline fri_float get_arr_max(fri_float *arr, int len) {

    unsigned int max_i =  get_arr_max_i(arr, len);

    return arr[max_i];

}
