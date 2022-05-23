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

#ifndef FRIRL_APP_HELPERS_H
#define FRIRL_APP_HELPERS_H

#include "frirl_types.h"
#include "FIVE.h"
#include "frirl.h"

#ifdef BUILD_VISUALIZATION
#include "gui.h"
#endif


///allocate states / actions based on the _dimension_desc struct .values_len member value
#define FRIRL_ALLOC_VALUES(x)   fri_float _local_ ## x ## _values[x.values_len] /*SIMD_ALIGN*/;     x.values   = _local_ ## x ## _values;
///allocate universe based on frirl_dimension_desc struct .universe_len member value
#define FRIRL_ALLOC_UNIVERSE(x) fri_float _local_ ## x ## _universe[x.universe_len] /*SIMD_ALIGN*/; x.universe = _local_ ## x ## _universe;

///allocate states / actions and universe, based on the frirl_dimension_desc struct .values_len and .universe_len member values
#define FRIRL_ALLOC_DIM(x)  FRIRL_ALLOC_VALUES(x);  FRIRL_ALLOC_UNIVERSE(x);

///generalte linear universe values based on frirl_dimension_desc struct values
#define FRIRL_GEN_FIXRES_UNIVERSE(x) frirl_gen_fixres_arr(x.universe, x.universe_len, x.universe_div);
///generalte linear state / action values based on frirl_dimension_desc struct values
#define FRIRL_GEN_FIXRES_VALUES(x)   frirl_gen_fixres_arr(x.values,   x.values_len,   x.values_div);

///generate linear state / action and universe values based on frirl_statedim_desc struct values
#define FRIRL_GEN_FIXRES_DIM(x)  FRIRL_GEN_FIXRES_VALUES(x);  FRIRL_GEN_FIXRES_UNIVERSE(x);


void frirl_gen_fixres_arr(fri_float *arr, int len, fri_float div);

void frirl_run(struct frirl_desc *frirl, int verbose);

void frirl_visualization_init(struct frirl_desc *frirl);
void frirl_visualization_deinit();

#endif //FRIRL_APP_HELPERS_H
