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


#ifndef FRIRL_UTILS_H
#define FRIRL_UTILS_H

#include "frirl_types.h"


void frirl_show_rb(struct frirl_desc *frirl);
void frirl_show_hex_rb(struct frirl_desc *frirl);
int frirl_save_rb_to_text_file(struct frirl_desc *frirl, const char *file_name);
int frirl_save_rb_to_bin_file(struct frirl_desc *frirl, const char *file_name);
int frirl_load_rb_from_bin_file(struct frirl_desc *frirl, const char *file_name);
void frirl_print_usage();
void frirl_parse_cmdline(struct frirl_desc *frirl, int argc, char **argv);

#endif //FRIRL_UTILS_H
