/*
** Fuzzy Rule Interpolation-based Reinforcement Learning (ANSI C / AVX version)
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** Author: David Vincze <david.vincze@webcode.hu>
**
** Various contributions by Daniel Palko <palko.daniel@uni-miskolc.hu>
**
** GUI code by Alex Toth <tothalex95@gmail.com>
**
** Imitation related contributions by Alex Martossy
**
** Mountain-car example
**  Based on original discrete version Q-learning demos programmed in MATLAB by: Jose Antonio Martin H. <jamartinh@fdi.ucm.es>
**
*/


#include "frirl.h"
#include "frirl_types.h"
#include "frirl_types_def.h"
#include "frirl_app_helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef BUILD_VISUALIZATION
#include "../../src/gui/gui.h"
#endif

#define USIZE 41
#define PI 3.14159265358979323846264338327


static inline void do_action(struct frirl_desc *frirl, fri_float action, fri_float *states, int states_len, fri_float *new_states) {
    // MountainCarDoAction: executes the action (a) into the mountain car environment
    // a: is the force to be applied to the car
    // x: is the vector containning the position and speed of the car
    // xp: is the vector containing the new position and velocity of the car

    // Parameters for the simulation

    fri_float position = states[0];
    fri_float speed    = states[1];

    const fri_float bpleft  = -1.5;
    const fri_float bsleft  = -0.07;
    const fri_float bsright = +0.07;

    // thermodynamic law, for a more real system with friction.
    fri_float speedt1 = (speed + (0.001 * action) + (-0.0025 * cos(3.0 * position))) * 0.999;
    fri_float post1;

    if (speedt1 < bsleft) {
        speedt1 = bsleft;
    }
    if (speedt1 > bsright) {
        speedt1 = bsright;
    }

    post1 = position + speedt1;

    if (post1 <= bpleft) {
        post1 = bpleft;
        speedt1 = 0.0;
    }

    DEBUG_MSG("speedt1: %.18f post1: %.18f\n", speedt1, post1);

    new_states[0] = post1;
    new_states[1] = speedt1;
}


static inline void get_reward(struct frirl_desc *frirl, fri_float *states, int states_len, struct frirl_reward_desc *reward) {

    fri_float position = states[0];

    fri_float r,f;

    const fri_float bpright = 0.45;

    r = -10;
    f = 0;

    if ( position >= bpright ) {
        r = 1000;
        f = 1;
    }

    reward->value   = r;
    reward->success = f;
}


static inline void quantize_observations(struct frirl_desc *frirl, fri_float *states, int states_len, fri_float *new_states) {

    int where;

    for (int i = 0; i < states_len; i++) {

        where = round((states[i] + fabs(frirl->statedims[i].values[0])) / frirl->statedims[i].values_div);
        /*DEBUG_MSG("frirl_quantize_observations: x[sno]: %.18f frirl->states[sno].statedim[0]: %.18f abs+: %.18f frirl->states[sno].statediv: %.18f befround: %.18f where: %d\n",
            states[i], frirl->statedim_arr[i].states[0], fabs(frirl->statedim_arr[i].states[0], frirl->statedim_arr[i].states_div,
            (states[i] + fabs(frirl->stated_arr[i].states[0])) / frirl->statedim_arr[i].states_div, where));*/

        if (where < 0) {
            where = 0;

            DEBUG_MSG("frirl_quantize_observations: where<0 hit\n");

        } else if (where > frirl->statedims[i].values_len - 1) {
            where = frirl->statedims[i].values_len - 1;

            DEBUG_MSG("frirl_quantize_observations: where>statedims hit\n");
        }

        new_states[i] = frirl->statedims[i].values[where];

        DEBUG_MSG("frirl_quantize_observations: e/s: %d/%d where: %d: %d value: %.3f\n", 0,0, i, where, frirl->statedims[i].values[where]);
    }

}


#ifdef BUILD_VISUALIZATION
static inline void draw(struct frirl_desc *frirl, fri_float *new_curr_state, double action, unsigned int steps) {
    // Ported from MATLAB to GLUT by Alex Toth <toth.alex@iit.uni-miskolc.hu>, 2018

    const double divider = 1;
    const double borderColor[] = {0, 0, 0};
    char text[32];
    double x = new_curr_state[0];

    sprintf(text, "Episode: %d Steps: %d", frirl->episode_num, steps);

    startDrawing();

    drawTitle(COLOR_BLACK, text);

    // Mountain
    double i;
    const double areaColor[] = {.1, .7, .1};
    double areaX1, areaY1, areaX2, areaY2;

    for (i = -1.6; i <= 0.6; i += 0.005)
    {
        areaX1 = i;
        areaY1 = sin(3 * areaX1);
        areaX2 = i + 0.005;
        areaY2 = sin(3 * areaX2);

        drawFilledRectangle(areaColor, areaX1 + .5, areaY1 / 2, areaX2 + .5, -1.0);
    }

    // Car
    const double carColor[] = {1, .7, .1};
    const double carX1 = x - 0.075;
    const double carY1 = sin(3 * (x - 0.075)) + 0.3;
    const double carX2 = x + 0.075;
    const double carY2 = sin(3 * (x + 0.075)) + 0.3;

    drawLine(carColor, carX1 + .5, carY1 / 2, carX2 + .5, carY2 / 2, 50);
    // drawFilledRectangle(carColor, carX1 + .5, carY1 / 2, carX2 + .5, carY2 / 2);
    // drawRectangle(borderColor, carX1 + .5, carY1 / 2, carX2 + .5, carY2 / 2, 1);

    // Wheels
    const double wheelColor[] = {.5, .5, .5};
    const double wheelX1 = x - 0.05;
    const double wheelY1 = sin(3 * (x - 0.05)) + 0.06;
    const double wheelX2 = x + 0.05;
    const double wheelY2 = sin(3 * (x + 0.05)) + 0.06;

    drawPoint(wheelColor, wheelX1 + .5, wheelY1 / 2, 50);
    drawPoint(wheelColor, wheelX2 + .5, wheelY2 / 2, 50);

    // Goal
    const double goalColor[] = {1, .7, .1};
    const double goalX = 0.45;
    const double goalY = sin(3 * 0.5) + 0.1;

    drawFilledAsterisk(goalColor, goalX + .5, goalY / 2, 5, .05);
    drawAsterisk(borderColor, goalX + .5, goalY / 2, 5, 1, .05);

    // Direction of the force
    const double forceColor[] = {0, 1, 0};
    double forceX, forceY;

    if (action < 0) {
        forceX = x - 0.08 - 0.1;
        forceY = sin(3 * (x - 0.05)) + 0.3;

        drawFilledTriangle(forceColor, forceX - 0.025 + .5, forceY / 2, forceX + 0.05 + .5, (forceY - 0.1) / 2, forceX + 0.05 + .5, (forceY + 0.1) / 2);
        drawTriangle(borderColor, forceX - 0.025 + .5, forceY / 2, forceX + 0.05 + .5, (forceY - 0.1) / 2, forceX + 0.05 + .5, (forceY + 0.1) / 2, 1);
    } else {
        if (action > 0) {
            forceX = x + 0.08 + 0.1;
            forceY = sin(3 * (x + 0.05)) + 0.3;

            drawFilledTriangle(forceColor, forceX + 0.025 + .5, forceY / 2, forceX - 0.05 + .5, (forceY - 0.1) / 2, forceX - 0.05 + .5, (forceY + 0.1) / 2);
            drawTriangle(borderColor, forceX + 0.025 + .5, forceY / 2, forceX - 0.05 + .5, (forceY - 0.1) / 2, forceX - 0.05 + .5, (forceY + 0.1) / 2, 1);
        }
    }

    drawAxis(-1, 1, divider, 1);

    endDrawing();
}
#endif // BUILD_VISUALIZATION


int main(int argc, char **argv) {

    struct frirl_desc frirl = frirl_desc_default;

    frirl_parse_cmdline(&frirl, argc, argv);

    frirl.reward_good_above = -5000.0;

    frirl.alpha   = 0.5;
    frirl.gamma   = 1.0;
    frirl.epsilon = 0.01;

    frirl.qdiff_pos_boundary    = 1.0;
    frirl.qdiff_neg_boundary    = -4.0;
    frirl.qdiff_final_tolerance = 500.0;

    frirl.skip_rules = 1;  // Compatibility with the original MATLAB version

    frirl.do_action_func  = do_action;
    frirl.get_reward_func = get_reward;
    frirl.quant_obs_func  = quantize_observations;
#ifdef BUILD_VISUALIZATION
    frirl.draw_func       = draw;
    frirl.original_learning = 0;
#endif

    frirl.construct_rb = 0;
    frirl.reduce_rb = 1;
    frirl.reduction_strategy = 1;

    struct frirl_dimension_desc statedim0 = {
        .values           = (fri_float []) {-1.5, -1.295, -1.09, -0.885, -0.68, -0.475, -0.27, -0.065, 0.14, 0.345},
        .values_len       = 10,
        .values_div       = 0.205,
        .values_steep     = 1.0,
        .values_def       = -0.5,
        .universe_len     = USIZE,
        .universe_div     = 0.1
    };
    FRIRL_ALLOC_UNIVERSE(statedim0);
    FRIRL_GEN_FIXRES_UNIVERSE(statedim0);

    struct frirl_dimension_desc statedim1 = {
        .values           = (fri_float []) {-0.07, -0.042, -0.014, 0.014, 0.042, 0.07},
        .values_len       = 6,
        .values_div       = 0.028,
        .values_steep     = 1.0,
        .values_def       = 0.0,
        .universe_len     = USIZE,
        .universe_div     = 0.005
    };
    FRIRL_ALLOC_UNIVERSE(statedim1);
    FRIRL_GEN_FIXRES_UNIVERSE(statedim1);

    frirl.statedims = (struct frirl_dimension_desc []) {statedim0, statedim1};
    frirl.statedims_len = 2;

    struct frirl_dimension_desc actiondim = {
        .values_len   = 3,
        .values_div   = 1.0,
        .universe_len = USIZE,
        .universe_div = 0.1
    };
    FRIRL_ALLOC_DIM(actiondim);
    FRIRL_GEN_FIXRES_DIM(actiondim);
    frirl.actiondim = actiondim;

    // Initialization
    frirl_init(&frirl);

    // help for imitation
    if (frirl.original_learning == 0) {
        printf("Press the \'u\' key for NO imitation. (Recommended for no GUI)\nOtherwise use 'a-s-d' for left-nothing-right actions.\nPress 'p' to exit and save at the end of current episode.\n");
    }

    frirl_run(&frirl, 1);

//    printf("\nIncrementally constructed rulebase:\n");
//    frirl_show_rb(&frirl);
    // or with precise hexadecimal float
    //frirl_show_rb_hex(&frirl);

    // dump the RB to files
    frirl_save_rb_to_bin_file(&frirl,"mountaincar.frirlrb.bin");
    frirl_save_rb_to_text_file(&frirl,"mountaincar.frirlrb.txt");

    // Destroy
    frirl_deinit(&frirl);

//#ifdef BUILD_VISUALIZATION
//    // Deinit visualization
//    deinitGLUT();
//#endif

}
