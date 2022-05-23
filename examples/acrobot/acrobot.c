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
**
** Acrobot example
**  Based on original discrete version Q-learning demos programmed in MATLAB by: Jose Antonio Martin H. <jamartinh@fdi.ucm.es>
**
*/


#include "frirl.h"
#include "frirl_types_def.h"
#include "frirl_app_helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define USIZE 41
#define PI 3.14159265358979323846264338327


static inline void do_action(struct frirl_desc *frirl, fri_float action, fri_float *states, int states_len, fri_float *new_states) {

    // constants
    const fri_float maxSpeed1 = 4 * PI;
    const fri_float maxSpeed2 = 9 * PI;
    const fri_float m1        = 1.0;
    const fri_float m2        = 1.0;
    const fri_float l1        = 1.0;
    const fri_float l2        = 1.0;
    const fri_float l1Square  = l1 * l1;
    const fri_float l2Square  = l2 * l2;
    const fri_float lc1       = 0.5;
    const fri_float lc2       = 0.5;
    const fri_float lc1Square = lc1 * lc1;
    const fri_float lc2Square = lc2 * lc2;
    const fri_float I1        = 1.0;
    const fri_float I2        = 1.0;
    const fri_float g         = 9.8;
    const fri_float delta_t   = 0.05;

    // state
    fri_float theta1     = states[0];
    fri_float theta2     = states[1];
    fri_float theta1_dot = states[2];
    fri_float theta2_dot = states[3];


    fri_float d1 = m1 * lc1Square + m2 * (l1Square + lc2Square + 2 * l1 * lc2 * cos(theta2)) + I1 + I2;
    fri_float d2 = m2 * (lc2Square + l1 * lc2 * cos(theta2)) + I2;

    fri_float phi2 = m2 * lc2 * g * cos(theta1 + theta2 - PI / 2);
    fri_float phi1 = -m2 * l1 * lc2 * theta2_dot * sin(theta2) * (theta2_dot - 2 * theta1_dot) + (m1 * lc1 + m2 * l1) * g * cos(theta1 - (PI / 2)) + phi2;

    // 'action' means 'torque'
    fri_float accel1;
    fri_float accel2 = (action + phi1 * (d2 / d1) - m2 * l1 * lc2 * theta1_dot * theta1_dot * sin(theta2) - phi2);
    accel2 = accel2 / (m2 * lc2Square + I2 - (d2 * d2 / d1));
    accel1 = -(d2 * accel2 + phi1) / d1;

    //  Adam White's beginning of loop
    for (int i = 0; i < 4; i++) {
        theta1_dot = theta1_dot + accel1 * delta_t;
        if (theta1_dot < -maxSpeed1)
            theta1_dot = -maxSpeed1;

        if (theta1_dot > maxSpeed1)
            theta1_dot = maxSpeed1;

        theta1     = theta1 + theta1_dot * delta_t;
        theta2_dot = theta2_dot + accel2 * delta_t;

        if (theta2_dot < -maxSpeed2)
            theta2_dot = -maxSpeed2;

        if (theta2_dot > maxSpeed2)
            theta2_dot = maxSpeed2;

        theta2 = theta2 + theta2_dot*delta_t;
    }


    /*
        while(theta1<-PI) {
             theta1 = theta1 + 2*PI; 
         }
        while(theta1>PI) {
             theta1 = theta1 - 2*PI; 
        }
        while(theta2<-PI) {
             theta2 = theta2 + 2*PI; 
        }
        while(theta2>PI) {
             theta2 = theta2 - 2*PI; 
        }
     */


    if (theta1 < -PI) {
        theta1 = -PI;
    }

    if (theta1 > PI) {
        theta1 = PI;
    }

    if (theta2 < -PI) {
        theta2 = -PI;
    }

    if (theta2 > PI) {
        theta2 = PI;
    }

    DEBUG_MSG("theta1: %.18f theta2: %.18f theta1_dot: %.18f theta2_dot: %.18f\n", theta1, theta2, theta1_dot, theta2_dot);

    new_states[0] = theta1;
    new_states[1] = theta2;
    new_states[2] = theta1_dot;
    new_states[3] = theta2_dot;
}


static inline void get_reward(struct frirl_desc *frirl, fri_float *states, int states_len, struct frirl_reward_desc *reward) {

    // GetReward returns the reward at the current state
    // s: a vector of acrobot state
    // r: the returned reward.
    // f: true if the acrobot reached the goal, otherwise f is false

    fri_float theta1 = states[0];
    fri_float theta2 = states[1];
    fri_float y_acrobot[3];
    fri_float r = -10; // y_acrobot[2];
    fri_float f = 0;
    fri_float goal;


    y_acrobot[0] = 0.0;
    y_acrobot[1] = y_acrobot[0] - cos(theta1);
    y_acrobot[2] = y_acrobot[1] - cos(theta2);

    goal = y_acrobot[0] + 1.0;

    if (y_acrobot[2] >= goal) {
        r = 1000; // 10*y_acrobot[2]
        f = 1;
    }

    reward->value = r;
    reward->success = f;

}


static inline void quantize_observations(struct frirl_desc *frirl, fri_float *states, int states_len, fri_float *new_states) {

    int where;

    for (int i = 0; i < states_len; i++) {

        where = round((states[i] + fabs(frirl->statedims[i].values[0])) / frirl->statedims[i].values_div);
        /*DEBUG_MSG("frirl_quantize_observations: x[sno]: %.18f frirl->states[sno].statedim[0]: %.18f abs+: %.18f frirl->states[sno].statediv: %.18f befround: %.18f where: %d\n",
            states[i], frirl->statedim_arr[i].states[0], fabs(frirl->statedim_arr[i].states[0], frirl->statedim_arr[i].states_div,
            (states[i] + fabs(frirl->statedim_arr[i].states[0])) / frirl->statedim_arr[i].states_div, where));*/

        if (where < 0) {
            where = 0;

            DEBUG_MSG("frirl_quantize_observations: where<0 hit\n");

        } else if (where > frirl->statedims[i].values_len - 1) {
            where = frirl->statedims[i].values_len - 1;

            DEBUG_MSG("frirl_quantize_observations: where>statedims hit\n");
        }

        new_states[i] = frirl->statedims[i].values[where];

        DEBUG_MSG("frirl_quantize_observations: e/s: %d/%d where: %d: %d value: %.3f\n", 0, 0, i, where, frirl->statedims[i].values[where]);
    }

}


static inline void draw(struct frirl_desc *frirl, fri_float *new_curr_state, double action, unsigned int steps) {
#ifdef BUILD_VISUALIZATION
       // Ported from MATLAB to GLUT by Alex Toth <toth.alex@iit.uni-miskolc.hu>, 2018

    const double divider = 3;
    char text[32];
    double theta1 = new_curr_state[0];
    double theta2 = new_curr_state[1];


    sprintf(text, "Episode: %d Steps: %d", frirl->episode_num, steps);

    startDrawing();

    drawTitle(COLOR_BLACK, text);

    drawAxis(-2, 2, divider, 1);

    // Acrobot
    const double acrobotX1 = 0 / divider;
    const double acrobotY1 = 0 / divider;
    const double acrobotX2 = acrobotX1 + sin(theta1) / divider;
    const double acrobotY2 = acrobotY1 - cos(theta1) / divider;
    const double acrobotX3 = acrobotX2 + sin(theta2) / divider;
    const double acrobotY3 = acrobotY2 - cos(theta2) / divider;

    const double color[] = {.7, .7, .7};

    drawLine(COLOR_BLACK, acrobotX1, acrobotY1, acrobotX2, acrobotY2, 1);
    drawLine(COLOR_BLACK, acrobotX2, acrobotY2, acrobotX3, acrobotY3, 1);

    drawPoint(color, acrobotX1, acrobotY1, 20);
    drawPoint(color, acrobotX2, acrobotY2, 20);
    drawPoint(COLOR_RED, acrobotX3, acrobotY3, 20);

    endDrawing();
#endif // BUILD_VISUALIZATION

}


int main(int argc, char **argv) {

    struct frirl_desc frirl = frirl_desc_default;

    frirl_parse_cmdline(&frirl, argc, argv);

    frirl.reward_good_above = 0.0;

    frirl.alpha   = 0.5;
    frirl.gamma   = 1.0;
    frirl.epsilon = 0.001;

    frirl.qdiff_pos_boundary    = 1.0;
    frirl.qdiff_neg_boundary    = -200.0;
    frirl.qdiff_final_tolerance = 50.0;

    frirl.skip_rules = 1;  // Compatibility with the original MATLAB version

    frirl.do_action_func  = do_action;
    frirl.get_reward_func = get_reward;
    frirl.quant_obs_func  = quantize_observations;
//#ifdef BUILD_VISUALIZATION
    frirl.draw_func       = draw;
//#endif

    struct frirl_dimension_desc statedim0 = {
        .values           = (fri_float []) {-1.570796326794897, -0.785398163397448, 0, 0.785398163397448, 1.570796326794897},
        .values_len       = 5,
        .values_div       = 0.785398163397448,
        .values_steep = 1.0,
        .values_def   = 0.0,
        .universe_len     = USIZE,
        .universe_div     = 0.1
    };
    FRIRL_ALLOC_UNIVERSE(statedim0);
    FRIRL_GEN_FIXRES_UNIVERSE(statedim0);

    struct frirl_dimension_desc statedim1 = {
        .values           = (fri_float []) {-1.570796326794897, -0.785398163397448, 0, 0.785398163397448, 1.570796326794897},
        .values_len       = 5,
        .values_div       = 0.785398163397448,
        .values_steep = 1.0,
        .values_def   = 0.0,
        .universe_len     = USIZE,
        .universe_div     = 0.1
    };
    FRIRL_ALLOC_UNIVERSE(statedim1);
    FRIRL_GEN_FIXRES_UNIVERSE(statedim1);

    struct frirl_dimension_desc statedim2 = {
//        .values           = (frirl_float []) {-0.785398163397448, 0, 0.785398163397448},
        .values_len       = 3,
        .values_div       = 0.785398163397448,
        .values_steep = 1.0,
        .values_def   = 0.0,
        .universe_len     = USIZE,
        .universe_div     = 0.05
    };
//    FRIRL_ALLOC_UNIVERSE(statedim2);
//    FRIRL_GEN_FIXRES_UNIVERSE(statedim2);
    FRIRL_ALLOC_DIM(statedim2);
    FRIRL_GEN_FIXRES_DIM(statedim2);

    struct frirl_dimension_desc statedim3 = {
//        .values           = (frirl_float []) {-0.785398163397448, 0, 0.785398163397448},
        .values_len       = 3,
        .values_div       = 0.785398163397448,
        .values_steep = 1.0,
        .values_def   = 0.0,
        .universe_len     = USIZE,
        .universe_div     = 0.05
    };
//    FRIRL_ALLOC_UNIVERSE(statedim3);
//    FRIRL_GEN_FIXRES_UNIVERSE(statedim3);
    FRIRL_ALLOC_DIM(statedim3);
    FRIRL_GEN_FIXRES_DIM(statedim3);

    frirl.statedims = (struct frirl_dimension_desc []) {statedim0, statedim1, statedim2, statedim3};
    frirl.statedims_len = 4;

    struct frirl_dimension_desc actiondim = {
        .values_len  = 3,
        .values_div  = 1.0,
        .universe_len = USIZE,
        .universe_div = 0.1
    };
    FRIRL_ALLOC_DIM(actiondim);
    FRIRL_GEN_FIXRES_DIM(actiondim);
    frirl.actiondim = actiondim;

    // Initialization
    frirl_init(&frirl);

    // do the hustle
    frirl_run(&frirl, 1);

    // display the constructed rule-base
//    printf("\nIncrementally constructed rulebase:\n");
//    frirl_show_rb(&frirl);
    // or with precise hexadecimal float
    //frirl_show_rb_hex(&frirl);

    // dump the RB to files
    frirl_save_rb_to_bin_file(&frirl,"acrobot.frirlrb.bin");
    frirl_save_rb_to_text_file(&frirl,"acrobot.frirlrb.txt");

    // deinit
    frirl_deinit(&frirl);

}
