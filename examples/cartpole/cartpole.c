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
** Cart-pole example
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

#ifdef BUILD_VISUALIZATION
#include "../../src/gui/gui.h"
#endif

#define USIZE 1001
#define PI 3.14159265358979323846264338327


static inline void do_action(struct frirl_desc *frirl, fri_float action, fri_float *states, int states_len, fri_float *new_states) {

    // Cart_Pole:  Takes an action (-10 to +10) and the current values of the
    // four state variables and updates their values by estimating the state
    // TAU seconds later.

    // Parameters for the simulation

    fri_float x          = states[0];
    fri_float x_dot      = states[1];
    fri_float theta      = states[2];
    fri_float theta_dot  = states[3];

    fri_float const g               = 9.8;      // Gravity
    fri_float const cart_mass       = 1.0;      // Mass of the cart is assumed to be 1Kg
    fri_float const pole_mass       = 0.1;      // Mass of the pole is assumed to be 0.1Kg
    fri_float const total_mass      = cart_mass + pole_mass;
    fri_float const len             = 0.5;      // Half of the length of the pole
    fri_float const pole_mass_len   = pole_mass * len; // max = 0.05
    fri_float const mag_force       = 10.0;
    fri_float const tau             = 0.02;     // Time interval for updating the values
    fri_float const fourthirds      = 4.0/3.0;

    fri_float force = action * mag_force; // max = 10

    fri_float temp     = (force + pole_mass_len * theta_dot * theta_dot * sin(theta)) / total_mass;
    fri_float thetaacc = (g * sin(theta) - cos(theta) * temp) / (len * (fourthirds - pole_mass * cos(theta) * cos(theta) / total_mass));
    fri_float xacc     = temp - pole_mass_len * thetaacc * cos(theta) / total_mass;

    // Update the four state variables, using Euler's method.

    DEBUG_MSG("x: %f x_dot: %f theta: %f theta_dot: %f\n",x,x_dot,theta,theta_dot);

    x         = x + tau * x_dot;
    x_dot     = x_dot + tau * xacc;
    theta     = theta + tau * theta_dot;
    theta_dot = theta_dot + tau * thetaacc;

    new_states[0] = x;
    new_states[1] = x_dot;
    new_states[2] = theta;
    new_states[3] = theta_dot;
}

static inline void get_reward(struct frirl_desc *frirl, fri_float *states, int states_len, struct frirl_reward_desc *reward) {

    double x         = states[0];
    //double x_dot     = states[1];
    double theta     = states[2];
    double theta_dot = states[3];

    double r, f;

    const double twelve_degrees     = PI/15;
    const double fourtyfive_degrees = PI/4;

    DEBUG_MSG("rew: x: %f theta: %f theta_dot: %f\n",x,theta,theta_dot);
    DEBUG_MSG("fabs(10*theta): %f pow(fabs(10*theta),2): %f\n",fabs(10*theta),pow(fabs(10*theta),2));

//    if ((x < -4.0) || (x > 4.0)  || (theta < -twelve_degrees) || (theta > twelve_degrees)) {
    if ((x < -4.0) || (x > 4.0)  || (theta < (-1*fourtyfive_degrees)) || (theta > fourtyfive_degrees)) {
//        r = -5000 - 50*fabs(x) - 100*fabs(theta);
        r = -10000 - 50*fabs(x) - 100*fabs(theta);      // max: -10000 -500 -20.9440 = -10521
//        r = -2000 - 25*fabs(x) - 50*abs(theta); // max -2260.5
        f = 1;
    } else {
//        r = 50 - 25*pow(abs(10*theta),2) - 10*abs(x) - 20*theta_dot;
//        r = 10 - 10*pow(fabs(10*theta),2) - 5*fabs(x) - 10*theta_dot; // max: -104.7884
        r =  10 - 1000*theta*theta - 5*fabs(x) - 10*theta_dot; // max: -104.7884
//        r = 5 - 2.5*pow(abs(theta),2) - 2.5*abs(x) - 2.5*theta_dot; // max: -28.3446
        f = 0;
    }

    DEBUG_MSG("reward: %f\n",r);

    reward->value = r;
    reward->success = f;
}

static inline void quantize_observations(struct frirl_desc *frirl, fri_float *states, int states_len, fri_float *new_states) {
    // perform quantization of the inputs, use only the allowed state values

    unsigned int sno, where;

    const double twelve_degrees = PI/15;
    const double three_degrees  = PI/60;
    double x[4];

    x[0] = states[0];
    x[1] = round(states[1]);
    x[2] = floor(states[2]/three_degrees)*three_degrees;
    x[3] = states[3];

    DEBUG_MSG("frirl_example_cartpole_quantize_observations: x[0]: %f x[1]: %f x[2]: %f x[3]: %f\n", x[0], x[1], x[2], x[3]);

    if (x[0] < 0) {
        x[0] = -1;
    }

    if (x[0] > 0) {
        x[0] = 1;
    }

    // x[0] == 0 -> zero stays zero

    if (x[1] < -1) {
        x[1] = -1;
    }

    if (x[1] > 1) {
        x[1] = 1;
    }

    if (x[2] > twelve_degrees) {
        x[2] = twelve_degrees;
    }

    if (x[2] < (-1 * twelve_degrees)) {
        x[2] = -1 * twelve_degrees;
    }

    if (x[3] < 0) {
        x[3] = -1;
    }

    if (x[3] > 0) {
        x[3] = 1;
    }

    new_states[0] = x[0];
    new_states[1] = x[1];
    new_states[2] = x[2];
    new_states[3] = x[3];
}


#ifdef BUILD_VISUALIZATION
static inline void draw(struct frirl_desc *frirl, fri_float *new_curr_state, double action, unsigned int steps) {
    // Ported from MATLAB to GLUT by Alex Toth <toth.alex@iit.uni-miskolc.hu>, 2018

    const double divider = 7; // MATLAB graph width is 12, here (opengl) only 2
    char text[32];

    double x = new_curr_state[0];
    double y = new_curr_state[1];
    double theta = new_curr_state[2];

    // Steps
    sprintf(text, "Episode: %d Steps: %d", frirl->episode_num, steps);

    startDrawing();

    drawTitle(COLOR_BLACK, text);

    drawAxis(-6, 6, divider, 2);

    // Car
    const double carColor[] = {.6, .6, .5};
    const double carX1 = (x - 1) / divider;
    const double carY1 = 0.25 / divider;
    const double carX2 = (x + 1) / divider;
    const double carY2 = 1.00 / divider;

    drawFilledRectangle(carColor, carX1, carY1, carX2, carY2);
    drawRectangle(COLOR_BLACK, carX1, carY1, carX2, carY2, 2);

    // Car Wheels
    const double wheelColor[] = {.5, .5, .5};
    const double wheelX1 = (x - .5) / divider;
    const double wheelY1 = 0.25 / divider;
    const double wheelX2 = (x + .5) / divider;
    const double wheelY2 = 0.25 / divider;
    const double wheelRadius = 0.40 / divider;

    drawFilledCircle(wheelColor, wheelX1, wheelY1, wheelRadius);
    drawCircle(COLOR_BLACK, wheelX1, wheelY1, wheelRadius, 2);
    drawFilledCircle(wheelColor, wheelX2, wheelY2, wheelRadius);
    drawCircle(COLOR_BLACK, wheelX2, wheelY2, wheelRadius, 2);

    // Pendulum
    const double pendulumX1 = x / divider;
    const double pendulumY1 = 1.00 / divider;
    const double pendulumX2 = (x + 1 * sin(theta)) / divider;
    const double pendulumY2 = (1.40 + 1.2 * cos(theta)) / divider;
    const double pendulumTopRadius = 0.30 / divider;

    drawLine(COLOR_RED, pendulumX1, pendulumY1, pendulumX2, pendulumY2, 5);
    drawPoint(COLOR_BLACK, pendulumX1, pendulumY1, 5);
    drawFilledCircle(COLOR_RED, pendulumX2, pendulumY2, pendulumTopRadius);
    drawCircle(COLOR_BLACK, pendulumX2, pendulumY2, pendulumTopRadius, 2);

    // Arrow
    double arrowFactorX = (action > 0 ? 1 : (action < 0 ? -1 : 0)) * 2.5;
    if (arrowFactorX > 0) {
        sprintf(text, "==>> %.0f", 10 * action);
    } else {
        if (arrowFactorX < 0) {
            sprintf(text, "%.0f <<==", 10 * action);
        } else {
            strcpy(text, "=0=");
            arrowFactorX = 0.25;
        }
    }

    drawText(COLOR_BLACK, (x + arrowFactorX - 0.5) / divider, 0.7 / divider, text);

    endDrawing();
}
#endif // BUILD_VISUALIZATION


int main(int argc, char **argv)
{
    struct frirl_desc frirl = frirl_desc_default;

    frirl_parse_cmdline(&frirl, argc, argv);

    frirl.reward_good_above = 0.0;

    frirl.alpha   = 0.3;
    frirl.gamma   = 1.0;
    frirl.epsilon = 0.001;

    frirl.qdiff_pos_boundary    = 1.0;
    frirl.qdiff_neg_boundary    = -200.0;
    frirl.qdiff_final_tolerance = 250.0;

    frirl.skip_rules = 1;  // Compatibility with the original MATLAB version

    frirl.do_action_func  = do_action;
    frirl.get_reward_func = get_reward;
    frirl.quant_obs_func  = quantize_observations;
#ifdef BUILD_VISUALIZATION
    frirl.draw_func       = draw;
#endif

    //statedim desc method 1
    struct frirl_dimension_desc statedim0 = {
        .values_len   = 2,
        .values_div   = 2.0,
        .values_steep = 1.0,
        .values_def   = 1.0,
        .universe_len = USIZE,
        .universe_div = 0.016
    };
    //alloc the .values member
//    FRIRL_ALLOC_VALUES(statedim0);
//    FRIRL_ALLOC_UNIVERSE(statedim0);
    //or a little simpler:
    FRIRL_ALLOC_DIM(statedim0);
    //set state and universe values, or in case of fixres:
//    FRIRL_GEN_FIXRES_VALUES(statedim0);
//    FRIRL_GEN_FIXRES_UNIVERSE(statedim0);
    //or a little simpler:
    FRIRL_GEN_FIXRES_DIM(statedim0);

    //statedim desc method 2
//    struct frirl_dimension_desc statedim0 = {
//        .values       = (frirl_float [2]) {0},
//        .values_len   = 2,
//        .values_div   = 2.0,
//        .values_steep = 1.0,
//        .values_def   = 1.0,
//        .universe     = (frirl_float [USIZE]) {0}, 
//        .universe_len = USIZE,
//        .universe_div = 0.016
//    };
//    FRIRL_GEN_FIXRES_VALUES(statedim0);
//    FRIRL_GEN_FIXRES_UNIVERSE(statedim0);

    //statedim desc method 3
//    struct frirl_dimension_desc statedim0 = {
//        .values       = (frirl_float []) {-1.0, 1.0 /*state values*/},
//        .values_len   = 2,
//        .values_div   = 2.0,
//        .values_steep = 1.0,
//        .values_def   = 1.0,
//        .universe     = (frirl_float []) { /*universe values*/ }, 
//        .universe_len = USIZE,
//        .universe_div = 0.016
//    };

    //statedim desc method 4
    //compact but not too readable
//    frirl_float s0[2]
//    frirl_float us0[USIZE];
//    struct frirl_dimension_desc statedim0 = {
//        s0, 2, 2.0, 1.0, 1.0,
//        us0, USIZE, 0.016
//    };
    //set state and universe values, or in case of fixres:
//    FRIRL_GEN_FIXRES_VALUES(statedim0);
//    FRIRL_GEN_FIXRES_UNIVERSE(statedim0);
    //or
//    FRIRL_GEN_FIXRES_DIM(statedim0);

    //statedim desc method 5
//    frirl_float s0[2]
//    frirl_float us0[USIZE];
//    struct frirl_dimension_desc statedim0;
//    statedim0.values       = s0;
//    statedim0.values_len   = 2;
//    statedim0.values_div   = 2.0;
//    statedim0.values_steep = 1.0;
//    statedim0.values_def   = 1.0;
//    statedim0.universe     = us0;
//    statedim0.universe_len = USIZE;
//    statedim0.universe_div = 0.016;
    //set state and universe values, or in case of fixres:
//    FRIRL_GEN_FIXRES_VALUES(statedim0);
//    FRIRL_GEN_FIXRES_UNIVERSE(statedim0);


    struct frirl_dimension_desc statedim1 = {
        .values_len   = 3,
        .values_div   = 1.0,
        .values_steep = 1.0,
        .values_def   = 0.0,
        .universe_len = USIZE,
        .universe_div = 0.032
    };
    FRIRL_ALLOC_DIM(statedim1);
    FRIRL_GEN_FIXRES_DIM(statedim1);

    // TODO: FRIRL_GEN_FIXRES_STATES gives wrong values
    struct frirl_dimension_desc statedim2 = {
        .values       = (fri_float []) {-0.2094, -0.1571, -0.1047, -0.0524, 0.0, 0.0524, 0.1047, 0.1571, 0.2094},
        .values_len   = 9,
        .values_div   = 0.0524,
        .values_steep = 21.485917317405871,
        .values_def   = 0.0,
        .universe_len = USIZE,
        .universe_div = 0.0031415926535897933
    };
    FRIRL_ALLOC_UNIVERSE(statedim2);
    FRIRL_GEN_FIXRES_UNIVERSE(statedim2);

    struct frirl_dimension_desc statedim3 = {
        .values_len   = 2,
        .values_div   = 2.0,
        .values_steep = 1.0,
        .values_def   = 0.0,
        .universe_len = USIZE,
        .universe_div = 0.016
    };
    FRIRL_ALLOC_DIM(statedim3);
    FRIRL_GEN_FIXRES_DIM(statedim3);

    frirl.statedims = (struct frirl_dimension_desc []) {statedim0, statedim1, statedim2, statedim3};
    frirl.statedims_len = 4;

    struct frirl_dimension_desc actiondim = {
        .values = (fri_float []) { -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.3999999999999999,
            -0.29999999999999992, -0.19999999999999995, -0.09999999999999998, 0.0,
            +0.09999999999999998, +0.19999999999999995, +0.29999999999999992,
            +0.3999999999999999, +0.5, +0.6, +0.7, +0.8, +0.9, +1.0 }, //(frirl_float [21]) {0},
        .values_len = 21,
        .values_div = 0.1,
        .universe_len = USIZE,
        .universe_div = 0.008
    };
//    FRIRL_ALLOC_VALUES(actiondim);
    FRIRL_ALLOC_UNIVERSE(actiondim);
//    FRIRL_ALLOC_DIM(actiondim);
//    FRIRL_GEN_FIXRES_VALUES(actiondim);
    FRIRL_GEN_FIXRES_UNIVERSE(actiondim);
//    FRIRL_GEN_FIXRES_DIM(actiondim);
    frirl.actiondim = actiondim;

    // Initialization
    frirl_init(&frirl);

    // do the hustle
    frirl_run(&frirl, 1);

    // display the constructed rule-base
//    printf("\nIncrementally constructed rulebase:\n");
//    frirl_show_rb(&frirl);
    // or with precise hexadecimal float
    //FRIRL_show_rb_hex(&frirl);

    // dump the RB to files
// TODO: timestamp
    frirl_save_rb_to_bin_file(&frirl,"cartpole.frirlrb.bin");
    frirl_save_rb_to_text_file(&frirl,"cartpole.frirlrb.txt");

    // Deinit
    frirl_deinit(&frirl);

//#ifdef BUILD_VISUALIZATION
//    // Deinit visualization
//    deinitGLUT();
//#endif

}

