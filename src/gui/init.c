#include "init.h"
#include <stdio.h>


void display () {}


void key_handler(unsigned char key, int x, int y) {

    if (key == 27) {
        deinitGLUT();
    }

}


void reshapeHandler(int width, int height) {

    glViewport(0, 0, width, height);
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    gluPerspective(0, (double)width / (double)height, 0.1, 100.0);

}


void key_up_handler(unsigned char key, int x, int y) { }


void initGLUT(int *argc, char **argv, int width, int height, const char *title) {

    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowSize(width, height);
    glutCreateWindow(title);
    glutDisplayFunc(display);

//    glViewport(0, 0, width, height);
    glutKeyboardFunc(key_handler);
//        glutKeyboardUpFunc(key_up_handler);
//    glutSpecialFunc(SpecialInput);
    glutReshapeFunc(reshapeHandler);

    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    // glEnable(GL_BLEND);

    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_FALSE);

    glClearColor(1, 1, 1, 1);

}


void deinitGLUT() {

    glutDestroyWindow(glutGetWindow());
    glutExit();

}
