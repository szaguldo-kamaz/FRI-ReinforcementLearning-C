#include "draw.h"


void startDrawing() {

    glClear(GL_COLOR_BUFFER_BIT);
    glClearColor(1, 1, 1, 1);

}


void endDrawing() {

    glutSwapBuffers();
    glutMainLoopEvent();

}


void drawPoint(const double *color, double x, double y, int pointSize) {

    glPointSize(pointSize);

    glColor3dv(color);

    glBegin(GL_POINTS);
    glVertex2d(x, y);
    glEnd();

}


void drawLine(const double *color, double x0, double y0, double x1, double y1, int lineWidth) {

    glLineWidth(lineWidth);

    glColor3dv(color);

    glBegin(GL_LINES);
    glVertex2d(x0, y0);
    glVertex2d(x1, y1);
    glEnd();

}


void drawRectangle(const double *color, double x0, double y0, double x1, double y1, int lineWidth) {

    glLineWidth(lineWidth);

    glColor3dv(color);

    glBegin(GL_LINES);
    glVertex2d(x0, y0);
    glVertex2d(x1, y0);
    glVertex2d(x1, y0);
    glVertex2d(x1, y1);
    glVertex2d(x1, y1);
    glVertex2d(x0, y1);
    glVertex2d(x0, y1);
    glVertex2d(x0, y0);
    glEnd();

}


void drawFilledRectangle(const double *color, double x0, double y0, double x1, double y1) {

    glColor3dv(color);
    glRectd(x0, y0, x1, y1);

}


void drawCircle(const double *color, double cx, double cy, double radius, int lineWidth) {

    int i;
    const int segments = 1000;

    glLineWidth(lineWidth);

    glColor3dv(color);

    glBegin(GL_LINE_LOOP);

    for (i = 0; i < segments; i++) {

        double theta = 2.0 * 3.1415926 * i / segments;

        double x = radius * cos(theta);
        double y = radius * sin(theta);

        glVertex2d(x + cx, y + cy);

    }

    glEnd();

}


void drawFilledCircle(const double *color, double cx, double cy, double radius) {

    int i;
    const int segments = 1000;

    glColor3dv(color);

    glBegin(GL_POLYGON);

    for(i = 0; i < segments; i++) {

        double theta = 2.0 * 3.1415926 * i / segments;

        double x = radius * cos(theta);
        double y = radius * sin(theta);

        glVertex2d(x + cx, y + cy);

    }

    glEnd();

}


void drawTriangle(const double *color, double x0, double y0, double x1, double y1, double x2, double y2, int lineWidth) {

    glLineWidth(lineWidth);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glColor3dv(color);
    glBegin(GL_TRIANGLES);
    glVertex2d(x0, y0);
    glVertex2d(x1, y1);
    glVertex2d(x2, y2);
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

}


void drawFilledTriangle(const double *color, double x0, double y0, double x1, double y1, double x2, double y2) {

    glColor3dv(color);

    glBegin(GL_TRIANGLES);
    glVertex2d(x0, y0);
    glVertex2d(x1, y1);
    glVertex2d(x2, y2);
    glEnd();

}


void drawAsterisk(const double *color, double x, double y, int branches, int lineWidth, double size) {
    const double degree = 360.0 / branches * .5;
    const double PI = 3.14159265359;
    const double PIper180 = PI / 180.0;
    double tmp;
    int i, j;

    glLineWidth(lineWidth);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glColor3dv(color);

    glBegin(GL_POLYGON);

    for (i = 0, j = branches - 1; i < branches * 2; i++) {
        if (i % 2 == 1)
        {
            tmp = ((degree * i) + 90 + 180) * PIper180;
            glVertex2d(x + cos(tmp) * size, y + sin(tmp) * size);
        }
        else
        {
            tmp = ((degree * i) + 90 + 180) * PIper180;
            glVertex2d(x + cos(tmp) * size * .4, y + sin(tmp) * size * .4);
        }
    }

    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

}


void drawFilledAsterisk(const double *color, double x, double y, int branches, double size) {

    const double degree = 360.0 / branches * .5;
    const double PI = 3.14159265359;
    const double PIper180 = PI / 180.0;
    double tmp;
    int i, j;

    glColor3dv(color);

    glBegin(GL_POLYGON);

    for (i = 0, j = branches - 1; i < branches * 2; i++) {

        if (i % 2 == 1) {
            tmp = ((degree * i) + 90 + 180) * PIper180;
            glVertex2d(x + cos(tmp) * size, y + sin(tmp) * size);
        } else {
            tmp = ((degree * i) + 90 + 180) * PIper180;
            glVertex2d(x + cos(tmp) * size * .4, y + sin(tmp) * size * .4);
        }

    }

    glEnd();

}


void drawText(const double *color, double x, double y, char *text) {

    glColor3dv(color);

    glRasterPos2d(x, y);

    glutBitmapString(GLUT_BITMAP_HELVETICA_12, text);
}


void drawTitle(const double *color, char *title) {

    glColor3dv(color);

    glRasterPos2d(0.0, 0.8);

    glutBitmapString(GLUT_BITMAP_HELVETICA_18, title);

}


void drawAxis(int min, int max, double divider, int steps) {

    int i;

    const double color[] = {0, 0, 0};
    char num[10];

    const double mind = min / divider;
    const double maxd = max / divider;
    double id;


    drawLine(color, mind, 0, maxd, 0, 3); // x axis
    drawLine(color, mind, mind, mind, maxd, 3); // y axis

    for (i = min; i <= max; i += steps) {
        id = i / divider;

        drawLine(color, id, 0.02, id, -0.03, 3); // x axis
        drawLine(color, mind + 0.02, id, mind - 0.02, id, 3); // y axis

        sprintf(num, "%d", i);
        drawText(color, id, -0.1, num);
        drawText(color, mind - 0.1, id, num);
    }

}
