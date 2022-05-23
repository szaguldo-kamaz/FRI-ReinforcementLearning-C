#ifndef DRAW_H
#define DRAW_H

#include <stdio.h>
#include <GL/freeglut.h>
#include <math.h>

void startDrawing();
void endDrawing();

void drawPoint(const double *color, double x, double y, int pointSize);
void drawLine(const double *color, double x0, double y0, double x1, double y1, int lineWidth);
void drawRectangle(const double *color, double x0, double y0, double x1, double y1, int lineWidth);
void drawFilledRectangle(const double *color, double x0, double y0, double x1, double y1);
void drawCircle(const double *color, double cx, double cy, double radius, int lineWidth);
void drawFilledCircle(const double *color, double cx, double cy, double radius);
void drawTriangle(const double *color, double x0, double y0, double x1, double y1, double x2, double y2, int lineWidth);
void drawFilledTriangle(const double *color, double x0, double y0, double x1, double y1, double x2, double y2);
void drawAsterisk(const double *color, double x, double y, int branches, int lineWidth, double size);
void drawFilledAsterisk(const double *color, double x, double y, int branches, double size);
void drawText(const double *color, double x, double y, char *text);
void drawTitle(const double *color, char *title);
void drawAxis(int min, int max, double divider, int steps);

#endif
