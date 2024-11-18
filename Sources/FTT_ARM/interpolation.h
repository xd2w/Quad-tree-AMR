#ifndef INTERPOLATION_H
#define INTERPOLATION_H

// extern void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
// extern void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);

extern void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
extern void splint(double xa[], double ya[], double y2a[], int n, float x, double *y);

#endif