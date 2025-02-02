#ifndef _LINALG_H_
#define _LINALG_H_

#include "nrutil.h"

extern void gaussj(float **a, int n, float **b, int m);
extern void svdcmp(float **a, int m, int n, float w[], float **v);
extern void svbksb(float **u, float w[], float **v, int m, int n, float b[], float x[]);

#endif