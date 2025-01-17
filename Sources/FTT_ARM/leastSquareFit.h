#ifndef _LEASTSQUAREFIT_H_
#define _LEASTSQUAREFIT_H_

#include "nrutil.h"

extern void lfit(float x[], float y[], float sig[], int ndat, float a[], int ia[],
                 int ma, float **covar, float *chisq, void (*funcs)(float, float[], int));

extern void covsrt(float **covar, int ma, int ia[], int mfit);

extern void svdfit(float x[], float y[], float sig[], int ndata, float a[], int ma,
                   float **u, float **v, float w[], float *chisq,
                   void (*funcs)(float, float[], int));

extern void svdvar(float **v, int ma, float w[], float **cvm);

#endif
