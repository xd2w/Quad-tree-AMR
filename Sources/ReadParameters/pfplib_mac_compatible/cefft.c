/* Four fft subroutines based on the algorithm used in fork (FGDP, p. 12) 
 * and modeled after Clayton's and Ottolini's subroutines
 *********************************************************************
 *								     *
 *		       lx-1					     *
 *	x(k) = scale * sum  [x(j)*exp(isign*2*pi*sqrt(-1)*j*k/lx]    *
 *		       j=0					     *
 *								     *
 *********************************************************************
 *
 *		refft (x,lx,isign,scale,mode)
 *		cefft (x,lx,isign,scale)
 *		rvfft (x,lx,nx,isign,scale,mode)
 *		cvfft (x,lx,nx,isign,scale)
 *
 * See each subroutine below for a description of the arguments.
 * Dave Hale, 3/17/82
 */
#include "Cmplxlib.h"
static double sintab[]
	   ={
		1.000000000000000e+00, /* sin(pi/2) */
		7.071067811865475e-01,
		3.826834323650897e-01,
		1.950903220161282e-01,
		9.801714032956060e-02,
		4.906767432741801e-02,
		2.454122852291228e-02,
		1.227153828571992e-02,
		6.135884649154475e-03,
		3.067956762965976e-03,
		1.533980186284765e-03,
		7.669903187427045e-04,
		3.834951875713956e-04,
		1.917475973107033e-04,
		9.587379909597734e-05,
		4.793689960306688e-05,
		2.396844980841821e-05,
		1.198422490506970e-05,
		5.992112452642428e-06,
		2.996056226334660e-06,
		1.498028113169011e-06	/* sin(pi/(2**21)) */
	   };

/* 	cefft	lx complex to lx complex
 *
 *	x	input/output array
 *	lx	length of complex array; must be a power of 2
 *	isign	sign of sqrt(-1)
 *	scale	scale factor; sqrt(1./lx) conserves energy
*/
cefft (x,lx,isign,scale)
complex *x;
int lx,isign; float scale;
   {
	register complex *px,*qx;
	complex *xplx;
	int m,j,k,step;
	float temp,real,imag;
	double cn,sn,cd,sd,dtemp,*psintab;

	/* bit reverse */
	xplx = x+lx;
	for(px=x, j=0; px<xplx; px++, j+=m)
	   {
		if(px < (qx=x+j))
		  {
			temp = qx->re; qx->re = px->re; px->re = temp;
			temp = qx->im; qx->im = px->im; px->im = temp;
		  }
		for (m=lx>>1; m>=1 && j>=m; j-=m, m>>=1);
	   }
	/* first butterfly and scaling */
	for(px=x, qx=x+1; px<xplx; px+=2, qx+=2)
	   {
		if (scale != 1.)
		  {
			px->re *= scale; px->im *= scale;
			qx->re *= scale; qx->im *= scale;
		  }
		temp = qx->re; qx->re = px->re-temp; px->re += temp;
		temp = qx->im; qx->im = px->im-temp; px->im += temp;
	   }
	/* remaining butterflies */
	for (j=2, psintab=sintab; j<lx; j=step)
	   {
		step = j<<1;
		sd = *psintab++;
		if (isign < 0) sd = -sd;
		dtemp = *psintab;
		cd = 2.*dtemp*dtemp;
		sn = 0.;
		cn = 1.;
		for(k=0; k<j; k++)
		   {
			for(px=x+k; px<xplx; px+=step)
			   {
				qx = px+j;
				real = cn*qx->re-sn*qx->im;
				imag = sn*qx->re+cn*qx->im;
				qx->re = px->re-real;
				qx->im = px->im-imag;
				px->re += real;
				px->im += imag;
			   }
			dtemp = cd*cn+sd*sn;
			sn += sd*cn-cd*sn;
			cn -= dtemp;
		   }
	   }
   }
