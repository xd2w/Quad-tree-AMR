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

/*
 * 	refft	lx real to lx/2 (or lx/2+1) complex or vice-versa
 *
 *	x	input/output array
 *	lx	length of real array; must be a power of 2
 *	isign	sign of sqrt(-1)
 *	scale	scale factor; sqrt(2./lx) conserves energy
 *	mode	 1 for lx real to lx/2 complex
 *		 2 for lx real to lx/2+1 complex
 *		-1 for lx/2 complex to lx real
 *		-2 for lx/2+1 complex to lx real
 */
refft(x,lx,isign,scale,mode)
register complex *x;
int lx,isign,mode; float scale;
   {
	register complex *xp,*xn;
	float real,imag,xsumre,xsumim,xdifre,xdifim;
	double aa,cn,sn,cd,sd,arg,sin();

	if(mode > 0) 				/* real to complex */
	  {
		cefft(x,lx/2,isign,scale);
		real = x->re+x->im; imag = x->re-x->im;
		x->re = real;
		if (mode == 1) 			/* pack Nyquist */
			x->im = imag;
		else				/* unpack Nyquist */
		  {
			xn = x+lx/2;
			xn->re = imag;
			xn->im = 0.; 
			x->im = 0.;
		  }
	  }
	cn = 1.; sn = 0.;			/* initial cosine and sine */
	arg = pi/lx;				/* = pi/lx */
	if (isign < 0) arg = -arg;
	aa = sin(arg);
	cd = 2.*aa*aa; sd = sin(arg+arg);	/* for cosine/sine recursion */
	for (xp = x+1, xn = x+lx/2-1; xp <= xn; xp++, xn--)
	  {
		aa = cd*cn+sd*sn;
		sn += sd*cn-cd*sn;		/* sin(2*arg*(xp-x)) */
		cn -= aa;			/* cos(2*arg*(xp-x)) */
		xsumre = 0.5*(xp->re+xn->re);
		xsumim = 0.5*(xp->im-xn->im);
		xdifre	= 0.5*(xp->re-xn->re);
		xdifim = 0.5*(xp->im+xn->im);
		real = sn*xdifre+cn*xdifim;
		imag = sn*xdifim-cn*xdifre;
		xp->re = xsumre+real;
		xp->im = imag+xsumim;
		xn->re = xsumre-real;
		xn->im = imag-xsumim;
	  }
	if(mode < 0) 				/* complex to real */
	  {
		if (mode == -2)			/* Nyquist not packed, */
		  {				/* so pack it */
			xn = x+lx/2;
			x->im = xn->re;
			xn->re = 0.;
		  }
		real = 0.5*(x->re+x->im);
		x->im = 0.5*(x->im-x->re);
		x->re = real;
		cefft(x,lx/2,isign,scale);
		for(xp=x,xn=x+lx/2; xp<xn; xp++)
			xp->im = -xp->im;
	   }
   }
