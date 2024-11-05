#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "variable.h"

/* ------------------------------------------------------------------- */
/* mx and mz define the normal direction; mx, mz and alpha,Real
  define the equation of the line. b the width of bloc */
Real VOL2(Real mx,Real mz,Real alpha,Real b)

{
  Real area, mb, alphaa;

  /* 0<mx or 0<mz */
  mb = mx*b;
  if (alpha <= 0.0 || b <= 0.0) 
    area = 0.0;
  else if (alpha >= (mb+mz))
    area = b;
  else {
    if (mb <= mz) {
      if (alpha < mb) 
	area = 0.5*alpha*alpha/(mx*mz);
      else if (alpha > mz) {
	alphaa = mb + mz - alpha;
	area = b - 0.5*alphaa*alphaa/(mx*mz);
      }
      else
	area = b*(alpha-0.5*mb)/mz;
    }
    else{
      if (alpha < mz) 
               area = 0.5*alpha*alpha/(mx*mz);
      else if (alpha > mb) {
	alphaa = mb + mz - alpha;
	area = b - 0.5*alphaa*alphaa/(mx*mz);
      }
      else
	area = (alpha-0.5*mz)/mx;
    }
  }

  return area;
}

