/* AKG attempt at mimicking the fetch function */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pfplib.h"

extern int opened;
extern int npar;
extern char *parname[100],*parvalue[100];
extern char *string,*parfile;

int dfetch(name,var)

char *name;
double *var;

{
  int i;
  float a;

  if (!opened)
    getpars();

  for(i=npar-1;(i>=0)&&(strcmp(name,parname[i])!=0);i--);

  if (i>=0)
  {
    sscanf(parvalue[i],"%f",&a);
    *var = a;
    return 1;
  }
  return 0;
}
