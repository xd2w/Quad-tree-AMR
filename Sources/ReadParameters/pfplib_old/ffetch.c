/* AKG attempt at mimicking the fetch function */

#include <stdio.h>

extern int opened;
extern int npar;
extern char *parname[100],*parvalue[100];
extern char *string,*parfile;

int ffetch(name,var)

char *name;
float *var;

{
  int i;

  if (!opened)
    getpars();

  for(i=npar-1;(i>=0)&&(strcmp(name,parname[i])!=0);i--);

  if (i>=0)
  {
    sscanf(parvalue[i],"%f",var);
    return(1);
  }
  return(0);
}
