#include <stdio.h>
#include <string.h>
#include "pfplib.h"

extern int opened;
extern int npar;
extern char *parname[100],*parvalue[100];
extern char *string,*parfile;

int ifetch(name,var)

char *name;
int *var;

{
  int i;
   
  if (!opened)
    getpars();

  for(i=npar-1;(i>=0)&&(strcmp(name,parname[i])!=0);i--);

  if (i>=0)
  {
    sscanf(parvalue[i],"%d",var);
    return(1);
  }
  return(0);
}
