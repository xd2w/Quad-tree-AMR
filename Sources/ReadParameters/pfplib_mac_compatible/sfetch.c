#include <stdio.h>
#include <string.h>
#include "pfplib.h"

extern int opened;
extern int npar;
extern char *parname[100],*parvalue[100];
extern char *string,*parfile;

int sfetch(name,var)

char *name;
char *var;

{
  int i,len;
   
  if (!opened)
    getpars();

  for(i=npar-1;(i>=0)&&(strcmp(name,parname[i])!=0);i--);

  if (i>=0)
  {
    if (strlen(parvalue[i])>0)
      sscanf(parvalue[i],"%s",var);

    if (var[0] == '"')
    {
      len = strlen(var);

      for (i=0;i<len-2;i++)
        var[i] = var[i+1];
 
      var[len-2] = '\0';
    }
    return(1);
  }
  return(0);
}
