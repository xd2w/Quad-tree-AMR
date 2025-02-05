/* AKG attempt at mimicking the fetch function */
/* Modified 24/04/93 by SZ to initialize correctly parname and parvalue */


#include <stdio.h>
#include <string.h>

int opened=0;
char *parname[100],*parvalue[100];
char string[80],parfile[80];
int npar=0;

extern int xargc;
extern char xargv[100][80];
extern char *alloc(int size);

void getpars(void)

{
  FILE *fp,*fopen();
  
  int i,namelen,foundparfile=0;
  char *equal;

  int *varint;
  float *varfloat;
  
  fp = fopen(xargv[1],"r"); 
 
  while (fscanf(fp, "%s",string) != EOF)
  {
    parname[npar]  = (char *) alloc(80);
    parvalue[npar] = (char *) alloc(80);
    
    bzero(parname[npar],80);
    bzero(parvalue[npar],80);

    equal = (char *) strchr(string,'=');

    if (equal != NULL)
    {
      namelen = equal-string;
      strncpy(parname[npar],string,namelen);
      strcpy(parvalue[npar],string+namelen+1);
      npar++;
    }
  }

  opened = 1;
}
