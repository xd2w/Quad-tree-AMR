/* AKG attempt at mimicking the snap function */

#include <stdio.h>

char *fopened[10];
int nopened=0;

void snap(filename,fn1,fn2,fn3,nbytes,ptr)

char *filename;
int fn1,fn2,fn3,nbytes;
char *ptr;

{
  FILE *hp,*dp,*fopen();
  char hfilename[40],dfilename[40];

  int found,i;

  sprintf(hfilename,"%s.h",filename);

  hp = fopen(hfilename,"w");

  fprintf(hp,"n1=%d\n",fn1);
  fprintf(hp,"n2=%d\n",fn2);
  fprintf(hp,"n3=%d\n",fn3);
/*  fprintf(hp,"sets next in=\"%s.h@\"\n",filename);*/

  fflush(hp);
  fclose(hp);

  sprintf(dfilename,"%s.h@",filename);

  found = 0;

  for (i=0;i<nopened;i++)
  {
    if (strcmp(filename,fopened[i]) == 0)
    {
      dp = fopen(dfilename,"a");
      found = 1;
    }
  }

  if (!found)
  {
    fopened[nopened] = (char *) alloc(40);
    strcpy(fopened[nopened],filename);
    nopened++;
    
    dp = fopen(dfilename,"w");
  }

  rite(dp,ptr,nbytes);

  fflush(dp);
  fclose(dp);
}
