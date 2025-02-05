/* pfplib program header */
/* AKG May, 1989 */
     
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pfplib.h"
     
int xargc;
char xargv[100][80];
FILE *infd,*outfd;
     
char datapath[160],outname[160];
char infilename[160],outfilename[160];
     
typedef struct { float re,im;} complex;
     
extern int getPars(int argc, char *argv[]);

/*
int main(int argc, char *argv[])
{
  getPars(argc, argv);
  return 1;
}
*/
int getPars(int argc, char *argv[])
{
  int i;
  FILE *datapathfile,*fopen();

  if(argc !=2)
  {
    printf("\nUsage: command  par.file\n\n");
    exit(1);
  }
  else
  {
    printf("\n---------------------------------------------\n");
    printf("Job command: %s  %s\n", argv[0], argv[1]);
    printf("---------------------------------------------\n");
  }

  xargc = argc;
     
  for (i=0;i<xargc;i++)
    strcpy(xargv[i],argv[i]);
     
  if(sfetch("in",infilename))
  {
    fprintf(stderr,"in=%s\n",infilename);
    if ((infd = fopen(infilename,"r")) == NULL)
    {
      fprintf(stderr,"pfplib: Unable to open infd %s\n",infilename);
      exit(0);
    }
  }
     
  if (isapipe(fileno(stdout)))
  {
    sprintf(outname,"/tmp/%d",getpid());
  }
  else
  {
    if (!sfetch("out",outfilename))
    {
      strcpy(outfilename,"/dev/null");
    }
     
    if (outfilename[0] != '/')
    {
      if ((datapathfile = fopen(".datapath","r")) == NULL)
      {
        fprintf(stderr,"pfplib: Unable to find .datapath\n");
        exit(0);
      }
      fscanf(datapathfile,"datapath=%s",datapath);
      sprintf(outname,"%s%s",datapath,outfilename);
    }
    else
    {
      strcpy(outname,outfilename);
    }
  }
  outfd = fopen(outname,"w");
  /*printf("in=%s\n",outname);*/
     
  return 1;
}

