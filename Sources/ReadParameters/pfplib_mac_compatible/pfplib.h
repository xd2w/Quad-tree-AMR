#ifndef  PFPLIB_H
#define  PFPLIB_H

extern int getPars(int argc, char *argv[]);
extern void getpars(void);
extern int ifetch(char *name, int *var);
extern int dfetch(char *name, double *var);
extern int sfetch(char *name, char *var);
extern int isapipe(int fd);

#endif


