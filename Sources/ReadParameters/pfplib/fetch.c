/* Modified 17/10/90 by S. Zaleski. Returns something now */

extern int opened;
extern int npar;
extern char *parname[100],*parvalue[100];
extern char *string,*parfile;

int fetch(name,vartype,var)

char *vartype;
char *name;
int *var;

{
  if (vartype[0] == 'd')
    return(ifetch(name,var));
  else if (vartype[0] == 'f')
    return(ffetch(name,var));
  else if (vartype[0] == 's')
    return(sfetch(name,var));
}
