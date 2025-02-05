/*
 *	error abortion subroutine
 *	a is a printf format string, while b-h are optional arguments
 *	When the subroutine is called, the program will stop after
 *	printing its message.
 *	Example:	err("Cannot divide %f by %f\n", x, y);
 *
 *      Also has Fortran companion routines err_(string) and erexit_(string)
 *
 * Revised  4/23/84    stew     echo to stdout (header) as well
 */
#include <stdio.h>

/*VARARGS1*/
void err( a,b,c,d,e,f,g,h )
char *a,*b,*c,*d,*e,*f,*g,*h; {
	extern char **xargv;
/*	fprintf( stderr, "%s: ", xargv[0] );   */
	fprintf( stderr, a,b,c,d,e,f,g,h );
	fflush(stderr);
/*	fprintf( stdout, "%s: ", xargv[0] );    */
	fprintf( stdout, a,b,c,d,e,f,g,h );
	fflush(stdout);
	exit(-1);
}
