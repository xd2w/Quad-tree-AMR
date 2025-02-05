#include <stdio.h>

#define MIN(A,B) ( ((A) < (B)) ? (A) : (B) )

rite (datafile,ptr,length)
FILE *datafile;
int length;
char *ptr;
/*
 *	buffered write with error detection
 *
 * modified 4/1/83  - added file descriptor to error messages - stew
 */
{
	extern int isapipe();
	int filedes,nwrite, total;
	int bufsiz = 65536;

        filedes = fileno(datafile);

	if(isapipe(filedes)) bufsiz = 4096;

	for (total=0; (total < length) && 
		(0 != (nwrite = write(filedes,ptr+total,MIN(length-total,bufsiz))));
		total += nwrite)
		{
		 if (nwrite<0)
		 {
		  perror("rite()");
		 }
		}
	return(total);
}



