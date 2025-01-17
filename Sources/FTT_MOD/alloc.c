/*
 *	allocation with error detection
 *
 * modified 1/26/83  S. Levin : added unsigned coersion to agree with malloc(3)
 *                              changed to NULL instead of 0
 * modified 3/19/83  S. Levin : keep internal 1024 byte block for small request
 *				include source for malloc inline to block missing
 *				cfree() linkage references.
 * modified 2/7/84   S. Levin : deleted calloc() call to reduce paging -- no
 *				longer is memory necessarily zeroed out of alloc.
 *
 * In order to allocate an array of floating point numbers, use
 * the following command in the calling routine:
 *
 *	float *x;
 *	x = (float *) alloc(nx*sizeof(float));
 *
 * nx is the number of elements needed in the array. 
 */
#include <stdio.h>
#define SMALLBLOCK -1
char *alloc (size)
int size;
{
	static char *myblock = NULL; static int bytesleft = 0;
	char *ptr; extern char *malloc();

	if(size > SMALLBLOCK)
	   {
	    if ((ptr = malloc ((unsigned)size)) == NULL)
		fprintf(stderr,"can't allocate %d bytes\n",size);
	    return (ptr);
	   }
	else
	  { 
	    if(bytesleft < size)
	       {
	        if ((myblock = malloc ((unsigned)SMALLBLOCK)) == NULL)
		   fprintf(stderr,"can't allocate %d bytes\n",size);
		bytesleft = SMALLBLOCK - size;
	       }
	   else
	       {
		bytesleft -= size;
	       }
	   ptr = myblock;
	   myblock += size;
	   return(ptr);
	  }
}
