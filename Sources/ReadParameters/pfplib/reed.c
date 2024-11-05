#include <stdio.h>

#define MIN(A,B) ( ((A) < (B)) ? (A) : (B) )

reed (datafile,ptr,length)
FILE *datafile;
int length;
char *ptr;
/*
 *        buffered read with error detection
 *
 * modified 2/22/83 - added return value
 * modified 3/6/83  - removed infinite loop on short data
 * modified 4/1/83  - added file descriptor to error messages - stew
 * modified 7/14/83 - added perror and err diagnostics - stew
 * modified 8/3/83 - increased BUFSIZ to 65536 if not a pipe - stew
 */
{
  extern int read();
  int nread, total, filedes;

  static int bufsiz = 65536;

  filedes = fileno(datafile);

  for (total=0; (total < length) && 
    (0 != (nread = read(filedes,ptr+total,MIN(length-total,bufsiz))));
    total += nread)
  {
    if (nread<0)
    {
      perror("reed()");
    }
  }
  return(total);
}
