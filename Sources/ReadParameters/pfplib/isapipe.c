#include <stdio.h>
#include <errno.h>
extern int errno;

int
isapipe(fd)
int fd;
{
  extern long lseek();
  register long rc;

  rc = lseek(fd,0L,1);
/*
  fprintf(stderr,"isapipe: fd=%d, lseek rc=%d, errno=%d\n",fd,rc,errno);
*/
  if(rc == -1)
    if(errno == ESPIPE || errno == EINVAL)
      return(1);

  return(0);
}
