      subroutine getpars()

      include 'common.h'

      npar = 0

100   format(a80)

10    continue
      npar = npar + 1
      read(5,100,end=20,err=40) parname(npar)
      goto 10

20    continue

      npar = npar - 1

      opened = 1

      return

 40   continue
      write(6,*) "getpars: error reading input file"
      end
c
c For compatibility with the ``Getpar'' library
c

      subroutine read_stdin()
      call getpars()
      return
      end
