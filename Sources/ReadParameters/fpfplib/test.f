      program test

      integer a,x
      real b
      doubleprecision c

      a = 10
      b = 10.0
      c = 10.0

      x = ifetch('andrew',a)
      x = ffetch('david',b)
      x = dfetch('john',c)

      write (6,*) ' a=',a
      write (6,*) ' b=',b
      write (6,*) ' c=',c

      stop
      end
