      integer function ifetch(name,var)

      include 'common.h'

      character *(*) name

      integer var,atoi,i,j,k

      if (opened.eq.0) then
        call getpars()
      endif

      ifetch = 0

      do i=1,npar
        j = index(parname(i),name)
        if (j.ne.0) then
          k = len(name)
          if (((j.eq.1).or.(parname(i)(j-1:j-1).eq.' '))
     $        .and.(parname(i)(k+j:k+j).eq.'=')) then
            var = atoi(parname(i)(k+j+1:len(parname(i))))
            ifetch = 1
          else if(parname(i)(j+k:j+k).eq.' ') then
             write(6,1000)
             write(6,*) parname(i)
             write(6,6666)
             stop 'ifetch killed you'
          endif
        endif
      enddo

      return

1000  format('ifetch: You are not allowed a blank between the',/
     $   'variable`s name and the = sign. Here is the offending line:')
6666  format('I will not let you continue in your silly ways.')

      end


      integer function atoi(string)

      include 'undefine.h'

      character *(*) string
      integer sign,first,var,c,c0,c9
      integer funny,i

      sign = 1
      first = 1
      var = 0
      funny = 0

      c0 = ichar('0')
      c9 = ichar('9')

      if (string(1:1).eq.'-') then
        sign = -1
        first = 2
      endif

      do i=first,len(string)
        c = ichar(string(i:i))
        if ((c.ge.c0).and.(c.le.c9).and.(funny.eq.0)) then
          var = 10*var + c - c0
        else
          funny = 1
        endif
      end do

      atoi = sign * var

      return
      end
