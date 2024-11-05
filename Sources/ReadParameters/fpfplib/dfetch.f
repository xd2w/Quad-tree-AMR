      integer function dfetch(name,var)

      include 'common.h'

      character *(*) name

      integer i,j,k,ffetch
      doubleprecision var,atod

      if (opened.eq.0) then
        call getpars()
      endif

      dfetch = 0

      do i=1,npar
        j = index(parname(i),name)
        if (j.ne.0) then
          k = len(name)

          if (((j.eq.1).or.(parname(i)(j-1:j-1).eq.' '))
     $        .and.(parname(i)(k+j:k+j).eq.'=')) then
            var = atod(parname(i)(k+j+1:len(parname(i))))
            dfetch = 1
          else if(parname(i)(j+k:j+k).eq.' ') then
             write(6,1000)
             write(6,*) parname(i)
             write(6,6666)
             stop 'dfetch killed you'
          endif

        endif
      enddo

      return

1000  format('dfetch: You are not allowed a blank between the',/
     $   'variable`s name and the = sign. Here is the offending line:')
6666  format('I will not let you continue in your silly ways.')

      end


      doubleprecision function atodsimp(string)

      include 'undefine.h'

      character *(*) string
      integer sign,first,c,c0,c9
      integer i,afterdec,funny

      doubleprecision postdec,predec,mul

      mul = 1.0
      sign = 1
      first = 1
      postdec = 0.0
      predec = 0.0
      afterdec = 0
      funny = 0

      c0 = ichar('0')
      c9 = ichar('9')

      if (string(1:1).eq.'-') then
        sign = -1
        first = 2
      endif

      do i=first,len(string)
        c = ichar(string(i:i))

        if (c.eq.ichar('.')) then
          afterdec = 1
        else if ((c.ge.c0).and.(c.le.c9)) then
          if ((afterdec.eq.0).and.(funny.eq.0)) then
            predec = 10.0d0*predec + 1.0d0*(c - c0)
          else if (funny.eq.0) then
            mul = mul * 0.1d0
            postdec = postdec + 1.0d0*(c - c0)*mul
          endif
        else
          funny = 1
        endif

      end do

      atodsimp = sign * (postdec+predec)

      return
      end



      doubleprecision function atod(string)

      include 'undefine.h'

      character *(*) string

      doubleprecision mantissa,atodsimp

      integer atoi,exponent
      integer hasexp

      hasexp = index(string,"e")
      if(hasexp.eq.0) then
        hasexp = index(string,"E")
      endif
      if(hasexp.eq.0) then
        hasexp = index(string,"d")
      endif
      if(hasexp.eq.0) then
        hasexp = index(string,"D")
      endif

      if(hasexp.eq.0) then
         atod = atodsimp(string)
      else
         mantissa = atodsimp(string(1:hasexp-1))
         exponent = atoi(string(hasexp+1:))
         atod     = mantissa*10.d0**(1.0d0*exponent)
      endif

      return
      end

