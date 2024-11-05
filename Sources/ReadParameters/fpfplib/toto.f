      program toto
      implicit none
      integer i
      character*8 parname

      parname = 'string'
      i=1
      print*, char(i+ichar('0'))
      print*, ichar('a')
      print*, ichar('z')
      print*, ichar('0')
      print*, ichar('9')
      print*,'intf.'//char(i+ichar('0'))
     *       //char(i+ichar('0'))//char(i+ichar('0'))
      print*,parname(1:7) // 'treatment' , len(parname)
      print*,parname(1:len(parname)+1) // 'treatment'

      end
