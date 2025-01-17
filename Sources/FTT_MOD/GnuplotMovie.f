      PROGRAM GnuplotMovie 
      implicit none       
      integer i,istart,istep,imax,numb
      character*3 suffix,itoa    

      write(6,*) 'read istart,istep,imax'
      read*,istart,istep,imax
      print*,istart,istep,imax

      OPEN(unit=18,file='gmv',status="unknown")  
      do i=istart,imax,istep
         numb = i
         suffix = itoa(numb)
         write(18,*) 
     *'plot [0:1.5][0:1]' //char(39)//'mesh.'//suffix//char(39)//'w l'

      enddo 

      stop
      end


      character*3 function itoa(i)
      implicit none
      INTEGER i, nmb0,rest
 
      nmb0 = ichar('0')
      rest = mod(i,10)
      i = i/10
      itoa(3:3) = char(rest+nmb0)
      rest = mod(i,10)
      i = i/10
      itoa(2:2) = char(rest+nmb0)
      rest = mod(i,10)
      i = i/10
      itoa(1:1) = char(rest+nmb0)
 
      return
      end
                             
