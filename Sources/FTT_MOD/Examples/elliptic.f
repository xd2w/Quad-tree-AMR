      program elliptic
      integer i, n
      real xc, yc, a, b, theta

      print*,'draw an elliptic'
      print*,'read n '
      read*, n
      print*,'read xc, yx'
      read*, xc, yc
      print*,'read a, b'
      read*, a, b

      theta = 2.*acos(-1.0)/n
       write(18,*) n
      do i=1,n
         angle=i*theta
         write(18,*) xc+a*cos(angle), yc+b*sin(angle)
      enddo

      end
