c     $Author: raus $
c     $Revision: 1.2 $
c     $Log: GP_get.f,v $
c Revision 1.2  1993/04/21  16:55:20  raus
c string size up to 128 chars
c
c Revision 1.1  1993/04/20  17:44:51  raus
c 1st save
c
c***********************************************************************
      CHARACTER*128 FUNCTION GETSTRING(i_line)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      INTEGER i_line
      INTEGER ibegin, i, iend
      CHARACTER*128 line1,line2
      CHARACTER*1 eq_char, blank
      DATA eq_char/'='/, blank/' '/
      INTRINSIC INDEX
c***
      line1 = lines_of_text(i_line)
      ibegin = 1 + INDEX(line1,eq_char)
      do i = 1, 128
         line2(i:i) = blank
      enddo
      line2(1:128 - ibegin + 1) = line1(ibegin:128)
      read(line2,'(A)') GETSTRING
c***
c---      write(*,*) line1(1:ibegin-1), '>>>', GETSTRING, '<<<'
c***
      return 
      end
c***********************************************************************
      DOUBLE PRECISION FUNCTION GETDOUBLE(i_line)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      INTEGER i_line, ibegin, i
      CHARACTER*128 line1,line2
      CHARACTER*1 eq_char, blank
      DATA eq_char/'='/, blank/' '/
      INTRINSIC INDEX
c***
      line1 = lines_of_text(i_line)
      ibegin = 1 + INDEX(line1,eq_char)
      do i = 1, 128
         line2(i:i) = blank
      enddo
      line2(ibegin:128) = line1(ibegin:128)
      read(line2(ibegin:128),'(BN,D20.15)') GETDOUBLE
c      read(line2,*) GETDOUBLE
c***
c      write(*,*) line1(1:ibegin-1), GETDOUBLE
c      write(*,*) line1(1:ibegin-1), line2(ibegin:128) 
c***
      return 
      end
c***********************************************************************
      REAL FUNCTION GETREAL(i_line)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      INTEGER i_line, ibegin, i
      CHARACTER*128 line1,line2
      CHARACTER*1 eq_char, blank
      DATA eq_char/'='/, blank/' '/
      INTRINSIC INDEX, len
c***
      line1 = lines_of_text(i_line)
      ibegin = 1 + INDEX(line1,eq_char)
      do i = 1, 128
         line2(i:i) = blank
      enddo
      line2(ibegin:128) = line1(ibegin:128)
      read(line2(ibegin:128),'(BN,F15.9)') GETREAL
c      read(line2,*) GETREAL
c***
c---      write(*,*) line1(1:ibegin-1), GETREAL
c***
      return 
      end
c***********************************************************************
      INTEGER FUNCTION GETINTEGER(i_line)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      INTEGER i_line, ibegin, i
      CHARACTER*128 line1,line2
      CHARACTER*1 eq_char, blank
      DATA eq_char/'='/, blank/' '/
      INTRINSIC INDEX
c***
      line1 = lines_of_text(i_line)
      ibegin = 1 + INDEX(line1,eq_char)
      do i = 1, 128
         line2(i:i) = blank
      enddo
      line2(ibegin:128) = line1(ibegin:128)
c      read(line2,*) GETINTEGER
      read(line2(ibegin:128),'(I8)') GETINTEGER
c***
c      write(*,*) line1(1:ibegin-1), GETINTEGER
c      write(*,*) line1(1:ibegin-1), line2(ibegin:128) 
c***
      return 
      end

