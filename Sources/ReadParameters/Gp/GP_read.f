c     $Author: raus $
c     $Revision: 1.3 $
c     $Log: GP_read.f,v $
c Revision 1.3  1993/07/20  15:03:20  raus
c fixed a bug: lines with ' ' after '=' are removed
c
c Revision 1.2  1993/04/21  16:56:30  raus
c string size up to 128 chars
c
c Revision 1.1  1993/04/20  17:46:13  raus
c 1st save
c
c***********************************************************************
      SUBROUTINE read_stdin
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      CHARACTER*128 line_of_text, line_aux
      CHARACTER*1 blank, eq_char, tab, com_1, com_2, char
      DATA blank/' '/, eq_char/'='/, tab/'\t'/, com_1/'#'/, com_2/'!'/
      INTEGER i, k, m, count, len_line
      INTRINSIC LEN, INDEX 
c***  
      count = 0
      do k = 1,nlines_max
         read (*,'(A128)',ERR=100,END=100) line_of_text 
         len_line = LEN(line_of_text)
         m = INDEX(line_of_text,eq_char)
c***  
c     lines with no '=' are discarded
c***  
         if (m .gt. 0) then
c***  
c     the text line is copied into the initialised auxiliary line,
c     all blanks and tabs are removed in the process
c***  
            do i=1,len_line
               line_aux(i:i) = blank
            enddo 
            m = 1
            do i = 1,len_line
               char = line_of_text(i:i)
               if (char .ne. blank .and. char .ne. tab) then
                  line_aux(m:m) = char
                  m = m + 1
               endif
            enddo 
c***  
c     discard line if 1st char is a '=', '#', '!' or ' ' and also
c     if 1st character after '=' is a '#', '!' or ' '
c***  
            m = 1
            char = line_aux(1:1)
            if (char .eq. eq_char .or. char .eq. com_1 .or.
     %          char .eq. com_2 .or. char .eq. blank) m = 0 
            i = INDEX(line_aux,eq_char)
            char = line_aux(i+1:i+1)   
            if (char .eq. com_1 .or. char .eq. com_2 .or.
     %          char .eq. blank) m = 0
c***  
c     eventual comments at end of line are discarded
c***  
            if (m .gt. 0) then
               m = INDEX(line_aux,com_1)
               if (m .gt. 0) then
                  do i = m,len_line
                     line_aux(i:i) = blank
                  enddo 
               endif 
               m = INDEX(line_aux,com_2)
               if (m .gt. 0) then
                  do i = m,len_line
                     line_aux(i:i) = blank
                  enddo 
               endif 
               line_of_text = line_aux
               len_line = LEN(line_of_text)
c***  
c     we've got it, finally !
c***  
cc wrong
               count = count + 1
               lines_of_text(count) = line_of_text
               name_var_len(count) = INDEX(line_of_text,eq_char) - 1
            endif 
         endif
      enddo 
c***  
  100 nlines = count
c***  
c---      call GPsort
c***
      return 
      end
c***********************************************************************
      SUBROUTINE read_file(name)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      CHARACTER*(*) name
      CHARACTER*128 line_of_text, line_aux
      CHARACTER*1 blank, eq_char, tab, com_1, com_2, char
      DATA blank/' '/, eq_char/'='/, tab/'\t'/, com_1/'#'/, com_2/'!'/
      INTEGER i, k, m, count, len_line
      INTRINSIC LEN, INDEX 
c***  
      count = nlines
      open(unit=10,file=name,status='old',ERR=100)  
      do k = 1,nlines_max
         read (10,'(A128)',ERR=100,END=100) line_of_text  
         len_line = LEN(line_of_text)
         m = INDEX(line_of_text,eq_char)
c***  
c     lines with no '=' are discarded
c***  
         if (m .gt. 0) then
c***  
c     the text line is copied into the initialised auxiliary line,
c     all blanks and tabs are removed in the process
c***  
            do i=1,len_line
               line_aux(i:i) = blank
            enddo 
            m = 1
            do i = 1,len_line
               char = line_of_text(i:i)
               if (char .ne. blank .and. char .ne. tab) then
                  line_aux(m:m) = char
                  m = m + 1
               endif
            enddo 
c***  
c     discard line if 1st char is a '=', '#', '!' or ' ' and also 
c     if 1st character after '=' is a '#', '!' or ' '
c***  
            m = 1
            char = line_aux(1:1)
            if (char .eq. eq_char .or. char .eq. com_1 .or.
     %          char .eq. com_2 .or. char .eq. blank) m = 0 
            i = INDEX(line_aux,eq_char)
            char = line_aux(i+1:i+1)   
            if (char .eq. com_1 .or. char .eq. com_2 .or.
     %          char .eq. blank) m = 0
c***  
c     eventual comments at end of line are discarded
c***  
            if (m .gt. 0) then
               m = INDEX(line_aux,com_1)
               if ( m .eq. 0) m = INDEX(line_aux,com_2)
               if (m .gt. 0) then
                  do i = m,len_line
                     line_aux(i:i) = blank
                  enddo 
               endif 
               line_of_text = line_aux
               len_line = LEN(line_of_text)
c***  
c     we've got it, finally !
c***  
               count = count + 1
               lines_of_text(count) = line_of_text
               name_var_len(count) = INDEX(line_of_text,eq_char) - 1
            endif 
         endif
      enddo 
c***  
  100 nlines = count
      close (10)
c***  
c---      call GPsort
c***
      return 
      end
