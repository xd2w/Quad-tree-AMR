c***********************************************************************
c     common space containing the input lines and the length of 
c     each variable name
c***
      INTEGER nlines, nlines_max
      PARAMETER (nlines_max=100)
      INTEGER name_var_len(nlines_max)
      CHARACTER*128 lines_of_text(nlines_max)
      COMMON /PARMRD/ nlines, name_var_len, lines_of_text 
