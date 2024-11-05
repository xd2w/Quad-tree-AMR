c     $Author: raus $
c     $Revision: 1.4 $
c     $Log: GP_fetch.f,v $
c Revision 1.4  1993/07/21  20:03:15  raus
c if the variable is read, now return 1
c
c Revision 1.3  1993/07/20  15:04:51  raus
c now 2 sets of fetch functions are available
c
c Revision 1.2  1993/04/21  16:52:26  raus
c size of the char function up to 128
c
c Revision 1.1  1993/04/20  17:41:37  raus
c 1st save; find_string to findstring
c
c***********************************************************************
      CHARACTER*128 FUNCTION GP_SFETCH(name_variable,def_val,message)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      CHARACTER*(*) name_variable, def_val, message
      CHARACTER*128 GETSTRING
      INTEGER line, FINDSTRING
      EXTERNAL FINDSTRING, GETSTRING
      INTRINSIC LEN
c***
      line = FINDSTRING(name_variable)
      if (line .eq. 0) then
         if (message .eq. 'DEFAULT') then
            GP_SFETCH = def_val
         elseif (message .eq. 'WARNING') then
            GP_SFETCH = def_val
            print*,'WARNING: using for ', name_variable,
     %           ' the default value ', def_val
         elseif (message .eq. 'SEVERE') then
            print*,'SEVERE: the value for ', name_variable,
     %           ' is not in the input files. STOP'
            stop
         else 
            print*,'allowed message values: DEFAULT, WARNING,',
     %           ' SEVERE.  STOP'
            stop
         endif 
      else 
         GP_SFETCH = GETSTRING(line)
      endif 
c***
c---       write (*,*) SFETCH
c***
      return 
      end
c***
c***********************************************************************
      INTEGER FUNCTION SFETCH(name_variable,string)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      CHARACTER*(*) name_variable, string
      CHARACTER*128 GETSTRING
      INTEGER line, FINDSTRING
      EXTERNAL FINDSTRING, GETSTRING
      INTRINSIC LEN
c***
      line = FINDSTRING(name_variable)
      if (line .eq. 0) then
         SFETCH = -1
         string = ' '
      else 
         SFETCH = 1
         string = GETSTRING(line)
      endif 
c***
c---       write (*,*) SFETCH
c***
      return 
      end
c***********************************************************************
      DOUBLE PRECISION FUNCTION GP_DFETCH(name_variable,def_val,message)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      DOUBLE PRECISION  def_val, GETDOUBLE
      CHARACTER*(*) name_variable, message
      INTEGER line, FINDSTRING
      EXTERNAL GETDOUBLE, FINDSTRING
c***
      line = FINDSTRING(name_variable)
      if (line .eq. 0) then
         if (message .eq. 'DEFAULT') then
            GP_DFETCH = def_val
         elseif (message .eq. 'WARNING') then
            GP_DFETCH = def_val
            print*,'WARNING: using for ', name_variable,
     %           ' the default value ', def_val
         elseif (message .eq. 'SEVERE') then
            print*,'SEVERE: the value for ', name_variable,
     %           ' is not in the input files. STOP'
            stop
         else 
            print*,'allowed message values: DEFAULT, WARNING,',
     %           ' SEVERE.  STOP'
            stop
         endif 
      else 
         GP_DFETCH = GETDOUBLE(line)
      endif 
c***
c---       write (*,*) GP_DFETCH
c***
      return
      end
c***********************************************************************
      INTEGER FUNCTION DFETCH(name_variable,dd)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      DOUBLE PRECISION  dd, GETDOUBLE
      CHARACTER*(*) name_variable
      INTEGER line, FINDSTRING
      EXTERNAL GETDOUBLE, FINDSTRING
c***
      line = FINDSTRING(name_variable)
      if (line .eq. 0) then
         DFETCH = -1
         dd = 0.d0
      else 
         DFETCH = 1
         dd = GETDOUBLE(line)
      endif 
c***
c---       write (*,*) DFETCH
c***
      return
      end
c***********************************************************************
      REAL FUNCTION GP_RFETCH(name_variable,def_val,message)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      REAL  def_val, GETREAL
      CHARACTER*(*) name_variable, message
      INTEGER line, FINDSTRING
      EXTERNAL GETREAL, FINDSTRING
c***
      line = FINDSTRING(name_variable)
      if (line .eq. 0) then
         if (message .eq. 'DEFAULT') then
            GP_RFETCH = def_val
         elseif (message .eq. 'WARNING') then
            GP_RFETCH = def_val
            print*,'WARNING: using for ', name_variable,
     %           ' the default value ', def_val
         elseif (message .eq. 'SEVERE') then
            print*,'SEVERE: the value for ', name_variable,
     %           ' is not in the input files. STOP'
            stop
         else 
            print*,'allowed message values: DEFAULT, WARNING,',
     %           ' SEVERE.  STOP'
            stop
         endif 
      else 
         GP_RFETCH = GETREAL(line)
      endif 
c***
c---       write (*,*) GP_RFETCH
c***
      return
      end
c***********************************************************************
      INTEGER FUNCTION RFETCH(name_variable,rr)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      REAL  rr, GETREAL
      CHARACTER*(*) name_variable
      INTEGER line, FINDSTRING
      EXTERNAL GETREAL, FINDSTRING
c***
      line = FINDSTRING(name_variable)
      if (line .eq. 0) then
         RFETCH = -1
         rr = 0.e0
      else 
         RFETCH = 1
         rr = GETREAL(line)
      endif 
c***
c---       write (*,*) RFETCH
c***
      return
      end
c***********************************************************************
      INTEGER FUNCTION GP_IFETCH(name_variable,def_val,message)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      INTEGER  def_val, GETINTEGER
      CHARACTER*(*) name_variable, message
      INTEGER line, FINDSTRING
      EXTERNAL GETINTEGER, FINDSTRING
c***
      line = FINDSTRING(name_variable)
      if (line .eq. 0) then
         if (message .eq. 'DEFAULT') then
            GP_IFETCH = def_val
         elseif (message .eq. 'WARNING') then
            GP_IFETCH = def_val
            print*,'WARNING: using for ', name_variable,
     %           ' the default value ', def_val
         elseif (message .eq. 'SEVERE') then
            print*,'SEVERE: the value for ', name_variable,
     %           ' is not in the input files. STOP'
            stop
         else 
            print*,'allowed message values: DEFAULT, WARNING,',
     %           ' SEVERE.  STOP'
            stop
         endif 
      else 
         GP_IFETCH = GETINTEGER(line)
      endif 
c***
c---       write (*,*) GP_IFETCH
c***
      return
      end
c***********************************************************************
      INTEGER FUNCTION IFETCH(name_variable,ii)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      INTEGER  ii, GETINTEGER
      CHARACTER*(*) name_variable
      INTEGER line, FINDSTRING
      EXTERNAL GETINTEGER, FINDSTRING
c***
      line = FINDSTRING(name_variable)
      if (line .eq. 0) then
         IFETCH = -1
         ii = 0
      else 
         IFETCH = 1
         ii = GETINTEGER(line)
      endif 
c***
c---       write (*,*) IFETCH
c***
      return
      end
