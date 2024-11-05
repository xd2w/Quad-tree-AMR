OBJ = GP_find.o GP_read.o GP_term.o GP_fetch.o GP_get.o GP_sort.o 


#MACHINE = $(shell uname)
#MACHINE = IRIX64

INCLUDES = GP_common.h
HOME=/home/jie
LIB = $(HOME)/lib
#AR = /usr/bin/ar 
FFLAGS = -O2
FC = gfortran

all: libgp$(MACHINE).a 

libgp$(MACHINE).a: $(OBJ)
	echo $(AR)
	$(AR) cr libgp$(MACHINE).a  $(OBJ)
	mv -f libgp$(MACHINE).a  $(LIB)/libgp$(MACHINE).a

clean:
	rm -f *.o libgp$(MACHINE).a

$(OBJ): $(INCLUDES)

