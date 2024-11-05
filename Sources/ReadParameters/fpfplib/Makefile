OBJ = getpars.o ffetch.o ifetch.o dfetch.o
#HOME=/export/home/jie
INCLUDES = common.h undefine.h

#MACHINE = $(shell uname)

#  Remplir eventuellement la ligne ci dessous et la decommenter
Lib = -L$(HOME)/lib -lfpfp$(MACHINE)
LIB = $(HOME)/lib 

FC = gfortran
FFLAGS = -O2


all: libfpfp$(MACHINE).a install test
#all: libfpfp$(MACHINE).a install test

libfpfp$(MACHINE).a: $(OBJ)
	ar cr libfpfp$(MACHINE).a $(OBJ)
#	ranlib libfpfp$(MACHINE).a

test: test.o 
	$(FC) -o test $(FFLAGS) test.o $(Lib)

install: libfpfp$(MACHINE).a
	mv -f libfpfp$(MACHINE).a $(LIB)
#	ranlib $(LIB)/libfpfp$(MACHINE).a

clean:
	rm -f *.o libfpfp$(MACHINE).a

$(OBJ): $(INCLUDES)

