# change these lines for your local machine
ROOTDIR= $(HOME)/Sources/
INCLUDE = $(HOME)/include
DEST  = $(HOME)/lib


#MACHINE = $(shell uname)
#MACHINE = IRIX64

#OBJ = ffetch.o getpars.o ifetch.o isapipe.o fetch.o err.o redin.o redout.o\
#      reed.o rite.o sfetch.o dfetch.o snap.o alloc.o nrutil.o pfplib.o
OBJ = ffetch.o getpars.o ifetch.o  isapipe.o  sfetch.o  alloc.o redin.o\
      dfetch.o pfplib.o nrutil.o

CFLAGS = -DCHECKMEMORY -O2 -I.
CFLAGS = -O2 -I.
CC = gcc-14
#CFLAGS = -g -I.
#CFLAGS = -DCHECKMEMORY -g -I.
Lib = -L$(HOME)/lib -lpfp

libpfp$(MACHINE): $(OBJ)
	@echo MACHINE = $(MACHINE)
	ar cr libpfp$(MACHINE).a $(OBJ)
	cp nrutil.h $(INCLUDE)/nrutil.h
	cp pfplib.h $(INCLUDE)/
	mv libpfp$(MACHINE).a $(HOME)/lib
	cp checkMemory.h $(INCLUDE)/
#	ranlib libpfp$(MACHINE).a

install: libpfp$(MACHINE)
	mv -f libpfp$(MACHINE).a $(DEST)/libpfp$(MACHINE).a
	ranlib $(DEST)/libpfp$(MACHINE).a
#cp pfplib.h $(INCLUDE)/pfplib.h
	cp nrutil.h $(INCLUDE)/nrutil.h
	cp checkMemory.h $(INCLUDE)/

clean: 
	rm -f *.o libpfp$(MACHINE).a

test: install  test.o
	cc $(CFLAGS) -o test test.o  $(Lib)
