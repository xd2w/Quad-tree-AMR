
#   These are the objects required for the ftt code
BALANCE = balanceFTT.o balanceFTTxP.o balanceFTTxN.o balanceFTTyP.o balanceFTTyN.o

# FTT_ARM = FTT_ARM.o ftt.o $(BALANCE) box.o splitCell.o fttStatistics.o cell_OctFlag.o \ 
# 	restrFlag.o restrField.o flagOct.o binCollection.o binCollectionAtLevel.o garbageCollection.o \
# 	seekCell.o checkOctTree.o establishNb.o checkNb.o octPropagation.o plotFTTInterf.o plotCellInterf.o \
# 	 plic.o vol2.o octTreeXSwp.o octTreeYSwp.o reMesh.o

FTT_ARM = FTT_HILB.o ftt.o $(BALANCE) box.o splitCell.o fttStatistics.o cell_OctFlag.o \
	restrFlag.o restrField.o flagOct.o binCollection.o binCollectionAtLevel.o garbageCollection.o \
	seekCell.o checkOctTree.o establishNb.o checkNb.o octPropagation.o stream.o plic.o octTreeXSwp.o \
	octTreeYSwp.o reMesh.o plotFTTInterf.o plotCellInterf.o interpolation.o


# ICBC = initialize.o initVOF.o computeVOF.o markVOF.o initFTT.o refineFTT.o bcOctTree.o
# ICBC = initialize.o initVOF.o computeVOF.o markVOF.o initFTT_Yung.o refineFTT_Yung.o bcOctTree_Yung.o
#ICBC = initialize.o initVOF.o computeVOF.o markVOF.o initFTT_Dodd.o refineFTT_Dodd.o bcOctTree_Dodd.o

ICBC = initialize.o initVOF.o initPotential.o computePotential.o computeVOF.o \
	initFTT_Yung.o refineFTT_Yung.o bcOctTree_Yung.o

IO  = readData.o plotFTTOct.o plotFTT.o plotSFC.o printFTTCell.o drawNgbCells.o drawChCells.o drawPrCells.o saveFTT.o readFTT.o initMemory.o computeCoord.o plotCircle.o

OBJECTS = $(FTT_ARM) $(ICBC) $(IO)
CFLAGS= -O2  -I$(HOME)/include
CFLAGS= -g  -I$(HOME)/include
#LIB= alloc.o   ffetch.o   ifetch.o   nrutil.o  redin.o dfetch.o  getpars.o  isapipe.o  pfplib.o  sfetch.o
CLIB=  -L$(HOME)/lib -lpfp

FTT_ARM: $(OBJECTS) $(LIB)
	$(CC) $(CFLAGS) -o FTT_HILB $(OBJECTS)  $(CLIB) -lm
	mv FTT_HILB Examples

intp_test:
	gcc -g  -I/homes/da533/include -o intp_test.o intp_test.c interpolation.o -L/homes/da533/lib -lpfp  -lm


tousb:
	rm -f /media/3FAF-7E76/FTT_ARM/*.h
	rm -f /media/3FAF-7E76/FTT_ARM/*.c
	rm -f /media/3FAF-7E76/FTT_ARM/Makefile
	cp *.c /media/3FAF-7E76/FTT_ARM/.         
	cp *.h /media/3FAF-7E76/FTT_ARM/.         
	cp Makefile /media/3FAF-7E76/FTT_ARM/.
	cp Examples/*.par /media/3FAF-7E76/FTT_ARM/Examples/.

fromusb:
	cp  /media/3FAF-7E76/FTT_ARM/*.c .
	cp  /media/3FAF-7E76/FTT_ARM/*.h . 
	cp /media/3FAF-7E76/FTT_ARM/Makefile . 
	cp /media/3FAF-7E76/FTT_ARM/Examples/*par Examples/.

clean:
	rm -f *.o core 
