
#   These are the objects required for the ftt code
BALANCE = balanceFTT.o balanceFTTxP.o balanceFTTxN.o balanceFTTyP.o balanceFTTyN.o

FTT_ARM = FTT_ARM.o ftt.o $(BALANCE) box.o splitCell.o fttStatistics.o cell_OctFlag.o \
	restrFlag.o restrField.o flagOct.o binCollection.o binCollectionAtLevel.o garbageCollection.o \
	seekCell.o checkOctTree.o establishNb.o checkNb.o octPropagation.o plotFTTInterf_curvature_top_down.o\
	plotCellInterf.o plotVOF.o plotCurvature.o plic_curvature_top_down.o vol2.o \
	octTreeXSwp_ellipse.o octTreeYSwp_ellipse.o reMesh.o miscUtils.o computeVelocity.o computeKappa.o

# FTT_ARM = FTT_ARM.o ftt.o $(BALANCE) box.o splitCell.o fttStatistics.o cell_OctFlag.o \
# 	restrFlag.o restrField.o flagOct.o binCollection.o binCollectionAtLevel.o garbageCollection.o \
# 	seekCell.o checkOctTree.o establishNb.o checkNb.o octPropagation.o plotFTTInterf_curvature_top_down.o\
# 	plotCellInterf.o plotVOF.o plotCurvature.o plic_curvature_top_down_HILB.o vol2.o \
# 	octTreeXSwp_ellipse.o octTreeYSwp_ellipse.o reMesh.o miscUtils.o computeVelocity.o computeKappa.o

# ICBC = initialize.o initVOF.o computeVOF.o markVOF.o initFTT.o refineFTT.o bcOctTree.o
# ICBC = initialize.o initVOF.o computeVOF.o markVOF.o initFTT_Yung.o refineFTT_Yung.o bcOctTree_Yung.o
#ICBC = initialize.o initVOF.o computeVOF.o markVOF.o initFTT_Dodd.o refineFTT_Dodd.o bcOctTree_Dodd.o
ICBC = initialize.o initVOF.o computeVOF_ellipse.o markVOF.o initFTT_Yung.o refineFTT_ellipse.o bcOctTree_Yung.o

IO  = readData.o plotFTTOct.o plotFTT.o plotSFC.o printFTTCell.o fGrad.o drawNgbCells.o \
	drawChCells.o drawPrCells.o saveFTT.o readFTT.o initMemory.o computeCoord.o plotThIntf.o

OBJECTS = $(FTT_ARM) $(ICBC) $(IO)
CFLAGS= -O2  -I$(HOME)/include
# CFLAGS= -g  -I$(HOME)/include
# LIB= alloc.o   ffetch.o   ifetch.o   nrutil.o  redin.o dfetch.o  getpars.o  isapipe.o  pfplib.o  sfetch.o
CLIB=  -L$(HOME)/lib -lpfp
# CLIB=  -L -lpfp

# CC=gcc-14

FTT_ARM: $(OBJECTS) #$(LIB)
	#$(CC) $(CFLAGS) -o FTT_ARM $(OBJECTS)  $(LIB) -lm
	$(CC) $(CFLAGS) -o FTT_ARM $(OBJECTS)  $(CLIB) -lm
	mv FTT_ARM Examples

run: FTT_ARM
	cd Examples && ./FTT_ARM ftt.par

plot : FTT_ARM
	cd Examples/DATA && gnuplot ~/include/plot_DATA.plt


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
