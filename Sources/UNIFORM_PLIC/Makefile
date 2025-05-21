
OBJS = start.o nrutil.o miscUtils.o uniform.o init.o potentialSolver.o plot.o plicx.o plicy.o plic.o 
OBJS = start.o nrutil.o miscUtils.o uniform.o init.o LBMSolver.o plot.o plicx.o plicy.o plic.o 
OBJS = start.o nrutil.o miscUtils.o uniform.o init.o PossionSolver.o plot.o plicx.o plicy.o plic.o 
CLIB=  -L$(HOME)/lib

CFLAGS = -I$(HOME)/include -O3

UNIFORM_PLIC: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(CLIB) -lm
	# @echo "Build complete. Run './UNIFORM_PLIC' to execute."
	# @echo "To clean up, run 'make clean'."

run: UNIFORM_PLIC
	@echo "Running the program..."
	@./UNIFORM_PLIC
	@echo "Execution complete. Check the output files in the DATA directory."

clean:
	rm -f $(OBJS) UNIFORM_PLIC
	@echo "Cleaned up object files and executable."
	@echo "To build again, run 'make'."
	@echo "To run the program, execute './UNIFORM_PLIC'."

plot:
	cd DATA && \
	gnuplot ../plot_DATA.plt

plot_vof:
	cd DATA && \
	gnuplot ../plot_VOF.plt

clear_data:
	rm DATA/*