FC = ifort
CC = icc
CCFLAGS = -c -openmp -lpthread -I${MKLROOT}/include 
MKLROOT  = /opt/intel/composer_xe_2013_sp1.3.174/mkl
FCFLAGS = -r8 -openmp -lpthread -I${MKLROOT}/include
OBJECTS = modules.o main.o step.o primvars.o primvars_.o residual.o compute_dt.o sigma.o roeflux.o reconstruct.o io.o calc_residual.o objective.o\
	initialize.o metrics.o sensitivity.o adjoint.o diffsizes.o roeflux_bq.o residual_bq.o objective_bq.o primvars_bq.o sigma_bq.o adBuffer.o adStack.o
BIN = bin/
SRC = src/

TPN             = /home/anandps/local/tapenade3.6/bin/tapenade
TPNFLAGS        = -r8 -backward -difffuncname "_bq"

single:$(OBJECTS)
	$(FC) -o single $(OBJECTS) -openmp -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm
	mv single $(BIN).
	cp $(BIN)single examples/.

adjoint:
	$(TPN) $(TPNFLAGS) src/modules.f90 src/residual.f90 src/sigma.f90 src/roeflux.f90 -outvars "res" -vars "q alpha" -head "residual"
	$(TPN) $(TPNFLAGS) src/modules.f90 src/primvars.f90 src/objective.f90 -outvars "obj" -vars "q alpha" -head "objective"
	@bash cleanup.bash
	mv *_bq* src/

%.o:$(SRC)%.f90
	$(FC) $(FCFLAGS) -c $<

%.o:$(SRC)%.f
	$(FC) $(FCFLAGS) -c $<

%.o:$(SRC)%.c
	$(CC) $(CCFLAGS) -c $<


clean:
	rm *.o bin/single *.mod 

clean_adjoint:
	rm src/*_bq*