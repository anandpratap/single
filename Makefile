FC = ifort
CC = icc
CCFLAGS = -c
FCFLAGS = -r8 -L/usr/lib/ -llapack -lblas -debug -traceback -check bounds -g
OBJECTS = modules.o main.o step.o primvars.o primvars_.o residual.o compute_dt.o sigma.o roeflux.o reconstruct.o io.o calc_residual.o objective.o\
	diffsizes.o roeflux_bq.o residual_bq.o objective_bq.o primvars_bq.o sigma_bq.o adBuffer.o adStack.o
BIN = bin/
SRC = src/

mflame:$(OBJECTS)
	$(FC) -o single $(OBJECTS) -L/usr/lib/ -llapack -lblas -g -traceback -debug
	mv single bin/.

TPN             = /home/anandps/local/tapenade3.6/bin/tapenade
TPNFLAGS        = -r8 -backward -difffuncname "_bq"

adjoint:
	$(TPN) $(TPNFLAGS) src/modules.f90 src/residual.f90 src/sigma.f90 src/roeflux.f90 -outvars "res" -vars "q alpha" -head "residual"
	$(TPN) $(TPNFLAGS) src/modules.f90 src/primvars.f90 src/objective.f90 -outvars "obj" -vars "q alpha" -head "objective"
	@bash cleanup.bash
	mv *_bq* src/
	rm *_cb*
	rm  *_cd*
	rm -rf tapenadehtml

%.o:$(SRC)%.f90
	$(FC) $(FCFLAGS) -c $<

%.o:$(SRC)%.f
	$(FC) $(FCFLAGS) -c $<

%.o:$(SRC)%.c
	$(CC) $(CCFLAGS) -c $<


clean:
	rm *.o bin/single *.mod *_bq* *_cb* src/*_bqv* src/*_bq* *_cd*
