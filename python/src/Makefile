FC=ifort

# Use same flags as SCM uses in non-debug mode
FFLAGS=-g -traceback -fpp -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -nowarn -sox -align array64byte -qno-opt-dynamic-align -real-size 64 -O2 -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3 -no-prec-div -no-prec-sqrt -march=core-avx2 -O2 -fPIC

LDFLAGS=-L$(netcdf_c_ROOT)/lib -lnetcdf -L$(netcdf_fortran_ROOT)/lib -lnetcdff
OBJS=physcons.o machine.o cu_gf_driver.o cu_gf_deep.o cu_gf_sh.o cu_gf_io.o

all: gf_test.exe

machine.o: machine.F90
	$(FC) -c $(FFLAGS) -o $@ $<

physcons.o: physcons.F90 machine.o
	$(FC) -c $(FFLAGS) -o $@ $<

cu_gf_io.o: cu_gf_io.F90
	$(FC) -c $(FFLAGS) -o $@ $<

cu_gf_deep.o: cu_gf_deep.F90
	 $(FC) -c $(FFLAGS) -o $@ $<

cu_gf_sh.o: cu_gf_sh.F90
	 $(FC) -c $(FFLAGS) -o $@ $<

cu_gf_driver.o: cu_gf_driver.F90 cu_gf_deep.o cu_gf_sh.o cu_gf_io.o
	$(FC) -c $(FFLAGS) -o $@ $<

gf_test.exe: $(OBJS) gf_test.F90
	$(FC) -o $@ $(OBJS) $(LDFLAGS) $(FFLAGS) gf_test.F90

clean:
	rm -f *.o *.mod *.exe
