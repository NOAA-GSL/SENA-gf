FC=gfortran
FFLAGS=-DNDEBUG -fbacktrace -ffp-contract=off -fno-sign-zero -g -I$(NETCDF)/include
LDFLAGS=-L$(NETCDF)/lib -lnetcdf -lnetcdff
OBJS=machine.o mt19937.o gf_utils.o cu_gf_driver.o cu_gf_deep.F90 cu_gf_sh.F90

all: gf_test.exe

machine.o: machine.F90
	$(FC) -c $(FFLAGS) -o $@ $<

mt19937.o: mt19937.F90
	$(FC) -c $(FFLAGS) -o $@ $<

gf_utils.o: gf_utils.F90
	$(FC) -c $(FFLAGS) -o $@ $<

cu_gf_deep.o: cu_gf_deep.F90
	 $(FC) -c $(FFLAGS) -o $@ $<

cu_gf_sh.o: cu_gf_sh.F90
	 $(FC) -c $(FFLAGS) -o $@ $<

cu_gf_driver.o: cu_gf_driver.F90 cu_gf_deep.o cu_gf_sh.o
	$(FC) -c $(FFLAGS) -o $@ $<

gf_test.exe: $(OBJS) gf_test.F90
	$(FC) -o $@ $(OBJS) $(LDFLAGS) gf_test.F90
clean:
	rm -f *.o *.mod *.exe
