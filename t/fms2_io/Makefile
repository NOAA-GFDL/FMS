INCLUDEADD = -I../..
LDADD= ../../libfms.a
RPATH = 

all: test_atmosphere_io.x

argparse.o: argparse.F90
	$(FC) $(FFLAGS) -o $@ -c $<

setup.o: setup.F90 argparse.o
	$(FC) $(FFLAGS) $(INCLUDEADD) -o $@ -c $<

test_atmosphere_io.x: test_atmosphere_io.F90 setup.o argparse.o
	$(FC) $(FFLAGS) $(INCLUDEADD) -o $@ $^ $(LDFLAGS) $(RPATH) $(LDADD)

check: test_atmosphere_io.x
	rm -f *.nc *.nc.*
	mpirun -n 12 ./test_atmosphere_io.x
	rm -f *.nc
	mpirun -n 18 ./test_atmosphere_io.x -d

clean:
	rm -f argparse.o argparse.mod setup.o setup.mod test_atmosphere_io.x
	rm -f logfile.000000.out *.nc
