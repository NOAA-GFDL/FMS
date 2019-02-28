FC = ftn
FFLAGS = -g -O0 -r8 -i4 -traceback -qopenmp
INCLUDE = -I../..
LDFLAGS = ../../libfms.a
RPATH = 

all: test_fms2_io.x

argparse.o: argparse.F90
	${FC} ${FFLAGS} -o $@ -c $<

test_fms2_io.x: test_fms2_io.F90 argparse.o create_atmosphere_domain.inc create_land_domain.inc \
                create_ocean_domain.inc atmosphere_restart_file_test.inc \
                land_unstructured_restart_file_test.inc land_compressed_restart_file_test.inc \
                ocean_restart_file_test.inc
	${FC} ${FFLAGS} ${INCLUDE} -o $@ $< argparse.o ${LDFLAGS} ${RPATH}

clean:
	rm -f argparse.o argparse.mod test_fms2_io.x
