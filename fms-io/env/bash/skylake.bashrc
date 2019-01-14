source $MODULESHOME/init/bash
module load intel/2018_up4
module load impi/2018_up4
module load hdf5/1.10.1
module load netcdf/4.6.1
function make() {
    /usr/bin/make -f ./Makefile.skylake "$@"
}
