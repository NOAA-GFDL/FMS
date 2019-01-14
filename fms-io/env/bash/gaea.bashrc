source $MODULESHOME/init/bash
module load cray-netcdf
module load cray-hdf5
function make() {
    /usr/bin/make -f Makefile.gaea "$@"
}
