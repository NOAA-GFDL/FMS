source $MODULESHOME/init/bash
module use -a /home/sdu/publicmodules
module load intel_compilers/15.0.0
module load netcdf-fortran/4.4.1
module load mpich/3.1.3
function make() {
    /usr/bin/make -f ./Makefile.gfdl-ws "$@"
}
