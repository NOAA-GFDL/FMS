\page install Installation
# libFMS Installation

## External Requirements

The following external libraries are required when building libFMS

* NetCDF C and Fortran (77/90) headers and libraries
* Fortran 2003 standard compiler
* Fortran compiler that supports Cray Pointer
* MPI C and Fortran headers and libraries (optional)
* Linux or Unix style system

## Supported Compilers

We strive to have libFMS built with as many C and Fortran compilers on as many
Unix/Linux type systems as possible.  However, internally, we only use libFMS
compiled with the Intel compilers.  Some groups have had success using libFMS
when compiled with the GNU C and Fortran compilers.

## MPI Support

The default way to build libFMS is with MPI support.  We have found that using
the MPI aware compiler, or the MPI compiler wrappers (mpif90, mpicc, etc.), in
general, offer the best results when building with MPI support.  If you decide
to not use an MPI aware compiler, you should pass the include and library
locations to the build system.

libFMS can be built without MPI support (sometimes called "no-comm mode").  To
build libFMS without MPI support, pass to `configure` the `--disable-mpi` flag.

## Supported Build Systems

### GNU Autoconf/Automake

In many cases, running the shell command `./configure && make && make install`
will build and install this package.  Since libFMS requires the netCDF libraries,
you will likely need to tell configure where to find the C and Fortran netCDF
headers and libraries via the `CPPFLAGS`, `FCFLAGS` and `LDFLAGS` variables

```bash
> ./configure CPPFLAGS="-I/path/to/netcdf/include" \
              FCFLAGS="-I/path/to/netcdff/include" \
              LDFLAGS="-L/path/to/netcdf/lib -L/path/to/netcdff/lib"
```

When building from a GitHub clone, the user must run `autoreconf -i` before
running the `./configure` script.

The `./configure` script will guess as many compiler options as required to
build libFMS.  In some cases you may be required to supply additional compiler
flags for your system.  Additional compiler options can be passed to the build
with the `CPPFLAGS`, `CFLAGS` and `FCFLAGS` variable.  The `./configure` option
`--disable-setting-flags` will not guess required compiler flags.  Using this
configure option will require the user to give all required flags.

Run `./configure --help` to see other available configure options.

### CMake

To build using CMake, follow the instructions in
[CMAKE_INSTRUCTIONS.mk](https://github.com/NOAA-GFDL/FMS/blob/master/CMAKE_INSTRUCTIONS.md).
Currently the CMake configuration is the most restrictive build option as the
compiler flags are immutable.

### GFDL MKMF

Make Makefile ([MKMF](https://github.com/NOAA-GFDL/mkmf)) is the GFDL developed
makefile generator application.  MKMF uses a list of file, or a list of
directories that contain Fortran or C file, a Make template file, and writes
a Makefile with the correct dependencies to build a library or executable.  To
use MKMF to build libFMS, do the following.

1. Clone the MKMF repository from GitHub
   `git clone https://github.com/NOAA-GFDL/mkmf.git`
2. Add the `mkmf/bin` directory to your PATH
   `export PATH=/path/to/mkmf/bin:${PATH}`
   or
   `setenv PATH /path/to/mkmf/bin:${PATH}`
3. Run the `list_paths` command, on the FMS directory
   `list_paths path/to/FMS`
   This will create a path_names file in the current directory
4. Pick a MKMF template from the `mkmf/templates` directory that most closely
   matches you system and compiler.  Make any modifications as need
5. Run `mkmf`
   ```shell
   mkmf -c "-Duse_netCDF -Duse_libMPI" \
         -I /path/to/netcdff/include \
         -I /path/to/netcdf/include \
         -t /path/to/mkmf/template. \
         -p libfms.a \
         path_names
   ```
6. Once built, place the `libfms.a`, `*.mod`, `include/\*` files in your install
   location.

While MKMF has been the build system used by GFDL to build all the GFDL models
for some time, this build method requires the most user intervention when
building libFMS on a new system.
