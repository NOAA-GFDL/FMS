# CI Information

The CI workflows and checks, and associated container environments for this repository
are listed below.
Actions run via Github-hosted runners unless otherwise noted.
Required CI for pull requests are listed first.

## Pull Request CI and checks

### Build libFMS with autotools
Required GNU build test for all pull requests/pushes.
Runs `make distcheck` after configuring via GNU autotools.

Container environment:
gcc            v7.3.0
mpich          v3.3a2
netcdf         v4.6.0
netcdf-fortran v4.4.4
autoconf       v2.69

container hosted at [noaagfdl/ubuntu_libfms_gnu:latest](https://hub.docker.com/r/noaagfdl/ubuntu_libfms_gnu)

`./configure` flags:
- `--enable-openmp`
- `--disable-openmp`
- `--enable-mixed-mode`
- `--disable-setting-flags`
- `--with-yaml`

### Build libfms with cmake
Required GNU build test for all pull requests/pushes.
Runs `make` after configuring via cmake.

Container environment:
gcc            v7.3.0
mpich          v3.3a2
netcdf         v4.6.0
netcdf-fortran v4.4.4
cmake          v3.22.0

container hosted at [noaagfdl/ubuntu_libfms_gnu:latest](https://hub.docker.com/r/noaagfdl/ubuntu_libfms_gnu)

cmake flags:
- `-DOPENMP=on`
- `-DOPENMP=on`
- `-DWITH_YAML=on`
- `-D64BIT=on`

### libFMS lint tests
Required test for all pull requests.
Checks code for line lengths, tabs, and trailing whitespace in accordance with
the project's [style guide](https://github.com/NOAA-GFDL/FMS/blob/main/CODE_STYLE.md).
The action is hosted on github [here](https://github.com/NOAA-GFDL/simple_lint).

## Parallelworks CI
The following CI workflows run on self-hosted runners through the parallelworks platform.
### Pull Request CI libFMS with intel
Optional(does not need to pass to merge) intel build test hosted on the parallelworks platform.
Runs `make check` with intel 18 and 21 compilers for all pull requests.

### Tag CI libFMS with AM4 regression
On alpha or beta tag creation, compiles and runs full AM4 model regression testing using the new FMS tag on parallelworks.
