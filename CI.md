# CI Information

The CI workflows and checks, and associated container environments for this repository
are listed below.
Actions run via Github-hosted runners unless otherwise noted.
Required CI for pull requests are listed first.

## Pull Request CI and checks

### Build libFMS with autotools using GCC

Required GNU build test for all pull requests/pushes.
Runs `make distcheck` after configuring via GNU autotools.

Runs on a container image with spack installed dependencies, on top a rocky linux base.

Dockerfile for image is stored in the [HPC-ME repository](github.com/noaa-gfdl/hpc-me).

Container environment:
gcc            v13.2.0
mpich          v4.0.2
netcdf-c       v4.9.2
netcdf-fortran v4.6.1
autoconf       v2.69
libyaml        v0.2.5

`./configure` flags tested:
- `--disable-openmp`
- `--enable-mixed-mode`
- `--with-mpi=no` (disables unit testing)
- `--disable-setting-flags`
- `--with-yaml`
- `--enable-test-input=/home/unit_tests_input`


### Build libfms with cmake using GCC
Required GNU build test for all pull requests/pushes.
Runs `make` after configuring via cmake.

cmake flags:
- `-DOPENMP=on`
- `-DOPENMP=on`
- `-DWITH_YAML=on`
- `-D64BIT=on`


### Build libfms with autotools using Intel Oneapi Compilers

Required build test for all pull requests. Workflow will build hdf5, netcdf, and libyaml and cache them for reuse. Cache can be used for a week before clearing.

Test Environment:
intel-oneapi		v2025.3.0
hdf5						v1.14.6
netcdf-c				v4.9.3
netcdf-fortran	v4.6.2
libyaml					v0.2.5

### Build libfms with autotools using Intel Classic Compilers

This test is triggered weekly on Sundays @ midnight and uses the intel 2023.1 classic compilers (ie. ifort/icc).

Test Environment:
intel-oneapi		v2025.3.0
hdf5						v1.12.2
netcdf-c				v4.8.1
netcdf-fortran	v4.6.0
libyaml					v0.2.5

### libFMS lint tests
Required test for all pull requests.
Checks code for line lengths, tabs, and trailing whitespace in accordance with
the project's [style guide](https://github.com/NOAA-GFDL/FMS/blob/main/CODE_STYLE.md).
The action is hosted on github [here](https://github.com/NOAA-GFDL/simple_lint).

## Parallelworks CI
The following CI workflows run on self-hosted runners through the parallelworks platform.

### AM5 testing

On all pull requests, a full scale model test is compiled using the `c96L65_am5f11d12r0_amip` experiment from AM5 and run on a parallelworks cluster.

This test will compile with debug flags and openmp enabled for a basic regression test.
