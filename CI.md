# CI Information

The CI workflows and checks, and associated container environments for this repository
are listed below.
Actions run via Github-hosted runners unless otherwise noted.
Required CI for pull requests are listed first.

## Pull Request CI and checks

### Build libFMS with autotools (GNU)

Required GNU build test for all pull requests/pushes.
Runs `make distcheck` after configuring via GNU autotools.

Runs on a container image with spack installed dependencies, on top a rocky linux base.

Dockerfile for image is stored at .github/workflows/Dockerfile.gnu for more specific information on the CI environment.

Container environment:
gcc            v12.3.0
mpich          v4.0.2
netcdf         v4.9.0
netcdf-fortran v4.6.0
autoconf       v2.69
libyaml        v0.2.5

`./configure` flags tested:
- `--disable-openmp`
- `--enable-mixed-mode`
- `--with-mpi=no` (disables unit testing)
- `--disable-setting-flags`
- `--with-yaml`
- `--enable-test-input=/home/unit_tests_input`

### Build libFMS with autotools (intel)

Required Intel build test for all pull requests.
Runs `make check` after configuring via autotools.

This runs on an intel-hosted container image from dockerhub.

To access the netcdf and libyaml dependencies, it builds and caches the resulting libaries for reuse.
The cache is cleared after not being used for a week.

`./configure` flags tested:
- `--disable-deprecated-io`
- `--enable-deprecated-io`


### Build libfms with cmake
Required GNU build test for all pull requests/pushes.
Runs `make` after configuring via cmake.

This uses the same container image as the GNU autotools CI.

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

### Build FMScoupler with FMS

Runs on pull requests. This workflow pulls in the main branch of the [FMScoupler repository](github.com/noaa-gfdl/FMScoupler),
and then run's its test script using the pull request's version of FMS.

The coupler's testing script runs a "null model" that uses simple placeholders for the components (ie. atmos_null, ice_null repositories).
It also uses [`mkmf`](github.com/noaa-gfdl/mkmf), a GFDL created build tool to generate it's Makefile's.

It uses the same docker image as the autotools GNU ci.

## Miscellaneous

### Documentation Site
The `github_doc_site.yml` workflow uses the program doxygen to parse our documentation and create a searchable site.
This site is then stored as an artifact, and deployed to github-pages to be accessible at [noaa-gfdl.github.io/FMS].


### Version updates
The `version.yml` workflow appends the -dev to the version number used by autotools upon the creation of a published release.

It will create a new pull request with the change automatically.

Since it's created by a bot account, the pull request will need to be closed and re-opened in order to be merged (due to above CI requirements).
