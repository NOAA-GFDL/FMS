#!/bin/sh

#***********************************************************************
#*                             Apache License 2.0
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* Licensed under the Apache License, Version 2.0 (the "License");
#* you may not use this file except in compliance with the License.
#* You may obtain a copy of the License at
#*
#*     http://www.apache.org/licenses/LICENSE-2.0
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
#* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#* PARTICULAR PURPOSE. See the License for the specific language
#* governing permissions and limitations under the License.
#***********************************************************************

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/fms2_io directory.

# Author: Ed Hartnett 6/10/20
#
# Set common test settings.
. ../test-lib.sh

# Create and enter output directory
output_dir

# make an input.nml for mpp_init to read
touch input.nml

# run the tests
test_expect_success "Test the filename_appendix functionality" '
  mpirun -n 1 ../test_file_appendix
'
test_expect_success "Simple IO test" '
  mpirun -n 6 ../test_io_simple
'
test_expect_success "Test the get_mosaic_tile_grid functionality" '
  mpirun -n 6 ../test_get_mosaic_tile_grid
'
test_expect_success "Test the get_valid is_valid functionality single PE" '
  mpirun -n 1 ../test_get_is_valid
'
test_expect_success "Test the get_valid is_valid functionality multiple PE" '
  mpirun -n 1 ../test_get_is_valid
'
test_expect_success "Test the unlimited compressed axis functionality" '
  mpirun -n 6 ../test_unlimit_compressed
'

test_expect_success "Test the chunksizes functionality with the default behavior" '
  mpirun -n 1 ../test_chunksizes
'

if test ! -z "$ncdump_skip" ; then
  test_expect_success "test_chunksizes_netcdf4.res.nc should be chunked" '
    ncdump -hsv var1 test_chunksizes_netcdf4.res.nc | grep "ChunkSizes"
  '
  test_expect_failure "test_chunksizes_classic.res.nc should not be chunked" '
    ncdump -hsv var1 test_chunksizes_classic.res.nc | grep "ChunkSizes"
  '
  test_expect_failure "test_chunksizes_64bit.res.nc should not be chunked" '
    ncdump -hsv var1 test_chunksizes_64bit.res.nc | grep "ChunkSizes"
  '
  test_expect_failure "test_chunksizes.res.nc should not be chunked" '
    ncdump -hsv var1 test_chunksizes.res.nc | grep "ChunkSizes"
  '
fi

cat <<_EOF > input.nml
&fms2_io_nml
  netcdf_default_format = "netcdf4"
/
_EOF
test_expect_success "Test the chunksizes functionality with netcdf4 as the default file format" '
  mpirun -n 1 ../test_chunksizes
'
if test ! -z "$ncdump_skip" ; then
  test_expect_success "test_chunksizes_netcdf4.res.nc should be chunked" '
    ncdump -hsv var1 test_chunksizes_netcdf4.res.nc | grep "ChunkSizes"
  '
  test_expect_failure "test_chunksizes_classic.res.nc should not be chunked" '
    ncdump -hsv var1 test_chunksizes_classic.res.nc | grep "ChunkSizes"
  '
  test_expect_failure "test_chunksizes_64bit.res.nc should not be chunked" '
    ncdump -hsv var1 test_chunksizes_64bit.res.nc | grep "ChunkSizes"
  '
  test_expect_success "test_chunksizes.res.nc should be chunked" '
    ncdump -hsv var1 test_chunksizes.res.nc | grep "ChunkSizes"
  '
fi
cat <<_EOF > input.nml
&fms2_io_nml
  deflate_level = 3
  shuffle = .true.
/
_EOF
test_expect_success "Test the deflate level and shuffle functionality with the default behavior" '
  mpirun -n 1 ../test_chunksizes
'
if test ! -z "$ncdump_skip" ; then
  test_expect_success "test_chunksizes_netcdf4.res.nc should be compressed" '
    ncdump -hsv var1 test_chunksizes_netcdf4.res.nc | grep "DeflateLevel"
  '
  test_expect_failure "test_chunksizes_classic.res.nc should not be compressed" '
    ncdump -hsv var1 test_chunksizes_classic.res.nc | grep "DeflateLevel"
  '
  test_expect_failure "test_chunksizes_64bit.res.nc should not be compressed" '
    ncdump -hsv var1 test_chunksizes_64bit.res.nc | grep "DeflateLevel"
  '
  test_expect_failure "test_chunksizes.res.nc should not be compressed" '
    ncdump -hsv var1 test_chunksizes.res.nc | grep "DeflateLevel"
  '
fi
cat <<_EOF > input.nml
&fms2_io_nml
  deflate_level = 3
  shuffle = .true.
  netcdf_default_format = "netcdf4"
/
_EOF
test_expect_success "Test the deflate level and shuffle functionality with netcdf4 as the default file format" '
  mpirun -n 1 ../test_chunksizes
'
if test ! -z "$ncdump_skip" ; then
  test_expect_success "test_chunksizes_netcdf4.res.nc should be compressed" '
    ncdump -hsv var1 test_chunksizes_netcdf4.res.nc | grep "DeflateLevel"
  '
  test_expect_failure "test_chunksizes_classic.res.nc should not be compressed" '
    ncdump -hsv var1 test_chunksizes_classic.res.nc | grep "DeflateLevel"
  '
  test_expect_failure "test_chunksizes_64bit.res.nc should not be compressed" '
    ncdump -hsv var1 test_chunksizes_64bit.res.nc | grep "DeflateLevel"
  '
  test_expect_success "test_chunksizes.res.nc should be compressed" '
    ncdump -hsv var1 test_chunksizes.res.nc | grep "DeflateLevel"
  '
fi
test_done
