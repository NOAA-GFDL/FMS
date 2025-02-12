#!/bin/sh

#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************
#
# Copyright (c) 2019-2021 Ed Hartnett, Uriel Ramirez, Seth Underwood

# Set common test settings.
. ../test-lib.sh

output_dir
rm -rf data_table data_table.yaml input.nml input_base.nml

for KIND in r4 r8
do
# Run tests with input if enabled
# skips if built with yaml parser(tests older behavior)
if test ! -z "$test_input_path" && test ! -z "$parser_skip"  ; then
  cat <<_EOF > input.nml
_EOF

  cp -r $test_input_path/data_override/INPUT .
  cat <<_EOF > diag_table
test_data_override
1 3 1 0 0 0

#output files
"test_data_override",  -1, "days", 1, "days", "time"

#output variables
"test_data_override_mod", "sst", "sst", "test_data_override",  "all", .false., "none", 2
"test_data_override_mod", "ice", "ice", "test_data_override",  "all", .false., "none", 2
_EOF
  cat <<_EOF > data_table
"ICE", "sst_obs",  "SST", "INPUT/sst_ice_clim.nc", .false., 300.0
"ICE", "sic_obs",  "SIC", "INPUT/sst_ice_clim.nc", .false., 300.0
"OCN", "sst_obs",  "SST", "INPUT/sst_ice_clim.nc", .false., 300.0
"LND", "sst_obs",  "SST", "INPUT/sst_ice_clim.nc", .false., 300.0
_EOF

  test_expect_success "data_override on cubic-grid with input (${KIND})" '
    mpirun -n 6 ../test_data_override_${KIND}
  '

cat <<_EOF > input.nml
&test_data_override_nml
   test_num=2
/
_EOF

  test_expect_success "data_override on latlon-grid with input (${KIND})" '
    mpirun -n 6 ../test_data_override_${KIND}
  '
fi
done
rm -rf INPUT *.nc # remove any leftover files to reduce size

test_done