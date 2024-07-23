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
[ ! -d "INPUT" ] && mkdir -p "INPUT"

cat <<_EOF > data_table.yaml
data_table:
- grid_name: OCN
  fieldname_in_model: runoff_obs
  override_file:
  - fieldname_in_file: runoff
    file_name: ./INPUT/bilinear_increasing.nc
    interp_method: bilinear
  factor: 1.0
- grid_name: OCN
  fieldname_in_model: runoff_obs_weights
  override_file:
  - fieldname_in_file: runoff
    file_name: ./INPUT/bilinear_increasing.nc
    interp_method: bilinear
    external_weights:
    - file_name: ./INPUT/remap_file.nc
      source: fregrid
  factor: 1.0
_EOF

cat <<_EOF > input.nml
&data_override_nml
  use_data_table_yaml = .True.
/

&test_data_override_ongrid_nml
  test_case = 4
  nlon = 5
  nlat = 6
  layout = 1, 2
/
_EOF

#The test only runs with yaml
if [ -z $parser_skip ]; then
  for KIND in r4 r8
  do
    rm -rf INPUT/.
    test_expect_success "test_data_override with and without weight files  -yaml (${KIND})" '
      mpirun -n 2 ../test_data_override_ongrid_${KIND}
      '
  done
fi

rm -rf INPUT
test_done
