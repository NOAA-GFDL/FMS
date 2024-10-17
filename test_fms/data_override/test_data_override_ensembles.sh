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

cat <<_EOF > data_table.ens_01.yaml
data_table:
 - grid_name: OCN
   fieldname_in_model: runoff
   override_file:
   - fieldname_in_file: runoff
     file_name: INPUT/runoff.daitren.clim.1440x1080.v20180328_ens_01.nc
     interp_method: none
   factor: 1.0
_EOF

cat <<_EOF > data_table.ens_02.yaml
data_table:
 - grid_name: OCN
   fieldname_in_model: runoff
   override_file:
   - fieldname_in_file: runoff
     file_name: INPUT/runoff.daitren.clim.1440x1080.v20180328_ens_02.nc
     interp_method: none
   factor: 1.0
_EOF

cat <<_EOF > input_base.nml
&data_override_nml
  use_data_table_yaml = .True.
/

&test_data_override_ongrid_nml
  test_case = 5
  write_only = .False.
/

&ensemble_nml
   ensemble_size = 2
/
_EOF

#The test only runs with yaml
if [ -z $parser_skip ]; then
  for KIND in r4 r8
  do
    rm -rf INPUT/.
    sed 's/write_only = .False./write_only = .True./g' input_base.nml > input.nml
    test_expect_success "Creating input files (${KIND})" '
      mpirun -n 12 ../test_data_override_ongrid_${KIND}
    '

    cp input_base.nml input.nml
    test_expect_success "test_data_override with two ensembles  -yaml (${KIND})" '
      mpirun -n 12 ../test_data_override_ongrid_${KIND}
      '
  done

cat <<_EOF > data_table.yaml
data_table:
 - grid_name: OCN
   fieldname_in_model: runoff
   override_file:
   - fieldname_in_file: runoff
     file_name: INPUT/runoff.daitren.clim.1440x1080.v20180328_ens_02.nc
     interp_method: none
   factor: 1.0
_EOF

    test_expect_failure "test_data_override with both data_table.yaml and data_table.ens_xx.yaml files" '
      mpirun -n 12 ../test_data_override_ongrid_${KIND}
      '
rm -rf INPUT
fi
test_done
