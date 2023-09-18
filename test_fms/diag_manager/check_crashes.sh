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

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/data_override directory.

# Set common test settings.

printf "&check_crashes_nml \n checking_crashes = .true. \n/" | cat > input.nml
sed '/tile/d' diag_table.yaml_base > diag_table.yaml
test_expect_failure "Missing tile when using the 'index' grid type" '
  mpirun -n 1 ../test_diag_yaml
'

sed '/new_file_freq: 6 hours/new_file_freq: 6/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "Missing new_file_freq_units when using new_file_freq_units" '
  mpirun -n 1 ../test_diag_yaml
'

sed 's/new_file_freq: 6 hours/new_file_freq: 6 mullions/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "new_file_freq_units is not valid" '
  mpirun -n 1 ../test_diag_yaml
'

sed '/file_duration: 12 hours/file_duration: 12/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "Missing file_duration_units when using file_duration" '
  mpirun -n 1 ../test_diag_yaml
'

sed 's/file_duration: 12 hours/file_duration: 12 mullions/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "file_duration_units is not valid" '
  mpirun -n 1 ../test_diag_yaml
'

sed 's/freq: 6 hours/freq: 6 mullions/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "freq units is not valid" '
  mpirun -n 1 ../test_diag_yaml
'

sed 's/freq: 6 hours/freq: -6 hours/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "freq is less than -1" '
  mpirun -n 1 ../test_diag_yaml
'

sed 's/kind: r4/kind: mullions/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "kind is not valid" '
  mpirun -n 1 ../test_diag_yaml
'

sed 's/reduction: average/reduction: mullions/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "reduction is not valid" '
  mpirun -n 1 ../test_diag_yaml
'

sed 's/reduction: average/reduction: diurnal0/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "diurnal samples is less than 0" '
  mpirun -n 1 ../test_diag_yaml
'

sed 's/reduction: average/reduction: diurnal99r/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "diurnal samples is not an integer" '
  mpirun -n 1 ../test_diag_yaml
'

sed 's/reduction: average/reduction: pow0/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "power value is less than 0" '
  mpirun -n 1 ../test_diag_yaml
'

sed 's/reduction: average/reduction: pow99r/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "power value is not an integer" '
  mpirun -n 1 ../test_diag_yaml
'

sed 's/grid_type: latlon/grid_type: ice_cream/g' diag_table.yaml_base > diag_table.yaml
test_expect_failure "the sub_region grid_type is not valid" '
  mpirun -n 1 ../test_diag_yaml
'
