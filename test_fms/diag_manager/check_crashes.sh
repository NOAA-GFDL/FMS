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
. ../test_common.sh

printf "&check_crashes_nml \n checking_crashes = .true. \n/" | cat > input.nml

echo "Test 27: Missing tile when using the 'index' grid type"
touch input.nml
sed '/tile/d' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since 'tile' was missing and the 'index' grid type was used"
  exit 3
fi

echo "Test 28: Missing new_file_freq_units when using new_file_freq_units"
touch input.nml
sed '/new_file_freq_units/d' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since 'new_file_freq_units' was missing and new_file_freq was used"
  exit 3
fi

echo "Test 29: new_file_freq_units is not valid"
touch input.nml
sed 's/new_file_freq_units: hours/new_file_freq_units: mullions/g' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since 'new_file_freq_units' is not valid"
  exit 3
fi

echo "Test 30: Missing file_duration_units when using file_duration"
touch input.nml
sed '/file_duration_units/d' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since 'file_duration_units' was missing and file_duration was used"
  exit 3
fi

echo "Test 31: file_duration_units is not valid"
touch input.nml
sed 's/file_duration_units: hours/file_duration_units: mullions/g' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since 'file_duration_units' is not valid"
  exit 3
fi

echo "Test 32: freq units is not valid"
touch input.nml
sed 's/freq_units: hours/freq_units: mullions/g' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since the freq units is not valid"
  exit 3
fi

echo "Test 33: freq is less than 0"
touch input.nml
sed 's/freq: 6/freq: -666/g' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since freq is not valid"
  exit 3
fi

echo "Test 34: realm is not valid"
touch input.nml
sed 's/realm: ATM/realm: mullions/g' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since realm is not valid"
  exit 3
fi

echo "Test 35: kind is not valid"
touch input.nml
sed 's/kind: float/kind: mullions/g' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since the kind is not valid"
  exit 3
fi

echo "Test 36: reduction is not valid"
touch input.nml
sed 's/reduction: average/reduction: mullions/g' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since the reduction method is not valid"
  exit 3
fi

echo "Test 37: diurnal samples is less than 0"
touch input.nml
sed 's/reduction: average/reduction: diurnal0/g' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since the number of diurnal samples is less than 0"
  exit 3
fi

echo "Test 38: diurnal samples is not an integer"
touch input.nml
sed 's/reduction: average/reduction: diurnal99r/g' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since the number of diurnal samples is not valid"
  exit 3
fi

echo "Test 39: power value is less than 0"
touch input.nml
sed 's/reduction: average/reduction: pow0/g' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since the power value is less than"
  exit 3
fi

echo "Test 40: power value is not an integer"
touch input.nml
sed 's/reduction: average/reduction: pow99r/g' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since the power value is not valid"
  exit 3
fi

echo "Test 41: the sub_region grid_type is not valid"
touch input.nml
sed 's/grid_type: latlon/grid_type: ice_cream/g' $top_srcdir/test_fms/diag_manager/diagTables/diag_table_yaml_26 > diag_table.yaml
run_test test_diag_yaml 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "Test should have failed since the sub_region grid_type"
  exit 3
fi
