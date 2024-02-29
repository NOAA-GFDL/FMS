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

# data_override with the default table (not setting namelist)
rm -rf input.nml data_table data_table.yaml
touch input.nml
touch data_table
test_expect_success "data_override_init with the default table" '
 mpirun -n 1 ../test_data_override_init
'

cat <<_EOF > input.nml
&data_override_nml
use_data_table_yaml=.false.
/
_EOF
test_expect_success "data_override_init setting use_data_table_yaml = .false." '
 mpirun -n 1 ../test_data_override_init
'

touch data_table.yaml
test_expect_failure "data_override_init both tables present" '
 mpirun -n 1 ../test_data_override_init
'
if [ ! -z $parser_skip ]; then
rm -rf data_table.yaml
cat <<_EOF > input.nml
&data_override_nml
use_data_table_yaml=.true.
/
_EOF
  test_expect_failure "data_override_init setting use_data_table_yaml = .true. but no compiling with yaml" '
    mpirun -n 1 ../test_data_override_init
  '
else
rm -rf data_table
cat <<_EOF > input.nml
&data_override_nml
use_data_table_yaml=.true.
/
_EOF
test_expect_success "data_override_init setting use_data_table_yaml = .true." '
 mpirun -n 1 ../test_data_override_init
'
fi
test_done