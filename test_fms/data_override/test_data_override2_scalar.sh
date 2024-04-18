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

if [ ! -z $parser_skip ]; then
  cat <<_EOF > input.nml
&data_override_nml
use_data_table_yaml=.False.
/
&test_data_override_ongrid_nml
  test_case = 3
/
_EOF
  printf '"OCN", "co2", "co2", "./INPUT/scalar.nc", "none" ,  1.0' | cat > data_table
else
cat <<_EOF > input.nml
&data_override_nml
use_data_table_yaml=.True.
/
&test_data_override_ongrid_nml
  test_case = 3
/
_EOF
cat <<_EOF > data_table.yaml
data_table:
 - gridname          : OCN
   fieldname_code    : co2
   fieldname_file    : co2
   file_name         : INPUT/scalar.nc
   interpol_method   : none
   factor            : 1.0
_EOF
fi

[ ! -d "INPUT" ] && mkdir -p "INPUT"
for KIND in r4 r8
do
rm -rf INPUT/*
test_expect_success "data_override scalar field (${KIND})" '
  mpirun -n 6 ../test_data_override_ongrid_${KIND}
'

done
rm -rf INPUT *.nc # remove any leftover files to reduce size

test_done