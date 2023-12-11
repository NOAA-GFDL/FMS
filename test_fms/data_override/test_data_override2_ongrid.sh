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
  cat <<_EOF > input_base.nml
&data_override_nml
use_data_table_yaml=.False.
/

&test_data_override_ongrid_nml
  nhalox=halo_size
  nhaloy=halo_size
/
_EOF
  printf '"OCN", "runoff", "runoff", "./INPUT/runoff.daitren.clim.1440x1080.v20180328.nc", "none" ,  1.0' | cat > data_table
else
cat <<_EOF > input_base.nml
&data_override_nml
use_data_table_yaml=.True.
/

&test_data_override_ongrid_nml
  nhalox=halo_size
  nhaloy=halo_size
/
_EOF
cat <<_EOF > data_table.yaml
data_table:
 - gridname          : OCN
   fieldname_code    : runoff
   fieldname_file    : runoff
   file_name         : INPUT/runoff.daitren.clim.1440x1080.v20180328.nc
   interpol_method   : none
   factor            : 1.0
_EOF
fi

[ ! -d "INPUT" ] && mkdir -p "INPUT"
for KIND in r4 r8
do
rm -rf INPUT/*
sed 's/halo_size/2/g' input_base.nml > input.nml
test_expect_success "data_override on grid with 2 halos in x and y (${KIND})" '
  mpirun -n 6 ../test_data_override_ongrid_${KIND}
'

rm -rf INPUT/*
sed 's/halo_size/0/g' input_base.nml > input.nml
test_expect_success "data_override on grid with 0 halos in x and y (${KIND})" '
  mpirun -n 6 ../test_data_override_ongrid_${KIND}
'

rm -rf INPUT/*
test_expect_success "data_override get_grid_v1 (${KIND})" '
  mpirun -n 1 ../test_get_grid_v1_${KIND}
'
done
rm -rf INPUT *.nc # remove any leftover files to reduce size

test_done