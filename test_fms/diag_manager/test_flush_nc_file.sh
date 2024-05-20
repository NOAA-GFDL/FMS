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

# Set common test settings.
. ../test-lib.sh

if [ -z "${skipflag}" ]; then
# create and enter directory for in/output files
output_dir

cat <<_EOF > diag_table.yaml
title: test_nc_flush
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_flush
  time_units: hours
  unlimdim: time
  freq: 1 hours
  varlist:
  - module: atmos
    var_name: ua
    reduction: average
    kind: r4
_EOF

# remove any existing files that would result in false passes during checks
rm -f *.nc

my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n flush_nc_files=.true. \n/" | cat > input.nml
test_expect_success "Running diag_manager with flush_nc_files = .true. (test $my_test_count)" '
  mpirun -n 1 ../test_var_masks
'
fi
test_done
