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
title: test_prepend_date
base_date: 1 1 1 0 0 0
diag_files:
- file_name: test_non_static
  time_units: hours
  unlimdim: time
  freq: 1 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: average
    kind: r4
- file_name: test_static
  time_units: hours
  unlimdim: time
  freq: -1 hours
  varlist:
  - module: ocn_mod
    var_name: var2
    reduction: none
    kind: r4
_EOF

# remove any existing files that would result in false passes during checks
rm -f *.nc
my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n/" | cat > input.nml
test_expect_success "Running diag_manager and checking that the date was prepended correctly (test $my_test_count)" '
  mpirun -n 1 ../test_prepend_date
'

cat <<_EOF > diag_table.yaml
title: test_prepend_date
base_date: 1 1 1 0 0 0
diag_files:
- file_name: test_non_static
  time_units: hours
  unlimdim: time
  freq: 1 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: average
    kind: r4
  - module: ocn_mod
    var_name: var1
    reduction: average
    kind: r4
_EOF

printf "&diag_manager_nml \n use_modern_diag=.true. \n/ \n &test_prepend_date_nml \n pass_diag_time=.false. \n /" | cat > input.nml

test_expect_failure "Running diag_manager with fields that have a different start time (test $my_test_count)" '
  mpirun -n 1 ../test_prepend_date
'

fi
test_done
