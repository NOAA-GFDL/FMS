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

# Copyright (c) 2019-2020 Ed Hartnett, Seth Underwood

# Set common test settings.
. ../test-lib.sh

if [ -z "${skipflag}" ]; then
# create and enter directory for in/output files
output_dir

cat <<_EOF > diag_table.ens_01.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_ens
  time_units: days
  unlimdim: time
  freq: 1 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: none
    kind: r8
_EOF

cat <<_EOF > diag_table.ens_02.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_ens
  time_units: days
  unlimdim: time
  freq: 2 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: none
    kind: r8
_EOF

cat <<_EOF > input.nml
&diag_manager_nml
  use_modern_diag = .True.
/

&ensemble_nml
   ensemble_size = 2
/
_EOF

my_test_count=1
test_expect_success "Running diag_manager with 2 ensembles (test $my_test_count)" '
  mpirun -n 2 ../test_ens_runs
'

cat <<_EOF > diag_table.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_0days
  time_units: days
  unlimdim: time
  freq: 0 days
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: none
    kind: r8
_EOF

my_test_count=`expr $my_test_count + 1`
test_expect_failure "Running diag_manager with both diag_table.yaml and diag_table.ens_xx.yaml files present (test $my_test_count)" '
  mpirun -n 2 ../test_ens_runs
'

fi
test_done
