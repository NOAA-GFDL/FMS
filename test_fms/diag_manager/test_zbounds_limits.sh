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

output_dir

touch input.nml
cat <<_EOF > diag_table
test_diag_manager_01
2 1 1 0 0 0

"file_1", 1, "hours", 1, "hours", "time"
 "lnd_mod", "var1", "var1", "file_1", "all", "max", "-1 -1 -1 -1 2 3", 2
_EOF

rm -f *.nc
my_test_count=1
test_expect_success "Running diag_manager (test $my_test_count)" '
  mpirun -n 6 ../test_zbounds_limits
'

# Repeat the test with the modern diag_manager
cat <<_EOF > input.nml
& diag_manager_nml
  use_modern_diag = .True.
/
_EOF

cat <<_EOF > diag_table.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: file_1
  time_units: hours
  unlimdim: time
  freq: 1 hours
  varlist:
  - module: lnd_mod
    var_name: var1
    reduction: max
    zbounds: 2 3
    kind: r4
_EOF

rm -rf diag_table
rm -f *.nc
my_test_count=2
test_expect_success "Running with modern diag manager (test $my_test_count)" '
  mpirun -n 6 ../test_zbounds_limits
'
test_done
