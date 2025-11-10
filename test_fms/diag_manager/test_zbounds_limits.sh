#!/bin/sh

#***********************************************************************
#*                             Apache License 2.0
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* Licensed under the Apache License, Version 2.0 (the "License");
#* you may not use this file except in compliance with the License.
#* You may obtain a copy of the License at
#*
#*     http://www.apache.org/licenses/LICENSE-2.0
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
#* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#* PARTICULAR PURPOSE. See the License for the specific language
#* governing permissions and limitations under the License.
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
test_expect_success "Test zbounds limits (legacy diag manager)" '
  mpirun -n 6 ../test_zbounds_limits
'

if [ -z "${parser_skip}" ]
then
  # Repeat the test with the modern diag_manager
  cat <<_EOF > input.nml
&diag_manager_nml
  use_modern_diag = .true.
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
  test_expect_success "Test zbounds limits (modern diag manager)" '
    mpirun -n 6 ../test_zbounds_limits
  '
fi

test_done
