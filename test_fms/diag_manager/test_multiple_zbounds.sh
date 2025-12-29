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

if [ -z "${parser_skip}" ]
then
  cat <<_EOF > input.nml
&diag_manager_nml
  use_modern_diag = .true.
/
_EOF

  cat <<_EOF > diag_table.yaml
title: test_multiple_zbounds
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_multiple_zbounds
  freq: 1 hours
  time_units: hours
  unlimdim: time
  module: atmos
  reduction: average
  kind: r4
  varlist:
  - var_name: ua_1
    zbounds: 3 5
  - var_name: ua_2
    zbounds: 1 1
_EOF

  test_expect_success "Test with multiple zbounds limits (modern diag manager)" '
    mpirun -n 1 ../test_multiple_zbounds
  '
fi

test_done
