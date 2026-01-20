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

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/horiz_interp directory.

# Copyright 2021 Seth Underwood

# Set common test settings.
. ../test-lib.sh

cat <<EOF > input.nml
&diag_integral_nml
      file_name  = 'diag_integral.out',
      time_units = 'days',
      output_interval = 1.0
      fields_per_print_line = 6
/
EOF

mkdir -p INPUT

test_expect_success "test_diag_integral r4" 'mpirun -n 1 ./test_diag_integral_r4'
test_expect_success "test_diag_integral r8" 'mpirun -n 1 ./test_diag_integral_r8'

rm -rf INPUT

test_done
