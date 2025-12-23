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

if [ -z "${parser_skip}" ]; then
# create and enter directory for in/output files
output_dir

###########################################################################
# CASE 1: In this example, a file is going to be created every 12 hours, it is going to use the
# end of the time_bounds (i.e 2 1 1 0 12 0 and 2 1 1 1 0 0) for the filename
cat <<_EOF > diag_table.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_multi_file%4yr%2mo%2dy%2hr
  time_units: days
  unlimdim: time
  freq: 2 hours
  new_file_freq: 12 hours
  module: atmos
  reduction: average
  kind: r4
  filename_time: end
  varlist:
  - var_name: ua
_EOF

# remove any existing files that would result in false passes during checks
rm -f *.nc
my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n/" | cat > input.nml
test_expect_success "Running diag_manager (test $my_test_count)" '
  mpirun -n 1 ../test_diag_multi_file
'

###########################################################################
# CASE 2: In this example, a file is going to be created starting from [2 1 1 12 0 0],
# every 4 hours, for 8 hours (so 2 files!)
rm -f *.nc
cat <<_EOF > diag_table.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_multi_file%4yr%2mo%2dy%2hr
  time_units: days
  unlimdim: time
  freq: 2 hours
  start_time: 2 1 1 12 0 0
  new_file_freq: 4 hours
  file_duration: 8 hours
  module: atmos
  reduction: average
  kind: r4
  filename_time: end
  varlist:
  - var_name: ua
_EOF

cat <<_EOF > input.nml
&diag_manager_nml
 use_modern_diag=.true.
/
&test_multi_file_nml
 test_case = 2
/
_EOF
test_expect_success "Running diag_manager (test $my_test_count)" '
  mpirun -n 1 ../test_diag_multi_file
'

###########################################################################
# CASE 3: In this example, a file is going to be created using data from [2 1 12 0 0]
# to [2 1 1 20 0 0] (8 hours!)
rm -f *.nc
cat <<_EOF > diag_table.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_multi_file
  time_units: days
  unlimdim: time
  freq: 2 hours
  start_time: 2 1 1 12 0 0
  file_duration: 8 hours
  module: atmos
  reduction: average
  kind: r4
  varlist:
  - var_name: ua
_EOF

cat <<_EOF > input.nml
&diag_manager_nml
 use_modern_diag=.true.
/
&test_multi_file_nml
 test_case = 3
/
_EOF
test_expect_success "Running diag_manager (test $my_test_count)" '
  mpirun -n 1 ../test_diag_multi_file
'

###########################################################################
# CASE 4: In this example, a file is going to be created every 4 hours for 8 hours
# (so 2 files), then every 2 hours for 2 hours (so 1 file), then every 12 hours for 12 hours
# (so 1 file)
rm -f *.nc
cat <<_EOF > diag_table.yaml
title: test_diag_manager_01
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_multi_file%4yr%2mo%2dy%2hr
  time_units: days
  unlimdim: time
  freq: 2 hours, 2 hours, 2 hours
  new_file_freq: 4 hours, 2 hours, 12 hours
  file_duration: 8 hours, 2 hours, 12 hours
  filename_time: end
  module: atmos
  reduction: average
  kind: r4
  varlist:
  - var_name: ua
_EOF

cat <<_EOF > input.nml
&diag_manager_nml
 use_modern_diag=.true.
/
&test_multi_file_nml
 test_case = 4
/
_EOF
test_expect_success "Running diag_manager (test $my_test_count)" '
  mpirun -n 1 ../test_diag_multi_file
'
fi
test_done
