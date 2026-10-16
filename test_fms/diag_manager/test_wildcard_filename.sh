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

# Tests diag_manager_nml::wildcard_filename_prefix/wildcard_filename_separator, which control how
# the time fields substituted into a wildcard (%) output file name are joined together.

# Set common test settings.
. ../test-lib.sh

if [ -z "${parser_skip}" ]; then
# create and enter directory for in/output files
output_dir

cat <<_EOF > diag_table.yaml
title: test_wildcard_filename
base_date: 2 1 1 0 0 0
diag_files:
- file_name: test_wildcard%4yr%2mo%2dy%2hr
  freq: 6 hours
  time_units: hours
  unlimdim: time
  new_file_freq: 12 hours
  file_duration: 18 hours
  varlist:
  - module: ocn_mod
    var_name: var0
    reduction: none
    kind: r4
_EOF

# remove any existing files that would result in false passes during checks
rm -f *.nc

my_test_count=1
printf "&diag_manager_nml \n use_modern_diag=.true. \n/" | cat > input.nml
test_expect_success "Default wildcard_filename_prefix/separator reproduce the historical '_' joined suffix (test $my_test_count)" '
  mpirun -n 1 ../test_wildcard_filename
'

rm -f *.nc
my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n wildcard_filename_prefix='.' \n wildcard_filename_separator='-' \n/ \n &test_wildcard_filename_nml \n expected_suffix1='.0002-01-01-06' \n expected_suffix2='.0002-01-01-15' \n/" | cat > input.nml
test_expect_success "Custom wildcard_filename_prefix/separator are used to join the substituted time fields (test $my_test_count)" '
  mpirun -n 1 ../test_wildcard_filename
'

rm -f *.nc
my_test_count=`expr $my_test_count + 1`
printf "&diag_manager_nml \n use_modern_diag=.true. \n wildcard_filename_separator='/' \n/" | cat > input.nml
test_expect_failure "wildcard_filename_separator containing '/' is rejected (test $my_test_count)" '
  mpirun -n 1 ../test_wildcard_filename
'

fi
test_done
