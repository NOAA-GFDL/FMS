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

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/mpp directory.

# Ryan Mulhall 02/2021

# Set common test settings.
. ../test-lib.sh

output_dir

# ensure input.nml file present
cat <<_EOF  > input.nml
&test_stdlog_nml
  test_num = 1
/
_EOF
# Run test with one processor
test_expect_success "test stdlog and stdwarn" '
    mpirun -n 2 ../test_stdlog
'
sed -i 's/1/2/' input.nml
test_expect_failure "test stdlog and stdwarn with fatal output" '
    mpirun -n 2 ../test_stdlog
'
# move file so we don't overwrite
mv warnfile.*.out warnfile.000000.out.old
sed -i 's/2/3/' input.nml
test_expect_success "check stdwarn output" '
    mpirun -n 1 ../test_stdlog
'
test_done
