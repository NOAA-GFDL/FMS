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
# execute tests in the test_fms/fms directory.

# Tom Robinson 03/02/2021

# Set common test settings.
. ../test-lib.sh

setup_test () {
  rm -rf RESTART && mkdir -p RESTART

# Create the base input.nml file needed for the tests
cat <<_EOF > input.nml
&test_fms_nml
/
_EOF
}

# Test the structured grid
setup_test
test_expect_success "test structured grid" '
  mpirun -n 6 ./test_fms
'
test_done
