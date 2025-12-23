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

# Ryan Mulhall 2/2021

# Set common test settings.
. ../test-lib.sh

echo "&test_mpp_chksum_nml" > input.nml
echo "test_num = 1" >> input.nml
echo "/" >> input.nml

test_expect_success "mpp_chksum simple functionality" '
    mpirun -n 4 ./test_mpp_chksum
'

sed -i 's/test_num = 1/test_num = 2/' input.nml
test_expect_success "mpp integer checksums" '
    mpirun -n 4 ./test_mpp_chksum
'

sed -i 's/test_num = 2/test_num = 3/' input.nml
test_expect_success "mpp_chksum with mixed precision" '
    mpirun -n 4 ./test_mpp_chksum
'
test_done
