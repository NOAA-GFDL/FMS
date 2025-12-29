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

# Tom Robinson 04/21/2020
# Ryan Mulhall 2/2021

# Set common test settings.
. ../test-lib.sh

# ensure input.nml file present
touch input.nml

export NUM_PES=1
test_expect_success "get number of PEs single processor" '
    mpirun -n 1 ./test_mpp_npes
'
export NUM_PES=2
test_expect_success "get number of PEs multiple processor" '
    mpirun -n 2 ./test_mpp_npes
'
test_done
