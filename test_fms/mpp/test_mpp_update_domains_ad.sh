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

# Jessica Liptak

# Set common test settings.
. ../test-lib.sh
# Run the test for one processor
touch input.nml

test_expect_success "update adjoint domains with 1 PE" '
    mpirun -n 1 ./test_mpp_update_domains_ad
'
test_expect_success "update adjoint domains with 4 PEs" '
    mpirun -n 4 ./test_mpp_update_domains_ad
'
test_done
