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

# TODO: Enable these tests once generalized indices work is complete
SKIP_TESTS="test_mpp_global_field.1 \
            test_mpp_global_field.2 \
            test_mpp_global_field.3 \
            test_mpp_global_field.4"

touch input.nml

for datatype in r4 r8 i4 i8
do
  test_expect_success "mpp global field functions ($datatype)" "
      mpirun -n 4 ./test_mpp_global_field_$datatype
  "
done

test_done
