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
# execute tests in the test_fms/monin_obukhov directory.

# Set common test settings.
. ../test-lib.sh

# Skipping these tests for now, to avoid CI failure due to constants_mod values
# which change according to the default real kind.
# TODO: Enable these tests after constants_mod and/or the CI has been updated
SKIP_TESTS="test_monin_obukhov2.1 test_monin_obukhov2.2"

# Run tests
for p in r4 r8
do
  cp ${top_srcdir}/test_fms/monin_obukhov/input.${p}.nml input.nml
  test_expect_success "test monin_obukhov_mod (${p})" "mpirun -n 1 ./test_monin_obukhov_${p}"
  rm -f input.nml
done

test_done
