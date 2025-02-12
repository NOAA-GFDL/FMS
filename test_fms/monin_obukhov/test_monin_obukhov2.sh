#!/bin/sh

#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
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
