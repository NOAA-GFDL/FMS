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
#
# Copyright (c) 2021 Seth Underwood

. ../test-lib.sh

# Copy and rename namelist file.
cat <<_EOF > input.nml
&test_affinity_nml
/
_EOF

if test "x$mpi_launcher" != "xsrun"; then
  SKIP_TESTS="$SKIP_TESTS $(basename $0 .sh).1"
fi

test_expect_success "FMS affinity places MPI processes correctly on slurm" '
  mpirun -n 6 ./test_affinity
'

test_done
