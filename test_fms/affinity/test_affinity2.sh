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
