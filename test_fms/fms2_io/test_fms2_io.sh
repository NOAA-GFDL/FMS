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
# execute tests in the test_fms/fms2_io directory.

# Authors: Raymond Menzel
# Jessica Liptak
#
# Set common test settings.
. ../test-lib.sh

# Create and enter output directory
output_dir

touch input.nml

# run the tests
test_expect_success "FMS2 IO Test" '
  mpirun -n 6 ../test_fms2_io
'

test_expect_success "Compressed writes tests" '
  mpirun -n 5 ../test_compressed_writes
'

cat <<_EOF > input.nml
&test_domain_io_nml
  layout = 1, 6
  io_layout = 1, 1
  filename = "test_simple_layout.nc"
/
_EOF
test_expect_success "Domain Read Write Tests with simple layout" '
  mpirun -n 6 ../test_domain_io
'

cat <<_EOF > input.nml
&test_domain_io_nml
  layout = 2, 8
  io_layout = 1, 2
  filename = "test_dist_layout.nc"
/
_EOF
test_expect_success "Domain Read Write Tests with 2 distributed files" '
  mpirun -n 16 ../test_domain_io
'

cat <<_EOF > input.nml
&test_domain_io_nml
  layout = 2, 8
  io_layout = 1, 2
  filename = "test_dist_layout.nc"
  use_edges = .true.
/
_EOF
test_expect_success "Domain Read Write Tests with 2 distributed files and EAST and NORTH axis" '
  mpirun -n 16 ../test_domain_io
'

cat <<_EOF > input.nml
&test_domain_io_nml
  nx = 33
  ny = 43
  layout = 4, 6
  io_layout = 2, 3
  filename = "test_non_uniform.nc"
/
_EOF
test_expect_success "Domain Read Write Tests with non uniform layouts" '
  mpirun -n 24 ../test_domain_io
'

cat <<_EOF > input.nml
&test_domain_io_nml
  layout = 3, 6
  io_layout = 1, 2
  mask_table = "mask_table"
  filename = "test_io_mask.nc"
/
_EOF

cat <<_EOF > mask_table
1
3,6
1,1
_EOF
test_expect_success "Domain Read Write Tests with a ocean mask" '
  mpirun -n 17 ../test_domain_io
'

rm input.nml
touch input.nml
test_expect_success "Write and read domain-decomposed data with generalized indices" '
  mpirun -n 16 ../test_generalized_indices
'

test_done
