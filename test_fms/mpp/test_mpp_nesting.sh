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

# Ed Hartnett 11/29/19
# Ryan Mulhall 2/2021

# Set common test settings.
. ../test-lib.sh

# create and enter directory for in/output
output_dir

# Create input_base.nml for test input
cat <<_EOF > input_base.nml
&test_mpp_domains_nml
nx=64
ny=64
nz=10
stackmax=10000000
debug=.false.
mpes = 3
whalo = 2
ehalo = 2
shalo = 2
nhalo = 2
x_cyclic_offset = 3
y_cyclic_offset = -4
warn_level = "fatal"
wide_halo_x = 0
wide_halo_y = 0
nx_cubic = 20
ny_cubic = 20
num_fields = 4
do_sleep = .false.
num_iter = 1
! NEST inputs
test_nest = .false.
test_nest_regional = .false.
num_nest = 3
tile_coarse =    1,  3,  7
tile_fine   =    7 , 8,  9
istart_coarse =  3,  3,  5
icount_coarse = 40,  5,  6
jstart_coarse =  3,  3,  6
jcount_coarse = 14,  6,  8
extra_halo = 0
ntiles_nest_all = 9
npes_nest_tile = 2, 2, 2, 2, 2, 2, 2, 1, 1
nest_level = 1, 1, 2
refine_ratio = 2, 2, 2
cyclic_nest = 'N'
! NEST inputs end
mix_2D_3D = .false.
ensemble_size = 1
layout_cubic = 0,0
layout_ensemble = 0,0
nthreads = 1
layout_tripolar = 0,0
wide_halo = .false.
/
_EOF

sed "s/test_nest = .false./test_nest = .true./" input_base.nml > input.nml
test_expect_success "update nest domain" '
    mpirun -n 16 ../test_mpp_nesting
'
## single face nest needs additional namelist changes
sed "s/tile_coarse =    1,  3,  7/tile_coarse =    1,  1,  2/" input_base.nml > input.nml
sed -i "s/tile_fine   =    7 , 8,  9/tile_fine   =    2 , 3,  4/" input.nml
sed -i "s/istart_coarse =  3,  3,  5/istart_coarse =  4,  3,  5/" input.nml
sed -i "s/icount_coarse = 40,  5,  6/icount_coarse = 12,  5,  6/" input.nml
sed -i "s/jstart_coarse =  3,  3,  6/jstart_coarse =  4,  3,  6/" input.nml
sed -i "s/jcount_coarse = 14,  6,  8/jcount_coarse = 12,  6,  8/" input.nml
sed -i "s/ntiles_nest_all = 9/ntiles_nest_all = 4/" input.nml
sed -i "s/npes_nest_tile = 2, 2, 2, 2, 2, 2, 2, 1, 1/npes_nest_tile = 2, 2, 2, 1/" input.nml
test_expect_success "update single face nest domain" '
    mpirun -n 7 ../test_mpp_nesting
'
sed -i "s/test_nest_regional = .false./test_nest_regional = .true./" input.nml
test_expect_success "nested with regional top level domain" '
    mpirun -n 7 ../test_mpp_nesting
'

test_done
