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
# execute tests in the test_fms/horiz_interp directory.

# Copyright 2021 Seth Underwood

# Set common test settings.
. ../test-lib.sh

# Prepare the directory to run the tests.
cat << EOF > input.nml
&sat_vapor_pres_nml
      construct_table_wrt_liq = .true.,
      construct_table_wrt_liq_and_ice = .true.,
      use_exact_qs = .true.,
      show_all_bad_values = .true.
/
EOF


#####
cat <<EOF > test_sat_vapor_pres.nml
&test_sat_vapor_pres_nml
  test1=.true.
  test2=.false.
  test3=.false.
  test4=.false.
  test5=.false.
 /
EOF
test_expect_success "test_compute_qs_r4" '
      mpirun -n 1 ./test_sat_vapor_pres_r4
  '
test_expect_success "test_compute_qs_r8" '
      mpirun -n 1 ./test_sat_vapor_pres_r8
  '

####
cat <<EOF > test_sat_vapor_pres.nml
&test_sat_vapor_pres_nml
  test1=.false.
  test2=.true.
  test3=.false.
  test4=.false.
  test5=.false.
 /
EOF
test_expect_success "test_compute_mrs_r4" '
      mpirun -n 1 ./test_sat_vapor_pres_r4
  '
test_expect_success "test_compute_mrs_r8" '
      mpirun -n 1 ./test_sat_vapor_pres_r8
  '

####
cat <<EOF > test_sat_vapor_pres.nml
&test_sat_vapor_pres_nml
  test1=.false.
  test2=.false.
  test3=.true.
  test4=.false.
  test5=.false.
 /
EOF
test_expect_success "test_lookup_es_des_r4" '
      mpirun -n 1 ./test_sat_vapor_pres_r4
  '
test_expect_success "test_lookup_es_des_r8" '
      mpirun -n 1 ./test_sat_vapor_pres_r8
  '

####
cat <<EOF > test_sat_vapor_pres.nml
&test_sat_vapor_pres_nml
  test1=.false.
  test2=.false.
  test3=.false.
  test4=.true.
  test5=.false.
 /
EOF
test_expect_success "test_lookup_es2_des2_r4" '
      mpirun -n 1 ./test_sat_vapor_pres_r4
  '
test_expect_success "test_lookup_es2_des2_r8" '
      mpirun -n 1 ./test_sat_vapor_pres_r8
  '

####
cat <<EOF > test_sat_vapor_pres.nml
&test_sat_vapor_pres_nml
  test1=.false.
  test2=.false.
  test3=.false.
  test4=.false.
  test5=.true.
 /
EOF
test_expect_success "test_lookup_es3_des3_r4" '
      mpirun -n 1 ./test_sat_vapor_pres_r4
  '
test_expect_success "test_lookup_es3_des3_r8" '
      mpirun -n 1 ./test_sat_vapor_pres_r8
  '

## test failures when out of range temps are used
cat <<EOF > test_sat_vapor_pres.nml
&test_sat_vapor_pres_nml
  test1=.false.
  test2=.false.
  test3=.true.
  test4=.false.
  test5=.false.
  test_show_all_bad = 0
 /
EOF

test_expect_failure "check bad temperature values 0d r4" '
      mpirun -n 2 ./test_sat_vapor_pres_r4
  '
test_expect_failure "check bad temperature values 0d r8" '
      mpirun -n 2 ./test_sat_vapor_pres_r8
  '

sed -i 's/test_show_all_bad = 0/test_show_all_bad = 1/' test_sat_vapor_pres.nml

test_expect_failure "check bad temperature values 1d r4" '
      mpirun -n 2 ./test_sat_vapor_pres_r4
  '
test_expect_failure "check bad temperature values 1d r8" '
      mpirun -n 2 ./test_sat_vapor_pres_r8
  '

sed -i 's/test_show_all_bad = 1/test_show_all_bad = 2/' test_sat_vapor_pres.nml

test_expect_failure "check bad temperature values 2d r4" '
      mpirun -n 2 ./test_sat_vapor_pres_r4
  '
test_expect_failure "check bad temperature values 2d r8" '
      mpirun -n 2 ./test_sat_vapor_pres_r8
  '

sed -i 's/test_show_all_bad = 2/test_show_all_bad = 3/' test_sat_vapor_pres.nml

test_expect_failure "check bad temperature values 3d r4" '
      mpirun -n 2 ./test_sat_vapor_pres_r4
  '
test_expect_failure "check bad temperature values 3d r8" '
      mpirun -n 2 ./test_sat_vapor_pres_r8
  '

test_done
