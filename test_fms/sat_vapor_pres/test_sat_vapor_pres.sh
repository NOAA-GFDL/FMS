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
# execute tests in the test_fms/horiz_interp directory.

# Copyright 2021 Seth Underwood

# Set common test settings.
. ../test-lib.sh

# Prepare the directory to run the tests.
cat << EOF > input.nml
&sat_vapor_pres_nml
      construct_table_wrt_liq = .true.,
      construct_table_wrt_liq_and_ice = .true.,
      use_exact_qs = .true.
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

test_done
