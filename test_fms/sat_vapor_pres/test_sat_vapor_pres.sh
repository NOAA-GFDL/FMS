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

cat <<EOF > test_sat_vapor_pres_in.nml
&test_sat_vapor_pres_nml
  test1=.false.
  test2=.false.
  test3=.false.
  test4=.false.
 /
EOF

# selects next test to run through nml
testNum=0
test_next()
{
  testNum=$((testNum + 1))
  sed "s/test$testNum *=.false./test$testNum =.true./" test_sat_vapor_pres_in.nml > test_sat_vapor_pres.nml
  # test #8 must set calendar type for #9 to pass
  test_expect_success "$1" '
      mpirun -n 1 ./test_sat_vapor_pres
  '
}

test_next "test_compute_qs"
test_next "test_mrs"
test_next "test_lookup_es_des"
test_next "test_lookup_es2_des2"
test_next "test_lookup_es3_des3"

test_done
