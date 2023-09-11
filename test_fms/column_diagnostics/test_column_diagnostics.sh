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
EOF


#####
test_expect_success "test_column_diagnostics r4" 'mpirun -n 1 ./test_column_diagnostics_r4'
test_expect_success "test_column_diagnostics r8" 'mpirun -n 1 ./test_column_diagnostics_r8'
test_done
