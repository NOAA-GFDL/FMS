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

cat <<EOF > input.nml
&diag_integral_nml
      file_name  = 'diag_integral.out',
      time_units = 'days',
      output_interval = 1.0
      fields_per_print_line = 6
/
EOF

mkdir -p INPUT

test_expect_success "test_diag_integral r4" 'mpirun -n 1 ./test_diag_integral_r4'
test_expect_success "test_diag_integral r8" 'mpirun -n 1 ./test_diag_integral_r8'

rm -rf INPUT

test_done
