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

# This is an automake file for the test_fms/fms directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the fms_mod.mod file.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/test_fms/fms/include -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_fms

# These are the sources for the tests.
test_fms_SOURCES = \
  test_fmsC.c \
  test_fms.F90 \
  include/test_fms.inc \
  include/test_fms_r4.fh \
  include/test_fms_r8.fh

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh

# Run the test programs.
TESTS = test_fms2.sh

# These will also be included in the distribution.
EXTRA_DIST = test_fms2.sh

CLEANFILES = input.nml *.out *.mod *.o *.dpi *.spi *.dyn *.spl

clean-local:
	rm -rf RESTART
