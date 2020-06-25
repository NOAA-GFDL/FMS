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

# This is part of FMS. It creates the doc Makefile.

# Ed Hartnett 6/25/20

# Run doxygen.
all:
	doxygen Doxyfile
	cat doxywarn.log
# Confirm warning log file is empty.
#	[ ! -s doxywarn.log ]

check: all

# Include these subdirs to include the documention files in the
# distribution.
#SUBDIRS = source images

CLEANFILES = *.log

EXTRA_DIST = Doxyfile.in
