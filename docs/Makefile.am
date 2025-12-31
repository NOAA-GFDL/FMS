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

# This is part of FMS. It creates the Makefile for the docs
# subdirectory.

# Ed Hartnett 6/25/20
# Ryan Mulhall 12/1/20

# Run doxygen.
all:
	doxygen Doxyfile
	cat doxywarn.log
# Confirm warning log file is empty.
#	[ ! -s doxywarn.log ]

check: all

EXTRA_DIST = Doxyfile.in grouping.h
