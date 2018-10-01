! -*-f90-*-*
!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

#ifndef __FMS_PLATFORM_
#define __FMS_PLATFORM_


!Set type kinds.
#ifdef PORTABLE_KINDS
use,intrinsic :: iso_fortran_env, only: real128
use,intrinsic :: iso_c_binding, only: c_double,c_float,c_int64_t, &
                                      c_int32_t,c_int16_t,c_intptr_t
#define QUAD_KIND real128
#define DOUBLE_KIND c_double
#define FLOAT_KIND c_float
#define LONG_KIND c_int64_t
#define INT_KIND c_int32_t
#define SHORT_KIND c_int16_t
#define POINTER_KIND c_intptr_t
#else
!These values are not necessarily portable.
#define QUAD_KIND 16
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#define SHORT_KIND 2
#define POINTER_KIND 8
!DEC$ MESSAGE:'Using 8-byte addressing'
#endif


!Control "pure" functions.
#ifdef NO_F95
#define _PURE
!DEC$ MESSAGE:'Not using pure routines.'
#else
#define _PURE pure
!DEC$ MESSAGE:'Using pure routines.'
#endif


!Control array members of derived types.
#ifdef NO_F2000
#define _ALLOCATABLE pointer
#define _NULL =>null()
#define _ALLOCATED associated
!DEC$ MESSAGE:'Using pointer derived type array members.'
#else
#define _ALLOCATABLE allocatable
#define _NULL
#define _ALLOCATED allocated
!DEC$ MESSAGE:'Using allocatable derived type array members.'
#endif


!Control use of cray pointers.
#ifdef NO_CRAY_POINTERS
#undef use_CRI_pointers
!DEC$ MESSAGE:'Not using cray pointers.'
#else
#define use_CRI_pointers
!DEC$ MESSAGE:'Using cray pointers.'
#endif


!Control size of integers that will hold address values.
!Appears for legacy reasons, but seems rather dangerous.
#ifdef _32bits
#define POINTER_KIND 4
!DEC$ MESSAGE:'Using 4-byte addressing'
#endif


!If you do not want to use 64-bit integers.
#ifdef no_8byte_integers
#define LONG_KIND INT_KIND
#endif


!If you do not want to use 32-bit floats.
#ifdef no_4byte_reals
#define FLOAT_KIND DOUBLE_KIND
#define NF_GET_VAR_REAL nf_get_var_double
#define NF_GET_VARA_REAL nf_get_vara_double
#define NF_GET_ATT_REAL nf_get_att_double
#undef OVERLOAD_R4
#undef OVERLOAD_C4
#endif


!If you want to use quad-precision.
#define NO_QUAD_PRECISION
#ifdef QUAD_PRECISION
#undef NO_QUAD_PRECISION
#else
#define QUAD_KIND DOUBLE_KIND
#endif


#endif
