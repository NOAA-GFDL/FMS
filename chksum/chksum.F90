!-----------------------------------------------------------------------
!         Checksums for bitwise comparison of floating point data
!
! AUTHOR: V. Balaji (vb@gfdl.gov)
!         SGI/GFDL Princeton University
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.  
!-----------------------------------------------------------------------

!these are used to determine hardware/OS/compiler

#ifdef __sgi
#ifdef _COMPILER_VERSION
#define sgi_mipspro
#else
#define sgi_generic
#endif
#endif

#if defined(_CRAY) || defined(sgi_mipspro)
#define SGICRAY
#endif

!compilers that support Cray pointers
#if defined(SGICRAY) || defined(__alpha)
#define use_CRI_pointers
#endif

!values of kind: double and long are 8-byte, float and int are 4-byte
#if defined(SGICRAY)
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#else
!these might be different on non-SGICRAY, I believe
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#endif
#ifdef sgi_generic
!this is for the Edinburgh n32/o32 compiler, which won't accept 8-byte ints
!at any price
#define LONG_KIND 4
#endif

module chksum_mod
!this module contains long integer function chksum.
!this does an exact sum of its argument as an integer
!(floating point checksums do not guarantee that all the bits match)
!chksum_mod uses Cray pointers to equivalence real/cmplx to int,
!then sums the resulting array as int.
!this works with IEEE words, where any bit sequence is a valid integer
  implicit none
  private
  character(len=256), private :: version='$Id: chksum.F90,v 1.2 1999/11/24 00:38:07 vb Exp $'

  interface chksum
     module procedure chksum_int_1d
     module procedure chksum_int_2d
     module procedure chksum_int_3d
     module procedure chksum_int_4d
#ifdef use_CRI_pointers
     module procedure chksum_r8_0d
     module procedure chksum_r8_1d
     module procedure chksum_r8_2d
     module procedure chksum_r8_3d
     module procedure chksum_r8_4d
     module procedure chksum_c8_0d
     module procedure chksum_c8_1d
     module procedure chksum_c8_2d
     module procedure chksum_c8_3d
     module procedure chksum_c8_4d
#endif
  end interface

  public :: chksum
  
  contains

    function chksum_int_1d(var)
      integer(LONG_KIND) :: chksum_int_1d
      integer(LONG_KIND), intent(in) :: var(:)
      chksum_int_1d = sum(var)
      return
    end function chksum_int_1d

    function chksum_int_2d(var)
      integer(LONG_KIND) :: chksum_int_2d
      integer(LONG_KIND), intent(in) :: var(:,:)
      chksum_int_2d = sum(var)
      return
    end function chksum_int_2d

    function chksum_int_3d(var)
      integer(LONG_KIND) :: chksum_int_3d
      integer(LONG_KIND), intent(in) :: var(:,:,:)
      chksum_int_3d = sum(var)
      return
    end function chksum_int_3d

    function chksum_int_4d(var)
      integer(LONG_KIND) :: chksum_int_4d
      integer(LONG_KIND), intent(in) :: var(:,:,:,:)
      chksum_int_4d = sum(var)
      return
    end function chksum_int_4d

#ifdef use_CRI_pointers
    function chksum_r8_0d(var)
      integer(LONG_KIND) :: chksum_r8_0d
      real(DOUBLE_KIND), intent(in) :: var
      integer(LONG_KIND) :: i(2)
      pointer( ptr, i )
      ptr = LOC(var)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          chksum_r8_0d = chksum_int_1d( i(1:2) )
      else
          chksum_r8_0d = chksum_int_1d( i(1:1) )
      endif
      return
    end function chksum_r8_0d

    function chksum_r8_1d(var)
      integer(LONG_KIND) :: chksum_r8_1d
      real(DOUBLE_KIND), intent(in) :: var(:)
      real(DOUBLE_KIND) :: g(size(var))
      integer(LONG_KIND) :: i(size(var)*2)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          chksum_r8_1d = chksum_int_1d( i(1:size(var)*2) )
      else
          chksum_r8_1d = chksum_int_1d( i(1:size(var)) )
      endif
      return
    end function chksum_r8_1d

    function chksum_r8_2d(var)
      integer(LONG_KIND) :: chksum_r8_2d
      real(DOUBLE_KIND), intent(in) :: var(:,:)
      real(DOUBLE_KIND) :: g(size(var,1),size(var,2))
      integer(LONG_KIND) :: i(size(var)*2)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          chksum_r8_2d = chksum_int_1d( i(1:size(var)*2) )
      else
          chksum_r8_2d = chksum_int_1d( i(1:size(var)) )
      endif
      return
    end function chksum_r8_2d

    function chksum_r8_3d(var)
      integer(LONG_KIND) :: chksum_r8_3d
      real(DOUBLE_KIND), intent(in) :: var(:,:,:)
      real(DOUBLE_KIND) :: g(size(var,1),size(var,2),size(var,3))
      integer(LONG_KIND) :: i(size(var)*2)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          chksum_r8_3d = chksum_int_1d( i(1:size(var)*2) )
      else
          chksum_r8_3d = chksum_int_1d( i(1:size(var)) )
      endif
      return
    end function chksum_r8_3d

    function chksum_r8_4d(var)
      integer(LONG_KIND) :: chksum_r8_4d
      real(DOUBLE_KIND), intent(in) :: var(:,:,:,:)
      real(DOUBLE_KIND) :: g(size(var,1),size(var,2),size(var,3),size(var,4))
      integer(LONG_KIND) :: i(size(var)*2)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          chksum_r8_4d = chksum_int_1d( i(1:size(var)*2) )
      else
          chksum_r8_4d = chksum_int_1d( i(1:size(var)) )
      endif
      return
    end function chksum_r8_4d

    function chksum_c8_0d(var)
      integer(LONG_KIND) :: chksum_c8_0d
      complex(DOUBLE_KIND), intent(in) :: var
      integer(LONG_KIND) :: i(4)
      pointer( ptr, i )
      ptr = LOC(var)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          chksum_c8_0d = chksum_int_1d( i(1:4) )
      else
          chksum_c8_0d = chksum_int_1d( i(1:2) )
      endif
      return
    end function chksum_c8_0d

    function chksum_c8_1d(var)
      integer(LONG_KIND) :: chksum_c8_1d
      complex(DOUBLE_KIND), intent(in) :: var(:)
      complex(DOUBLE_KIND) :: g(size(var))
      integer(LONG_KIND) :: i(size(var)*4)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          chksum_c8_1d = chksum_int_1d( i(1:size(var)*4) )
      else
          chksum_c8_1d = chksum_int_1d( i(1:size(var)*2) )
      endif
      return
    end function chksum_c8_1d

    function chksum_c8_2d(var)
      integer(LONG_KIND) :: chksum_c8_2d
      complex(DOUBLE_KIND), intent(in) :: var(:,:)
      complex(DOUBLE_KIND) :: g(size(var,1),size(var,2))
      integer(LONG_KIND) :: i(size(var)*4)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          chksum_c8_2d = chksum_int_1d( i(1:size(var)*4) )
      else
          chksum_c8_2d = chksum_int_1d( i(1:size(var)*2) )
      endif
      return
    end function chksum_c8_2d

    function chksum_c8_3d(var)
      integer(LONG_KIND) :: chksum_c8_3d
      complex(DOUBLE_KIND), intent(in) :: var(:,:,:)
      complex(DOUBLE_KIND) :: g(size(var,1),size(var,2),size(var,3))
      integer(LONG_KIND) :: i(size(var)*4)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          chksum_c8_3d = chksum_int_1d( i(1:size(var)*4) )
      else
          chksum_c8_3d = chksum_int_1d( i(1:size(var)*2) )
      endif
      return
    end function chksum_c8_3d

    function chksum_c8_4d(var)
      integer(LONG_KIND) :: chksum_c8_4d
      complex(DOUBLE_KIND), intent(in) :: var(:,:,:,:)
      complex(DOUBLE_KIND) :: g(size(var,1),size(var,2),size(var,3),size(var,4))
      integer(LONG_KIND) :: i(size(var)*4)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          chksum_c8_4d = chksum_int_1d( i(1:size(var)*4) )
      else
          chksum_c8_4d = chksum_int_1d( i(1:size(var)*2) )
      endif
      return
    end function chksum_c8_4d
#endif

  end module chksum_mod

#ifdef test_chksum  
  program test
    use chksum_mod

    real :: a, b(2)

!this series of checksums demonstrates why FP checksums do not verify all
!the bits of an answer: adding a number<epsilon to 1 results in the same
!FP answer but different integer answer.
    b(1) = 1.

    b(2) = 0.
    a = b(1) + b(2)
    print '(x,a,2z18)', 'chksum(1)          real, int=', chksum(a), chksum(b)

    b(2) = epsilon(1.)
    a = b(1) + b(2)
    print '(x,a,2z18)', 'chksum(1+eps)      real, int=', chksum(a), chksum(b)

    b(2) = 0.25*epsilon(1.)
    a = b(1) + b(2)
    print '(x,a,2z18)', 'chksum(1+0.25*eps) real, int=', chksum(a), chksum(b)
  end program test
#endif
