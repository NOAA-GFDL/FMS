
module fft_mod

!-----------------------------------------------------------------------
!these are used to determine hardware/OS/compiler

#ifdef __sgi
#  ifdef _COMPILER_VERSION
!the MIPSPro compiler defines _COMPILER_VERSION
#    define sgi_mipspro
#  else
#    define sgi_generic
#  endif
#endif

!fft uses the SCILIB on SGICRAY, and the NAG library otherwise
#if defined(_CRAY) || defined(sgi_mipspro)
#  define SGICRAY
#endif

use platform_mod, only: R8_KIND, R4_KIND
use      fms_mod, only: write_version_number,  &
                        error_mesg, FATAL
#ifndef SGICRAY
#ifndef NAGFFT
use    fft99_mod, only: fft991, set99
#endif
#endif

implicit none
private

!----------------- interfaces --------------------

public :: fft_init, fft_end, fft_grid_to_fourier, fft_fourier_to_grid

interface fft_grid_to_fourier
  module procedure fft_grid_to_fourier_float_2d, fft_grid_to_fourier_double_2d, &
                   fft_grid_to_fourier_float_3d, fft_grid_to_fourier_double_3d
end interface

interface fft_fourier_to_grid
  module procedure fft_fourier_to_grid_float_2d, fft_fourier_to_grid_double_2d, &
                   fft_fourier_to_grid_float_3d, fft_fourier_to_grid_double_3d
end interface

!---------------------- private data -----------------------------------

! tables for trigonometric constants and factors
! (not all will be used)
real(R8_KIND), allocatable, dimension(:) :: table8
real(R4_KIND), allocatable, dimension(:) :: table4
real         , allocatable, dimension(:) :: table99
integer      , allocatable, dimension(:) :: ifax

logical :: do_log =.true.
integer :: leng, leng1, leng2, lenc    ! related to transform size

logical :: module_is_initialized=.false.

!  cvs version and tag name
character(len=128) :: version = '$Id: fft.F90,v 1.3 2002/07/16 22:55:18 fms Exp $'
character(len=128) :: tagname = '$Name: havana $'

!-----------------------------------------------------------------------
!
!                    WRAPPER FOR FFT
!
!   Provides fast fourier transtorm (FFT) between real grid
!   space and complex fourier space.
!
!   The complex fourier components are passed in the following format.
!
!        fourier (1)     = cmplx ( a(0), b(0) )
!        fourier (2)     = cmplx ( a(1), b(1) )
!            :              :
!            :              :
!        fourier (n/2+1) = cmplx ( a(n/2), b(n/2) )
!
!   where n = length of each real transform
!
!   fft uses the SCILIB on SGICRAY, otherwise the NAG library or
!   a standalone version of Temperton's fft is used
!     SCFFTM and CSFFTM are used on Crays
!     DZFFTM and ZDFFTM are used on SGIs
!   The following NAG routines are used: c06fpf, c06gqf, c06fqf.
!   These routine names may be slightly different on different 
!   platforms.
!
!-----------------------------------------------------------------------

contains

!#######################################################################

 function fft_grid_to_fourier_float_2d (grid) result (fourier)

!-----------------------------------------------------------------------

   real   (R4_KIND), intent(in),  dimension(:,:)  :: grid
   complex(R4_KIND), dimension(lenc,size(grid,2)) :: fourier

!-----------------------------------------------------------------------
!
!  input
!  -----
!   grid = Multiple transforms in grid point space, the first dimension
!          must be n+1 (where n is the size of each real transform).
!
!  returns
!  -------
!    Multiple transforms in complex fourier space, the first dimension
!    must equal n/2+1 (where n is the size of each real transform).
!    The remaining dimensions must be the same size as the input
!    argument "grid".
!
!-----------------------------------------------------------------------
#ifdef SGICRAY
#  ifdef _CRAY
!  local storage for cray fft
   real(R4_KIND), dimension((2*leng+4)*size(grid,2)) :: work
#  else
!  local storage for sgi fft
   real(R4_KIND), dimension(leng2) :: work
#  endif
#else
#  ifdef NAGFFT
!  local storage for nag fft
   real(R4_KIND), dimension(size(grid,2),leng) :: data, work
#  else
!  local storage for temperton fft
   real, dimension(leng2,size(grid,2)) :: data
   real, dimension(leng1,size(grid,2)) :: work
#  endif   
#endif   

   real(R4_KIND) :: scale
   integer :: j, k, num, len_grid, ifail

!-----------------------------------------------------------------------

      if (.not.module_is_initialized) &
        call error_handler ('fft_grid_to_fourier',  &
                            'fft_init must be called.')

!-----------------------------------------------------------------------

      len_grid = size(grid,1)
#ifdef SGICRAY
      if (len_grid /= leng1) call error_handler ('fft_grid_to_fourier',  &
                        'size of first dimension of input data is wrong')
#else
      if (len_grid < leng) call error_handler ('fft_grid_to_fourier',  &
                                   'length of input data too small.')
#endif
!-----------------------------------------------------------------------
!----------------transform to fourier coefficients (+1)-----------------

      num   = size(grid,2)    ! number of transforms

#ifdef SGICRAY
!  Cray/SGI fft
      scale = 1./real(leng)
#  ifdef _CRAY
      call scfftm (-1,leng,num,scale, grid,leng1, fourier,lenc,  &
                   table4, work, 0)
#  else
      call scfftm (-1,leng,num,scale, grid,leng1, fourier,lenc,  &
                   table4, work, 0)
#  endif
#else
#  ifdef NAGFFT
!  NAG fft
!  will not allow float kind for NAG
      call error_handler ('fft_grid_to_fourier',  &
                          'float kind not supported for nag fft')
      do j=1,size(grid,2)
         data(j,1:leng) = grid(1:leng,j)
      enddo
!!!!! call c06fpe ( num, leng, data, 's', table4, work, ifail )
      scale = 1./sqrt(float(leng))
      data = data * scale
      fourier(1,:) = cmplx( data(:,1), 0. )
      do k=2,lenc-1
         fourier(k,:) = cmplx( data(:,k), data(:,leng-k+2) )
      enddo
      fourier(lenc,:) = cmplx( data(:,lenc), 0. )
#  else
!  Temperton fft
      do j=1,num
        data(1:leng,j) = grid(1:leng,j)
      enddo
      call fft991 (data,work,table99,ifax,1,leng2,leng,num,-1)
      do j=1,size(grid,2)
      do k=1,lenc
        fourier(k,j) = cmplx( data(2*k-1,j), data(2*k,j) )
      enddo
      enddo
#  endif
#endif
!-----------------------------------------------------------------------

 end function fft_grid_to_fourier_float_2d

!#######################################################################

 function fft_fourier_to_grid_float_2d (fourier) result (grid)

!-----------------------------------------------------------------------

   complex(R4_KIND),  intent(in),  dimension(:,:)     :: fourier
   real   (R4_KIND), dimension(leng1,size(fourier,2)) :: grid

!-----------------------------------------------------------------------
!
!  input
!  -----
!  fourier = Multiple transforms in complex fourier space, the first 
!            dimension must equal n/2+1 (where n is the size of each
!            real transform).
!
!  returns
!  -------
!    Multiple transforms in grid point space, the first dimension
!    must be n+1 (where n is the size of each real transform).
!    The remaining dimensions must be the same size as the input
!    argument "fourier".
!
!-----------------------------------------------------------------------
#ifdef SGICRAY
#  ifdef _CRAY
!  local storage for cray fft
   real(R4_KIND), dimension((2*leng+4)*size(fourier,2)) :: work
#  else
!  local storage for sgi fft
   real(R4_KIND), dimension(leng2) :: work
#  endif
#else
#  ifdef NAGFFT
!  local storage for nag fft
   real(R4_KIND), dimension(size(fourier,2),leng) :: data, work
#  else
!  local storage for temperton fft
   real, dimension(leng2,size(fourier,2)) :: data
   real, dimension(leng1,size(fourier,2)) :: work
#  endif   
#endif   

   real(R4_KIND) :: scale
   integer :: j, k, num, len_fourier, ifail

!-----------------------------------------------------------------------

   if (.not.module_is_initialized) &
      call error_handler ('fft_grid_to_fourier',  &
                          'fft_init must be called.')

!-----------------------------------------------------------------------

      len_fourier = size(fourier,1)
      num         = size(fourier,2)    ! number of transforms

#ifdef SGICRAY
      if (len_fourier /= lenc) call error_handler ('fft_fourier_to_grid', &
               'size of first dimension of input data is wrong')
#else
      if (len_fourier < lenc) call error_handler ('fft_fourier_to_grid',  &
                               'length of input data too small.')
#endif
!-----------------------------------------------------------------------
!----------------inverse transform to real space (-1)-------------------

#ifdef SGICRAY
!  Cray/SGI fft
      scale = 1.0
#  ifdef _CRAY
      call csfftm (+1,leng,num,scale, fourier,len_fourier,  &
                     grid,leng1, table4, work, 0)
#  else
      call csfftm (+1,leng,num,scale, fourier,len_fourier,  &
                     grid,leng1, table4, work, 0)
#  endif
#else
#  ifdef NAGFFT
!  NAG fft
!  will not allow float kind for nag
      call error_handler ('fft_fourier_to_grid',  &
                          'float kind not supported for nag fft')

  ! save input complex array in real format (herm.)
      do k=1,lenc
         data(:,k) = real(fourier(k,:))
      enddo
      do k=2,lenc-1
         data(:,leng-k+2) = aimag(fourier(k,:))
      enddo

!!!!! call c06gqe ( num, leng, data, ifail )
!!!!! call c06fqe ( num, leng, data, 's', table4, work, ifail )

  ! scale and transpose data
      scale = sqrt(real(leng))
      do j=1,num
         grid(1:leng,j) = data(j,1:leng)*scale
      enddo
#  else
!  Temperton fft
      do j=1,num
      do k=1,lenc
         data(2*k-1,j) = real (fourier(k,j))
         data(2*k  ,j) = aimag(fourier(k,j))
      enddo
      enddo
      call fft991 (data,work,table99,ifax,1,leng2,leng,num,+1)
      do j=1,num
         grid(1:leng,j) = data(1:leng,j)
      enddo
#  endif
#endif

!-----------------------------------------------------------------------

 end function fft_fourier_to_grid_float_2d

!#######################################################################

 function fft_grid_to_fourier_double_2d (grid) result (fourier)

!-----------------------------------------------------------------------

   real   (R8_KIND), intent(in),  dimension(:,:)  :: grid
   complex(R8_KIND), dimension(lenc,size(grid,2)) :: fourier

!-----------------------------------------------------------------------
!
!  input
!  -----
!   grid = Multiple transforms in grid point space, the first dimension
!          must be n+1 (where n is the size of each real transform).
!
!  returns
!  -------
!    Multiple transforms in complex fourier space, the first dimension
!    must equal n/2+1 (where n is the size of each real transform).
!    The remaining dimensions must be the same size as the input
!    argument "grid".
!
!-----------------------------------------------------------------------
#ifdef SGICRAY
#  ifdef _CRAY
!  local storage for cray fft
   real(R8_KIND), dimension((2*leng+4)*size(grid,2)) :: work
#  else
!  local storage for sgi fft
   real(R8_KIND), dimension(leng2) :: work
#  endif
#else
#  ifdef NAGFFT
!  local storage for nag fft
   real(R8_KIND), dimension(size(grid,2),leng) :: data, work
#  else
!  local storage for temperton fft
   real, dimension(leng2,size(grid,2)) :: data
   real, dimension(leng1,size(grid,2)) :: work
#  endif   
#endif   

   real(R8_KIND) :: scale
   integer :: j, k, num, len_grid, ifail

!-----------------------------------------------------------------------

      if (.not.module_is_initialized) &
        call error_handler ('fft_grid_to_fourier',  &
                            'fft_init must be called.')

!-----------------------------------------------------------------------

      len_grid = size(grid,1)
#ifdef SGICRAY
      if (len_grid /= leng1) call error_handler ('fft_grid_to_fourier',  &
                        'size of first dimension of input data is wrong')
#else
      if (len_grid < leng) call error_handler ('fft_grid_to_fourier',  &
                                   'length of input data too small.')
#endif
!-----------------------------------------------------------------------
!----------------transform to fourier coefficients (+1)-----------------

      num   = size(grid,2)    ! number of transforms
#ifdef SGICRAY
!  Cray/SGI fft
      scale = 1./float(leng)
#  ifdef _CRAY
      call scfftm (-1,leng,num,scale, grid,leng1, fourier,lenc,  &
                   table8, work, 0)
#  else
      call dzfftm (-1,leng,num,scale, grid,leng1, fourier,lenc,  &
                   table8, work, 0)
#  endif
#else
#  ifdef NAGFFT
!  NAG fft
      do j=1,size(grid,2)
         data(j,1:leng) = grid(1:leng,j)
      enddo
      call c06fpf ( num, leng, data, 's', table8, work, ifail )
      scale = 1./sqrt(float(leng))
      data = data * scale
      fourier(1,:) = cmplx( data(:,1), 0. )
      do k=2,lenc-1
         fourier(k,:) = cmplx( data(:,k), data(:,leng-k+2) )
      enddo
      fourier(lenc,:) = cmplx( data(:,lenc), 0. )
#  else
!  Temperton fft
      do j=1,num
        data(1:leng,j) = grid(1:leng,j)
      enddo
      call fft991 (data,work,table99,ifax,1,leng2,leng,num,-1)
      do j=1,size(grid,2)
      do k=1,lenc
        fourier(k,j) = cmplx( data(2*k-1,j), data(2*k,j) )
      enddo
      enddo
#  endif
#endif
!-----------------------------------------------------------------------

 end function fft_grid_to_fourier_double_2d

!#######################################################################

 function fft_fourier_to_grid_double_2d (fourier) result (grid)

!-----------------------------------------------------------------------

   complex(R8_KIND),  intent(in),  dimension(:,:)     :: fourier
   real   (R8_KIND), dimension(leng1,size(fourier,2)) :: grid

!-----------------------------------------------------------------------
!
!  input
!  -----
!  fourier = Multiple transforms in complex fourier space, the first 
!            dimension must equal n/2+1 (where n is the size of each
!            real transform).
!
!  returns
!  -------
!    Multiple transforms in grid point space, the first dimension
!    must be n+1 (where n is the size of each real transform).
!    The remaining dimensions must be the same size as the input
!    argument "fourier".
!
!-----------------------------------------------------------------------
#ifdef SGICRAY
#  ifdef _CRAY
!  local storage for cray fft
   real(R8_KIND), dimension((2*leng+4)*size(fourier,2)) :: work
#  else
!  local storage for sgi fft
   real(R8_KIND), dimension(leng2) :: work
#  endif
#else
#  ifdef NAGFFT
!  local storage for nag fft
   real(R8_KIND), dimension(size(fourier,2),leng) :: data, work
#  else
!  local storage for temperton fft
   real, dimension(leng2,size(fourier,2)) :: data
   real, dimension(leng1,size(fourier,2)) :: work
#  endif   
#endif   

   real(R8_KIND) :: scale
   integer :: j, k, num, len_fourier, ifail

!-----------------------------------------------------------------------

      if (.not.module_is_initialized) &
        call error_handler ('fft_grid_to_fourier',  &
                            'fft_init must be called.')

!-----------------------------------------------------------------------

      len_fourier = size(fourier,1)
      num         = size(fourier,2)    ! number of transforms

#ifdef SGICRAY
      if (len_fourier /= lenc) call error_handler ('fft_fourier_to_grid', &
               'size of first dimension of input data is wrong')
#else
      if (len_fourier < lenc) call error_handler ('fft_fourier_to_grid',  &
                               'length of input data too small.')
#endif
!-----------------------------------------------------------------------
!----------------inverse transform to real space (-1)-------------------

#ifdef SGICRAY
!  Cray/SGI fft
      scale = 1.0
#  ifdef _CRAY
      call csfftm (+1,leng,num,scale, fourier,len_fourier,  &
                     grid,leng1, table8, work, 0)
#  else
      call zdfftm (+1,leng,num,scale, fourier,len_fourier,  &
                     grid,leng1, table8, work, 0)
#  endif
#else
#  ifdef NAGFFT
!  NAG fft

    ! save input complex array in real format (herm.)
      do k=1,lenc
         data(:,k) = real(fourier(k,:))
      enddo
      do k=2,lenc-1
         data(:,leng-k+2) = aimag(fourier(k,:))
      enddo

      call c06gqf ( num, leng, data, ifail )
      call c06fqf ( num, leng, data, 's', table8, work, ifail )

    ! scale and transpose data
      scale = sqrt(real(leng))
      do j=1,num
         grid(1:leng,j) = data(j,1:leng)*scale
      enddo
#  else
!  Temperton fft
      do j=1,num
      do k=1,lenc
         data(2*k-1,j) = real (fourier(k,j))
         data(2*k  ,j) = aimag(fourier(k,j))
      enddo
      enddo
      call fft991 (data,work,table99,ifax,1,leng2,leng,num,+1)
      do j=1,num
         grid(1:leng,j) = data(1:leng,j)
      enddo
#  endif
#endif

!-----------------------------------------------------------------------

 end function fft_fourier_to_grid_double_2d

!#######################################################################
!                   interface overloads
!#######################################################################

 function fft_grid_to_fourier_float_3d (grid) result (fourier)

!-----------------------------------------------------------------------
   real   (R4_KIND),    intent(in),  dimension(:,:,:) :: grid
   complex(R4_KIND), dimension(lenc,size(grid,2),size(grid,3)) :: fourier
   integer :: n
!-----------------------------------------------------------------------

    do n = 1, size(grid,3)
      fourier(:,:,n) = fft_grid_to_fourier_float_2d (grid(:,:,n))
    enddo

!-----------------------------------------------------------------------

 end function fft_grid_to_fourier_float_3d

!#######################################################################

 function fft_fourier_to_grid_float_3d (fourier) result (grid)

!-----------------------------------------------------------------------
   complex(R4_KIND),  intent(in),  dimension(:,:,:) :: fourier
   real   (R4_KIND), dimension(leng1,size(fourier,2),size(fourier,3)) :: grid
   integer :: n
!-----------------------------------------------------------------------

    do n = 1, size(fourier,3)
      grid(:,:,n) = fft_fourier_to_grid_float_2d (fourier(:,:,n))
    enddo

!-----------------------------------------------------------------------

 end function fft_fourier_to_grid_float_3d

!#######################################################################

 function fft_grid_to_fourier_double_3d (grid) result (fourier)

!-----------------------------------------------------------------------
   real   (R8_KIND),    intent(in),  dimension(:,:,:) :: grid
   complex(R8_KIND), dimension(lenc,size(grid,2),size(grid,3)) :: fourier
   integer :: n
!-----------------------------------------------------------------------

    do n = 1, size(grid,3)
      fourier(:,:,n) = fft_grid_to_fourier_double_2d (grid(:,:,n))
    enddo

!-----------------------------------------------------------------------

 end function fft_grid_to_fourier_double_3d

!#######################################################################

 function fft_fourier_to_grid_double_3d (fourier) result (grid)

!-----------------------------------------------------------------------
   complex(R8_KIND),  intent(in),  dimension(:,:,:) :: fourier
   real   (R8_KIND), dimension(leng1,size(fourier,2),size(fourier,3)) :: grid
   integer :: n
!-----------------------------------------------------------------------

    do n = 1, size(fourier,3)
      grid(:,:,n) = fft_fourier_to_grid_double_2d (fourier(:,:,n))
    enddo

!-----------------------------------------------------------------------

 end function fft_fourier_to_grid_double_3d

!#######################################################################
!#######################################################################

 subroutine fft_init (n)

!-----------------------------------------------------------------------
   integer, intent(in) :: n
!-----------------------------------------------------------------------
!
!   n = size (length) of each transform
!
!-----------------------------------------------------------------------
#ifdef SGICRAY
   real   (R4_KIND) ::  dummy4(1)
   complex(R4_KIND) :: cdummy4(1)
   real   (R8_KIND) ::  dummy8(1)
   complex(R8_KIND) :: cdummy8(1)
   integer :: isys(0:1)
#else
#  ifdef NAGFFT
   real(R8_KIND) :: data8(n), work8(n)
   real(R4_KIND) :: data4(n), work4(n)
   integer       :: ifail4, ifail8
#  endif
#endif
!-----------------------------------------------------------------------
!   --- fourier transform initialization ----

      if (module_is_initialized) &
      call error_handler ('fft_init', 'attempted to reinitialize fft')

!  write version and tag name to log file
   if (do_log) then
      call write_version_number (version, tagname)
      do_log = .false.
   endif

!  variables that save length of transform
      leng = n; leng1 = n+1; leng2 = n+2; lenc = n/2+1

#ifdef SGICRAY
#  ifdef _CRAY
!  initialization for cray
!  float kind may not apply for cray
      allocate (table4(100+2*leng), table8(100+2*leng))   ! size may be too large?
      call scfftm (0,leng,1,0.0, dummy4, 1, cdummy4, 1, table4, dummy4, 0)
      call scfftm (0,leng,1,0.0, dummy8, 1, cdummy8, 1, table8, dummy8, 0)
#  else
!  initialization for sgi
      allocate (table4(leng+256), table8(leng+256))
      isys(0) = 1
      call scfftm (0,leng,1,0.0, dummy4, 1, cdummy4, 1, table4, dummy8, isys)
      call dzfftm (0,leng,1,0.0, dummy8, 1, cdummy8, 1, table8, dummy8, isys)
#  endif
#else
#  ifdef NAGFFT
!  initialization for nag fft
      ifail8 = 0
      allocate (table8(100+2*leng))   ! size may be too large?
      call c06fpf ( 1, leng, data8, 'i', table8, work8, ifail8 )

!  will not allow float kind for nag
      ifail4 = 0
!!!!! allocate (table4(100+2*leng))
!!!!! call c06fpe ( 1, leng, data4, 'i', table4, work4, ifail4 )

      if (ifail4 /= 0 .or. ifail8 /= 0) then
          call error_handler ('fft_init', 'nag fft initialization error')
      endif
#  else
!  initialization for Temperton fft
      allocate (table99(3*leng/2+1))
      allocate (ifax(10))
      call set99 ( table99, ifax, leng )
#  endif
#endif

      module_is_initialized = .true.

!-----------------------------------------------------------------------

 end subroutine fft_init

!#######################################################################

 subroutine fft_end

!-----------------------------------------------------------------------
!
!   unsets transform size and deallocates memory
!
!-----------------------------------------------------------------------
!   --- fourier transform un-initialization ----

      if (.not.module_is_initialized) &
        call error_handler ('fft_end', &
          'attempt to un-initialize fft that has not been initialized')

      leng = 0; leng1 = 0; leng2 = 0; lenc = 0

      if (allocated(table4))  deallocate (table4)
      if (allocated(table8))  deallocate (table8)
      if (allocated(table99)) deallocate (table99)

      module_is_initialized = .false.

!-----------------------------------------------------------------------

 end subroutine fft_end

!#######################################################################
! wrapper for handling errors

 subroutine error_handler ( routine, message )
 character(len=*), intent(in) :: routine, message

   call error_mesg ( routine, message, FATAL )

!  print *, 'ERROR: ',trim(routine)
!  print *, 'ERROR: ',trim(message)
!  stop 111

 end subroutine error_handler

!#######################################################################

end module fft_mod

#ifdef test_fft
program test
use fft_mod
integer, parameter :: lot = 2
real   , allocatable :: ain(:,:), aout(:,:)
complex, allocatable :: four(:,:)
integer :: i, j, m, n
integer :: ntrans(2) = (/ 60, 90 /)

! test multiple transform lengths
  do m = 1,2

  ! set up input data
    n = ntrans(m)
    allocate (ain(n+1,lot),aout(n+1,lot),four(n/2+1,lot))
    call random_number (ain(1:n,:))
    aout(1:n,:) = ain(1:n,:)

    call fft_init (n)
  ! transform grid to fourier and back
    four = fft_grid_to_fourier (aout)
    aout = fft_fourier_to_grid (four)

  ! print original and transformed
    do j=1,lot
    do i=1,n
      write (*,'(2i4,3(2x,f15.9))') j, i, ain(i,j), aout(i,j), aout(i,j)-ain(i,j)
    enddo
    enddo

    call fft_end
    deallocate (ain,aout,four)
  enddo

end program test
#endif

