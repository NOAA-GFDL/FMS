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

module fft_mod

!-----------------------------------------------------------------------

use utilities_mod, only: error_mesg, FATAL

implicit none
private

!----------------- interfaces --------------------

public   fft_init, fft_grid_to_fourier, fft_fourier_to_grid

interface fft_grid_to_fourier
  module procedure fft_grid_to_fourier_2d, fft_grid_to_fourier_3d
end interface
  private          fft_grid_to_fourier_2d, fft_grid_to_fourier_3d

interface fft_fourier_to_grid
  module procedure fft_fourier_to_grid_2d, fft_fourier_to_grid_3d
end interface
  private          fft_fourier_to_grid_2d, fft_fourier_to_grid_3d

!---------------------- private data -----------------------------------

integer,parameter :: r8_kind = selected_real_kind(15,307)

real(r8_kind), allocatable, dimension(:)   :: table

logical :: do_init=.true.
integer :: len,lenp1,lenc

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
!fft uses the SCILIB on SGICRAY, and the NAG library otherwise
!   SCFFTM and CSFFTM are used on Crays
!   DZFFTM and ZDFFTM are used on SGIs
!   The following NAG routines are used: c06fpf, c06gqf, c06fqf.
!   These routine names may be slightly different on different 
!   platforms.
!
!-----------------------------------------------------------------------

contains

!#######################################################################

 function fft_grid_to_fourier_2d (grid) result (fourier)

!-----------------------------------------------------------------------

   real,    intent(in),  dimension(:,:)  :: grid
   complex, dimension(lenc,size(grid,2)) :: fourier

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
!      ----------- local storage -------------
#ifdef SGICRAY
   real(r8_kind), dimension((2*len+4)*size(grid,2)) :: work
#else
   real(r8_kind), dimension(size(grid,2),len) :: data,work
#endif   
   real     scale
   integer  j,k, num, len_grid, ierr

!-----------------------------------------------------------------------

      if (do_init) call error_mesg ('fft_grid_to_fourier',  &
                                   'fft_init must be called.', FATAL)

!-----------------------------------------------------------------------

      len_grid = size(grid,1)
#ifdef SGICRAY
      if (len_grid /= lenp1) call error_mesg ('fft_grid_to_fourier',  &
              'size of first dimension of input data is wrong', FATAL)
#else
      if (len_grid < len) call error_mesg ('fft_grid_to_fourier',  &
                               'length of input data too small.', FATAL)
#endif
!-----------------------------------------------------------------------
!----------------transform to fourier coefficients (+1)-----------------

      num   = size(grid,2)
#ifdef SGICRAY
      scale = 1./float(len)
#  ifdef _CRAY
      call scfftm (-1,len,num,scale, grid,lenp1, fourier,lenc,  &
                   table, work, 0)
#  else
      call dzfftm (-1,len,num,scale, grid,lenp1, fourier,lenc,  &
                   table, work, 0)
#  endif
#else
      do j=1,size(grid,2)
         data(j,1:len) = grid(1:len,j)
      enddo
      call c06fpf ( num, len, data, 's', table, work, ierr )
      scale = 1./sqrt(float(len))
      data = data * scale
      
      
      fourier(1,:) = cmplx( data(:,1), 0. )
      do k=2,lenc-1
         fourier(k,:) = cmplx( data(:,k), data(:,len-k+2) )
      enddo
      fourier(lenc,:) = cmplx( data(:,lenc), 0. )
#endif
!-----------------------------------------------------------------------

 end function fft_grid_to_fourier_2d

!#######################################################################

 function fft_fourier_to_grid_2d (fourier) result (grid)

!-----------------------------------------------------------------------

   complex,  intent(in),  dimension(:,:)  :: fourier
   real, dimension(lenp1,size(fourier,2)) :: grid

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
!   ----------- local storage -------------
#ifdef SGICRAY
   real(r8_kind), dimension((2*len+4)*size(fourier,2)) :: work
#else
   real(r8_kind), dimension(size(fourier,2),len) :: data,work
#endif
   real     scale
   integer  j,k, num, len_fourier, ierr

!-----------------------------------------------------------------------

      if (do_init) call error_mesg ('fft_grid_to_fourier',  &
                                    'fft_init must be called.', FATAL)

!-----------------------------------------------------------------------

      len_fourier  = size(fourier,1)
#ifdef SGICRAY
      if (len_fourier /= lenc) call error_mesg ('fft_fourier_to_grid', &
               'size of first dimension of input data is wrong', FATAL)
#else
      if (len_fourier < lenc) call error_mesg ('fft_fourier_to_grid',  &
                               'length of input data too small.', FATAL)
#endif
!-----------------------------------------------------------------------
!----------------inverse transform to real space (-1)-------------------

      num   = size(fourier,2)
#ifdef SGICRAY
      scale = 1.0
#  ifdef _CRAY
      call csfftm (+1,len,num,scale, fourier,len_fourier,  &
                     grid,lenp1, table, work, 0)
#  else
      call zdfftm (+1,len,num,scale, fourier,len_fourier,  &
                     grid,lenp1, table, work, 0)
#  endif
#else
!-----------------------------------------------------------------------
!  save input complex array in real format (herm.)

      do k=1,lenc
         data(:,k) = real(fourier(k,:))
      enddo
      do k=2,lenc-1
         data(:,len-k+2) = aimag(fourier(k,:))
      enddo

      call c06gqf ( num, len, data, ierr )
      call c06fqf ( num, len, data, 's', table, work, ierr )

!---------- scale and transpose data --------------

      scale = sqrt(float(len))

      do j=1,size(fourier,2)
         grid(1:len,j) = data(j,1:len)*scale
      enddo
!-----------------------------------------------------------------------
#endif
 end function fft_fourier_to_grid_2d

!#######################################################################

 function fft_grid_to_fourier_3d (grid) result (fourier)

!-----------------------------------------------------------------------
   real,    intent(in),  dimension(:,:,:) :: grid
   complex, dimension(lenc,size(grid,2),size(grid,3)) :: fourier
   integer :: n
!-----------------------------------------------------------------------

    do n = 1, size(grid,3)
      fourier(:,:,n) = fft_grid_to_fourier_2d (grid(:,:,n))
    enddo

!-----------------------------------------------------------------------

 end function fft_grid_to_fourier_3d

!#######################################################################

 function fft_fourier_to_grid_3d (fourier) result (grid)

!-----------------------------------------------------------------------
   complex,  intent(in),  dimension(:,:,:) :: fourier
   real, dimension(lenp1,size(fourier,2),size(fourier,3)) :: grid
   integer :: n
!-----------------------------------------------------------------------

    do n = 1, size(fourier,3)
      grid(:,:,n) = fft_fourier_to_grid_2d (fourier(:,:,n))
    enddo

!-----------------------------------------------------------------------

 end function fft_fourier_to_grid_3d

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
   real       dummy(1)
   complex   cdummy(1)
#else
   real(r8_kind) :: data(n), work(n)
   integer       :: ierr
#endif
!-----------------------------------------------------------------------

!   --- fourier transform initialization ----

      len = n; lenp1 = len + 1; lenc = len / 2 + 1
      allocate (table(100+2*len))
#ifdef SGICRAY
#  ifdef _CRAY
      call scfftm (0,len,1,0.0, dummy, 1, cdummy, 1, table, dummy, 0)
#  else
      call dzfftm (0,len,1,0.0, dummy, 1, cdummy, 1, table, dummy, 0)
#  endif
#else
      ierr = 0
      call c06fpf ( 1, len, data, 'i', table, work, ierr )

      if (ierr /= 0) then
          call error_mesg ('fft_init', 'nag fft initialization error', &
                           FATAL)
      endif
#endif
      do_init = .false.

!-----------------------------------------------------------------------

 end subroutine fft_init

!#######################################################################

end module fft_mod

