!!#include <fms_platform.h>

module drifters_io_mod
  implicit none  
  private

  public :: drifters_io_type, drifters_io_new, drifters_io_del, drifters_io_set_time_units
  public :: drifters_io_set_position_names, drifters_io_set_position_units, drifters_io_set_field_names
  public :: drifters_io_set_field_units, drifters_io_write

  ! Globals
  integer, parameter, private   :: MAX_STR_LEN = 128
  ! Include variable "version" to be written to log file.
#include<file_version.h>

  real :: drfts_eps_t = 10.*epsilon(1.)
  

  type drifters_io_type
     real                 :: time
     integer              :: it ! time index
     integer              :: it_id  ! infinite axis index
     integer              :: ncid
     integer              :: nc_positions, nc_fields, nc_ids, nc_time, nc_index_time    
     logical              :: enddef
  end type drifters_io_type

contains

!###############################################################################
  subroutine drifters_io_new(self, filename, nd, nf, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(in)  :: filename
    integer, intent(in)           :: nd  ! number of dims
    integer, intent(in)           :: nf  ! number of fields
    character(len=*), intent(out) :: ermesg

    integer ier, nc_it_id, nc_nd, nc_nf
    integer :: size1(1), size2(2)
    include 'netcdf.inc'

    ermesg=''
    self%enddef = .FALSE.

    ier = nf_create(filename, NF_CLOBBER, self%ncid)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_create ('//filename//') '//nf_strerror(ier)

    ! global attributes
    ier = nf_put_att_text(self%ncid, NF_GLOBAL, 'version', len_trim(version), trim(version))
    

    ! dimensions
    ier = nf_def_dim(self%ncid, 'it_id', NF_UNLIMITED, nc_it_id)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_dim (it_id) '//nf_strerror(ier)

    ier = nf_def_dim(self%ncid, 'nf', nf, nc_nf)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_dim (nf) '//nf_strerror(ier)

    ier = nf_def_dim(self%ncid, 'nd', nd, nc_nd)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_dim (nd) '//nf_strerror(ier)

    ! variables
    size1 = (/nc_it_id/)
    ier = nf_def_var(self%ncid, 'index_time', NF_INT, 1, size1, self%nc_index_time)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_var (index_time)'//nf_strerror(ier)

    ier = nf_def_var(self%ncid, 'time', NF_DOUBLE, 1, size1, self%nc_time)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_var (time)'//nf_strerror(ier)

    ier = nf_def_var(self%ncid, 'ids', NF_INT, 1, size1, self%nc_ids)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_var (ids)'//nf_strerror(ier)

    size2 = (/nc_nd, nc_it_id/)
    ier = nf_def_var(self%ncid, 'positions', NF_DOUBLE, 2, size2, self%nc_positions)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_var (positions)'//nf_strerror(ier)

    size2 = (/nc_nf, nc_it_id/)
    ier = nf_def_var(self%ncid, 'fields', NF_DOUBLE, 2, size2, self%nc_fields)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_var (fields)'//nf_strerror(ier)
    
    self%time  = -huge(1.)
    self%it    = -1
    self%it_id = 1

  end subroutine drifters_io_new

!###############################################################################
  subroutine drifters_io_del(self, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(out) :: ermesg

    integer ier
    include 'netcdf.inc'

    ermesg = ''
    
    ier = nf_close(self%ncid)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_del::nf_close '//nf_strerror(ier)
    
  end subroutine drifters_io_del

!###############################################################################
  subroutine drifters_io_set_time_units(self, name, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: ermesg

    integer ier
    include 'netcdf.inc'

    ermesg = ''
    ier = nf_put_att_text(self%ncid, NF_GLOBAL, &
         & 'time_units', len_trim(name), trim(name))
    if(ier/=NF_NOERR) &
         & ermesg = 'drifters_io_set_time_units::failed to add time_units attribute ' &
         & //nf_strerror(ier)

  end subroutine drifters_io_set_time_units

!###############################################################################
  subroutine drifters_io_set_position_names(self, names, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(in)  :: names(:)
    character(len=*), intent(out) :: ermesg

    integer n, ier, i
    character(len=128) :: attname
    include 'netcdf.inc'

    n = size(names)
    ermesg = ''

    do i = 1, n
       write(attname, '(i6)' ) i
       attname = 'name_'//adjustl(attname)
       ier = nf_put_att_text(self%ncid, self%nc_positions, &
            & trim(attname), len_trim(names(i)), trim(names(i)))
       if(ier/=NF_NOERR) &
            & ermesg = 'drifters_io_set_position_names::failed to add name attribute to positions '//nf_strerror(ier)
    enddo

  end subroutine drifters_io_set_position_names

!###############################################################################
  subroutine drifters_io_set_position_units(self, names, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(in)  :: names(:)
    character(len=*), intent(out) :: ermesg

    integer n, ier, i
    character(len=128) :: attname
    include 'netcdf.inc'

    n = size(names)
    ermesg = ''

    do i = 1, n
       write(attname, '(i6)' ) i
       attname = 'unit_'//adjustl(attname)
       ier = nf_put_att_text(self%ncid, self%nc_positions, &
            & trim(attname), len_trim(names(i)), trim(names(i)))
       if(ier/=NF_NOERR) &
            & ermesg = 'drifters_io_set_position_names::failed to add unit attribute to positions '//nf_strerror(ier)
    enddo
        
  end subroutine drifters_io_set_position_units

!###############################################################################
  subroutine drifters_io_set_field_names(self, names, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(in)  :: names(:)
    character(len=*), intent(out) :: ermesg

    integer n, ier, i
    character(len=128) :: attname
    include 'netcdf.inc'

    n = size(names)
    ermesg = ''

    do i = 1, n
       write(attname, '(i6)' ) i
       attname = 'name_'//adjustl(attname)
       ier = nf_put_att_text(self%ncid, self%nc_fields, &
            & trim(attname), len_trim(names(i)), trim(names(i)))
       if(ier/=NF_NOERR) &
            & ermesg = 'drifters_io_set_field_names::failed to add name attribute to fields '//nf_strerror(ier)
    enddo

  end subroutine drifters_io_set_field_names

!###############################################################################
  subroutine drifters_io_set_field_units(self, names, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(in)  :: names(:)
    character(len=*), intent(out) :: ermesg

    integer n, ier, i
    character(len=128) :: attname
    include 'netcdf.inc'

    n = size(names)
    ermesg = ''

    do i = 1, n
       write(attname, '(i6)' ) i
       attname = 'unit_'//adjustl(attname)
       ier = nf_put_att_text(self%ncid, self%nc_fields, &
            & trim(attname), len_trim(names(i)), trim(names(i)))
       if(ier/=NF_NOERR) &
            & ermesg = 'drifters_io_set_field_units::failed to add unit attribute to fields '//nf_strerror(ier)
    enddo
    
  end subroutine drifters_io_set_field_units
!###############################################################################

  subroutine drifters_io_write(self, time, np, nd, nf, ids, positions, fields, ermesg)
    type(drifters_io_type)        :: self
    real, intent(in)              :: time
    integer, intent(in)           :: np    ! number of dirfters
    integer, intent(in)           :: nd    ! number of dimensions
    integer, intent(in)           :: nf    ! number of fields
    integer, intent(in)           :: ids(np)          ! of size np
    real, intent(in)              :: positions(nd,np) ! nd times np
    real, intent(in)              :: fields(nf,np)    ! nf times np
    character(len=*), intent(out) :: ermesg

    integer ier, i
    integer :: start1(1), len1(1), start2(2), len2(2)
    integer :: it_indices(np)
    real    :: time_array(np)
    include 'netcdf.inc'

    ermesg = ''
    
    if(.not. self%enddef) then
       ier = nf_enddef(self%ncid)
       if(ier/=NF_NOERR) then 
            ermesg = 'drifters_io_write::nf_enddef failure. No data will be written. '//nf_strerror(ier)
            return
       endif
       self%enddef = .TRUE.
    endif

    if(abs(time - self%time) > drfts_eps_t) then
       self%it = self%it + 1
       self%time = time
    endif

    start1(1) = self%it_id
    len1(1)   = np

    it_indices = (/(self%it,i=1,np)/)
    ier = nf_put_vara_int( self%ncid, self%nc_index_time, start1, len1, it_indices )
    if(ier/=NF_NOERR) &
         & ermesg = 'drifters_io_write::failed to write index_time: ' //nf_strerror(ier)

    time_array = (/(time,i=1,np)/)
    ier = nf_put_vara_double( self%ncid, self%nc_time, start1, len1, time_array )
    if(ier/=NF_NOERR) &
         & ermesg = 'drifters_io_write::failed to write time: ' //nf_strerror(ier)

    ier = nf_put_vara_int(self%ncid, self%nc_ids, start1, len1, ids)
    if(ier/=NF_NOERR) &
         & ermesg = 'drifters_io_write::failed to write ids: '//nf_strerror(ier)

    start2(1) = 1
    start2(2) = self%it_id

    len2(1)   = nd
    len2(2)   = np

    ier = nf_put_vara_double(self%ncid, self%nc_positions, start2, len2, positions)
    if(ier/=NF_NOERR) &
         & ermesg = 'drifters_io_write::failed to write positions: '//nf_strerror(ier)

    len2(1)   = nf
    len2(2)   = np    
    
    ier = nf_put_vara_double(self%ncid, self%nc_fields, start2, len2, fields)
    if(ier/=NF_NOERR) &
         & ermesg = 'drifters_io_write::failed to write fields: '//nf_strerror(ier)

    self%it_id = self%it_id + np
    
  end subroutine drifters_io_write

end module drifters_io_mod
 !###############################################################################
 !###############################################################################
#ifdef _TEST_DRIFTERS_IO
! set FC=pgf95
! set FOPTS='-r8 -g -Mdclchk -Minform=warn'
! set INCS='-I/usr/local/include'
! set LIBS='-L/usr/local/lib -lnetcdf'
! $FC $INCS $FOPTS -D_TEST_DRIFTERS_IO drifters_io.F90 $LIBS
program test
  use drifters_io_mod
  implicit none
  type(drifters_io_type) :: drfts_io
  character(len=128) :: ermesg
  character(len=31) :: filename
  integer :: np, nd, nf, nt, i, j, k, npmax
  real :: dt, time, xmin, xmax, ymin, ymax, u, v, dr, x, y
  integer, allocatable :: ids(:)
  real, allocatable :: positions(:,:), fields(:,:)

  ! number of dimensions
  nd = 3
  ! number of fields 
  nf = 2
  ! max number of dirfters 
  npmax = 20
  ! number of time steps
  nt = 50
  ! starting time
  time = 0.

  ! domain boundary. (drifters outside domain will not be written to file.)
  xmin = 0.
  ymin = 0.
  xmax = 1.
  ymax = 1.

  ! constant velocity
  u = (xmax-xmin)*sqrt(2.)
  v = (ymax-ymin)*sqrt(2.)
  dt = 1/real(nt)

  ! open file
  
  filename = 'test.nc'
  call drifters_io_new(drfts_io, filename, nd, nf, ermesg)
  if(ermesg/='') print *,'ERROR after drifters_io_new: ', ermesg

  ! set attributes

  call drifters_io_set_position_names(drfts_io, (/'x','y','z'/), ermesg)
  if(ermesg/='') print *,'ERROR after drifters_io_position_names: ', ermesg

  ! note the trailing blanks in the first field, which are added here to 
  ! ensure that "salinity" will not be truncated (all names must have the 
  ! same length)
  call drifters_io_set_field_names(drfts_io, (/'temp    ','salinity'/), ermesg)
  if(ermesg/='') print *,'ERROR after drifters_io_field_names: ', ermesg

  call drifters_io_set_position_units(drfts_io, (/'deg east ','deg north','meters'/), ermesg)
  if(ermesg/='') print *,'ERROR after drifters_io_position_units: ', ermesg
  
  call drifters_io_set_field_units(drfts_io, (/'deg K ','ppm'/), ermesg)
  if(ermesg/='') print *,'ERROR after drifters_io_field_units: ', ermesg

  allocate(positions(nd, npmax), ids(npmax), fields(nf, npmax))
  dr = sqrt( (xmax-xmin)**2 + (ymax-ymin)**2 )/real(npmax)


  ! x
  positions(1, :) = +(/ (i*dr,i=0,npmax-1) /)/sqrt(2.)
  ! y
  positions(2, :) = -(/ (i*dr,i=0,npmax-1) /)/sqrt(2.)
  ! z
  positions(3, :) = 0.

  ! drifters' identity array (can be any integer number)
  ids = (/ (i, i=1, npmax) /)
  
  ! set fields as a function of space time
  fields(1, :) = sqrt( (positions(1,:)-xmin)**2 + (positions(2,:)-ymin)**2 )
  fields(2, :) = positions(1,:)-u*time + positions(2,:)-v*time ! invariant

  ! write to disk only drifters inside domain
  do i = 1, npmax
     x = positions(1,i)
     y = positions(2,i)
     if(x>=xmin .and. x<=xmax .and. y>=ymin .and. y<=ymax) then
        call drifters_io_write(drfts_io, time, np=1, nd=nd, nf=nf, &
             & ids=ids(i), positions=positions(:,i), fields=fields(:,i), ermesg=ermesg)
        if(ermesg/='') print *,'ERROR after drifters_io_write: ', ermesg
     endif
  enddo

  ! advect
  
  do j = 1, nt
     time = time + dt
     positions(1, :) = positions(1, :) + u*dt
     positions(2, :) = positions(2, :) + v*dt
     fields(1, :) = sqrt( (positions(1,:)-xmin)**2 + (positions(2,:)-ymin)**2 )
     fields(2, :) = positions(1,:)-u*time + positions(2,:)-v*time ! invariant

     do i = 1, npmax
        x = positions(1,i)
        y = positions(2,i)
        if(x>=xmin .and. x<=xmax .and. y>=ymin .and. y<=ymax) then
           call drifters_io_write(drfts_io, time, np=1, nd=nd, nf=nf, &
                & ids=ids(i), positions=positions(:,i), fields=fields(:,i), ermesg=ermesg)
           if(ermesg/='') print *,'ERROR after drifters_io_write: ', ermesg
        endif
     enddo
     
  enddo

  deallocate(positions, ids, fields)

  call drifters_io_del(drfts_io, ermesg)
  if(ermesg/='') print *,'ERROR after drifters_io_del: ', ermesg

end program test
#endif 
! _TEST_DRIFTERS_IO
