#include <fms_platform.h>


module drifters_input_mod
  implicit none
  private

  public :: drifters_input_type, drifters_input_new, drifters_input_del, drifters_input_save, assignment(=)

  ! Globals
  integer, parameter, private   :: MAX_STR_LEN = 128
  ! Include variable "version" to be written to log file.
#include<file_version.h>
  character, parameter, private :: SEPARATOR = ' '

  type drifters_input_type
     ! Be sure to update drifters_input_new, drifters_input_del and drifters_input_copy_new
     ! when adding members
     character(len=MAX_STR_LEN), _ALLOCATABLE :: position_names(:) _NULL
     character(len=MAX_STR_LEN), _ALLOCATABLE :: position_units(:) _NULL
     character(len=MAX_STR_LEN), _ALLOCATABLE :: field_names(:)    _NULL
     character(len=MAX_STR_LEN), _ALLOCATABLE :: field_units(:)    _NULL
     character(len=MAX_STR_LEN), _ALLOCATABLE :: velocity_names(:) _NULL
     real                      , _ALLOCATABLE :: positions(:,:) _NULL
     integer                   , _ALLOCATABLE :: ids(:)         _NULL
     character(len=MAX_STR_LEN)               :: time_units
     character(len=MAX_STR_LEN)               :: title
     character(len=MAX_STR_LEN)               :: version
  end type drifters_input_type

  interface assignment(=)
     module procedure drifters_input_copy_new
  end interface
 

  contains

!===============================================================================

  subroutine drifters_input_new(self, filename, ermesg)
    type(drifters_input_type)    :: self
    character(len=*), intent(in) :: filename
    character(len=*), intent(out):: ermesg

    ! Local 
    integer :: ier, ncid, nd, nf, np, ipos, j, id, i, isz
    character(len=MAX_STR_LEN) :: attribute
    include 'netcdf.inc'
    
    ermesg = ''

    ier = nf_open(filename, NF_NOWRITE, ncid)
    if(ier/=NF_NOERR) then
       ermesg = 'drifters_input: ERROR could not open netcdf file '//filename
       return
    endif

    ! version
    ier = NF_PUT_ATT_TEXT(NCID, NF_GLOBAL, 'version', len(version), version)

    ier = NF_INQ_DIMID(NCID, 'nd', id)
    if(ier/=NF_NOERR) then
       ermesg = 'drifters_input: ERROR could not find "nd" (number of dimensions)'
       ier = nf_close(ncid)
       return
    endif
    ier = NF_INQ_DIMLEN(NCID, id, nd)

    ! determine number of fields (nf)
    attribute = ''
    ier = nf_get_att_text(ncid, NF_GLOBAL, 'field_names', attribute)
    isz = min(len(attribute), len(trim(attribute))+1)
    attribute(isz:isz) = ' '
    ipos = 1
    nf = 0
    do i = 1, isz
       if(attribute(i:i)==SEPARATOR) then
          nf = nf + 1
       endif
    enddo
      
    ier = NF_INQ_DIMID(NCID, 'np', id)
    if(ier/=NF_NOERR) then
       ermesg = 'drifters_input: ERROR could not find "np" (number of particles)'
       ier = nf_close(ncid)
       return
    endif
    ier = NF_INQ_DIMLEN(NCID, id, np)

    allocate(self%position_names(nd))
    allocate(self%position_units(nd))
    allocate(self%field_names(nf))
    allocate(self%field_units(nf))
    allocate(self%velocity_names(nd))
    allocate(self%ids(np))
    allocate(self%positions(nd, np))

    ier = NF_INQ_VARID(NCID, 'ids', id)
    if(ier/=NF_NOERR) then
       ermesg = 'drifters_input: ERROR could not find "ids"'
       ier = nf_close(ncid)
       return
    endif
    ier = NF_GET_VAR_INT(NCID, id, self%ids)

    ier = NF_INQ_VARID(NCID, 'positions', id)
    if(ier/=NF_NOERR) then
       ermesg = 'drifters_input: ERROR could not find "positions"'
       ier = nf_close(ncid)
       return
    endif
    ier = NF_GET_VAR_DOUBLE(NCID, id, self%positions)

    attribute = ''
    ier = nf_get_att_text(ncid, NF_GLOBAL, 'version', attribute)
    self%version = trim(attribute)
   
    attribute = ''
    ier = nf_get_att_text(ncid, NF_GLOBAL, 'time_units', attribute)
    self%time_units = trim(attribute)

    attribute = ''
    ier = nf_get_att_text(ncid, NF_GLOBAL, 'title', attribute)
    self%title = trim(attribute)

    attribute = ''
    ier = nf_get_att_text(ncid, id, 'names', attribute)
    isz = min(len(attribute), len(trim(attribute))+1)
    attribute(isz:isz) = ' '
    ipos = 1
    j = 1
    do i = 1, isz
       if(attribute(i:i)==SEPARATOR) then
          self%position_names(j)  = trim(adjustl(attribute(ipos:i-1)))
          ipos = i+1
          j = j + 1
          if(j > nd) exit
       endif
    enddo

    attribute = ''
    ier = nf_get_att_text(ncid, id, 'units', attribute)
    isz = min(len(attribute), len(trim(attribute))+1)
    attribute(isz:isz) = ' '
    ipos = 1
    j = 1
    do i = 1, isz
       if(attribute(i:i)==SEPARATOR) then
          self%position_units(j)  = trim(adjustl(attribute(ipos:i-1)))
          ipos = i+1
          j = j + 1
          if(j > nd) exit
       endif
    enddo

    attribute = ''
    ier = nf_get_att_text(ncid, NF_GLOBAL, 'field_names', attribute)
    isz = min(len(attribute), len(trim(attribute))+1)
    attribute(isz:isz) = ' '
    ipos = 1
    j = 1
    do i = 1, isz
       if(attribute(i:i)==SEPARATOR) then
          self%field_names(j)  = trim(adjustl(attribute(ipos:i-1)))
          ipos = i+1
          j = j + 1
          if(j > nf) exit
       endif
    enddo

    attribute = ''
    ier = nf_get_att_text(ncid, NF_GLOBAL, 'field_units', attribute)
    isz = min(len(attribute), len(trim(attribute))+1)
    attribute(isz:isz) = ' '
    ipos = 1
    j = 1
    do i = 1, isz
       if(attribute(i:i)==SEPARATOR) then
          self%field_units(j)  = trim(adjustl(attribute(ipos:i-1)))
          ipos = i+1
          j = j + 1
          if(j > nf) exit
       endif
    enddo

    attribute = ''
    ier = nf_get_att_text(ncid, NF_GLOBAL, 'velocity_names', attribute)
    isz = min(len(attribute), len(trim(attribute))+1)
    attribute(isz:isz) = ' '
    ipos = 1
    j = 1
    do i = 1, isz
       if(attribute(i:i)==SEPARATOR) then
          self%velocity_names(j)  = trim(adjustl(attribute(ipos:i-1)))
          ipos = i+1
          j = j + 1
          if(j > nd) exit
       endif
    enddo

  end subroutine drifters_input_new

!===============================================================================
  subroutine drifters_input_del(self, ermesg)
    type(drifters_input_type)    :: self
    character(len=*), intent(out):: ermesg

    integer :: iflag

    ermesg = ''

    deallocate(self%position_names, stat=iflag)
    deallocate(self%position_units, stat=iflag)
    deallocate(self%field_names, stat=iflag)
    deallocate(self%field_units, stat=iflag)
    deallocate(self%velocity_names, stat=iflag)
    deallocate(self%ids, stat=iflag)
    deallocate(self%positions, stat=iflag)
    
  end subroutine drifters_input_del

!===============================================================================
  subroutine drifters_input_copy_new(new_instance, old_instance)

    type(drifters_input_type), intent(inout) :: new_instance
    type(drifters_input_type), intent(in)    :: old_instance

    allocate(new_instance%position_names( size(old_instance%position_names) ))
    allocate(new_instance%position_units( size(old_instance%position_units) ))
    allocate(new_instance%field_names( size(old_instance%field_names) ))
    allocate(new_instance%field_units( size(old_instance%field_units) ))
    allocate(new_instance%velocity_names( size(old_instance%velocity_names) ))
    new_instance%position_names = old_instance%position_names
    new_instance%position_units = old_instance%position_units 
    new_instance%field_names    = old_instance%field_names
    new_instance%field_units    = old_instance%field_units
    new_instance%velocity_names = old_instance%velocity_names
    new_instance%time_units     = old_instance%time_units
    new_instance%title          = old_instance%title
    new_instance%version        = old_instance%version
    allocate(new_instance%positions( size(old_instance%positions,1),size(old_instance%positions,2) ))
    new_instance%positions      = old_instance%positions
    allocate(new_instance%ids(size(old_instance%ids)))
    new_instance%ids            = old_instance%ids

  end subroutine drifters_input_copy_new

!===============================================================================
  subroutine drifters_input_save(self, filename, geolon, geolat, ermesg)
    ! save state in netcdf file. can be used as restart file.
    type(drifters_input_type)    :: self
    character(len=*), intent(in ):: filename
    real, intent(in), optional   :: geolon(:), geolat(:)
    character(len=*), intent(out):: ermesg

    integer ncid, nc_nd, nc_np, ier, nd, np, nf, nc_pos, nc_ids, i, j, n
    integer nc_lon, nc_lat
    character(len=MAX_STR_LEN) :: att

    include 'netcdf.inc'

    ermesg = ''
    
    ier = nf_create(filename, NF_CLOBBER, ncid)
    if(ier/=NF_NOERR) then 
       ermesg = 'drifters_input: ERROR cannot create '//filename
       return
    endif

    nd = size(self%positions, 1)
    np = size(self%positions, 2)
    nf = size(self%field_names)

    ! dimensions
    ier = nf_def_dim(ncid, 'nd', nd, nc_nd)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR creating dim "nd" '//nf_strerror(ier)

    ier = nf_def_dim(ncid, 'np', np, nc_np)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR creating dim "np" '//nf_strerror(ier)

    ! global attributes
    ier = nf_put_att_text(ncid, NF_GLOBAL, 'title', len_trim(self%title), self%title)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR setting global att "title" ' &
         & //nf_strerror(ier)
    
    ier = nf_put_att_text(ncid, NF_GLOBAL, 'time_units', len_trim(self%time_units), self%time_units)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR setting global att "time_units" ' &
         & //nf_strerror(ier)

    att = ''
    j = 1
    do i = 1, nf
       n = len_trim(self%field_units(i))
       att(j:j+n+1) = trim(self%field_units(i)) // ' '
       j = j + n + 1
    enddo
    ier = nf_put_att_text(ncid, NF_GLOBAL, 'field_units',   len_trim(att), &
         & att)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR setting global att "field_units" ' &
         & //nf_strerror(ier)

    att = ''
    j = 1
    do i = 1, nf
       n = len_trim(self%field_names(i))
       att(j:j+n+1) = trim(self%field_names(i)) // ' '
       j = j + n + 1
    enddo
    ier = nf_put_att_text(ncid, NF_GLOBAL, 'field_names',   len_trim(att), &
         & att)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR setting global att "field_names" ' &
         & //nf_strerror(ier)

    att = ''
    j = 1
    do i = 1, nd
       n = len_trim(self%velocity_names(i))
       att(j:j+n+1) = trim(self%velocity_names(i)) // ' '
       j = j + n + 1
    enddo
    ier = nf_put_att_text(ncid, NF_GLOBAL, 'velocity_names',   len_trim(att), &
         & att)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR setting global att "velocity_names" ' &
         & //nf_strerror(ier)

    ! variables
    ier = nf_def_var(ncid, 'positions', NF_DOUBLE, 2, (/nc_nd, nc_np/), nc_pos)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR creating var "positions" '//nf_strerror(ier)

    ier = nf_def_var(ncid, 'ids', NF_INT, 1, (/nc_np/), nc_ids)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR creating var "ids" '//nf_strerror(ier)

    ! optional: longitudes/latitudes in deg
    if(present(geolon)) then
       ier = nf_def_var(ncid, 'longitude', NF_DOUBLE, 1, (/nc_np/), nc_lon)
       if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR creating var "longitude" ' &
            & //nf_strerror(ier)
       att = 'degrees_east'
       ier = nf_put_att_text(ncid, nc_lon, 'units', len(trim(att)), trim(att))
       if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR setting att "units" to "longitude" ' &
         & //nf_strerror(ier)
    endif
    if(present(geolat)) then
       ier = nf_def_var(ncid, 'latitude', NF_DOUBLE, 1, (/nc_np/), nc_lat)
       if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR creating var "latitude" ' &
            & //nf_strerror(ier)
       att = 'degrees_north'
       ier = nf_put_att_text(ncid, nc_lat, 'units', len(trim(att)), trim(att))
       if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR setting att "units" to "latitude" ' &
         & //nf_strerror(ier)
    endif
    
    ! variable attributes
    
    att = ''
    j = 1
    do i = 1, nd
       n = len_trim(self%position_units(i))
       att(j:j+n+1) = trim(self%position_units(i)) // ' '
       j = j + n + 1
    enddo
    ier = nf_put_att_text(ncid, nc_pos, 'units',   len_trim(att), &
         & att)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR setting att "units" to "positions" ' &
         & //nf_strerror(ier)
    
    att = ''
    j = 1
    do i = 1, nd
       n = len_trim(self%position_names(i))
       att(j:j+n+1) = trim(self%position_names(i)) // ' '
       j = j + n + 1
    enddo
    ier = nf_put_att_text(ncid, nc_pos, 'names',   len_trim(att), &
         & att)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR setting att "names" to "positions" ' &
         & //nf_strerror(ier)

    ! end of define mode
    ier = nf_enddef(ncid)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR could not end define mode ' &
         & //nf_strerror(ier)

    ! data
    ier = nf_put_var_double(ncid, nc_pos, self%positions)    
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR could not write "positions" ' &
         & //nf_strerror(ier)
    
    ier = nf_put_var_int(ncid, nc_ids, self%ids)    
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR could not write "ids" ' &
         & //nf_strerror(ier)

    if(present(geolon)) then
       ier = nf_put_var_double(ncid, nc_lon, geolon)    
       if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR could not write "geolon" ' &
            & //nf_strerror(ier)       
    endif
    if(present(geolat)) then
       ier = nf_put_var_double(ncid, nc_lat, geolat)    
       if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR could not write "geolat" ' &
            & //nf_strerror(ier)       
    endif
    
    
    ier = nf_close(ncid)
    if(ier/=NF_NOERR) ermesg = 'drifters_input_save: ERROR could not close file ' &
         & //nf_strerror(ier)
    
  end subroutine drifters_input_save

end module drifters_input_mod

!===============================================================================
!===============================================================================
#ifdef _TEST_DRIFTERS_INPUT
program test
  use drifters_input_mod
  implicit none
  character(len=128) :: ermesg
  integer :: i
  
  type(drifters_input_type) :: obj
  
  call drifters_input_new(obj, 'input.nc', ermesg)
  if(ermesg/='') print *,'ERROR: ', ermesg

  print *,'field_names:'
  do i = 1, size(obj%field_names)
     print *,trim(obj%field_names(i))
  enddo

  print *,'velocity_names:'
  do i = 1, size(obj%velocity_names)
     print *,trim(obj%velocity_names(i))
  enddo

  print *,'ids = ', obj%ids

  print *,'positions: '
  do i = 1, size(obj%positions, 2)
     print *,obj%positions(:,i)
  enddo

  call drifters_input_del(obj, ermesg)
end program test

#endif
