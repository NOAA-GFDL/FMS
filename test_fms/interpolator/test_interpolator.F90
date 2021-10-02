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

!***********************************************************************
!* Requires namelist but it can be empty. Requires a diag_table such as
!* below, requires the input files: (aerosol.climatology.nc and
!* o3.climatology.nc)

!* diag_table
!* test_diag_manager_01
!* 1 3 1 0 0 0

!* #output files
!*  "diag_test_01",  1, "days", 1, "days", "time"

!* #output variables
!* "test_diag_manager_mod", "dat1", "dat1", "diag_test_01",  "all", .false., "none", 2
!***********************************************************************

program test_interpolator

use mpp_mod, only: input_nml_file

use mpp_mod
use mpp_domains_mod
use fms_mod
use time_manager_mod
use diag_manager_mod!, only : diag_axis_init, file_exist, MPP_NPES, &
                    !  MPP_PE, REGISTER_DIAG_FIELD, SEND_DATA, SET_DATE,&
                    !  SET_TIME

use interpolator_mod
!use sulfate_mod
!use ozone_mod
use constants_mod !, only : grav, constants_init, PI
use time_interp_mod, only : time_interp_init

implicit none
integer, parameter :: nsteps_per_day = 8, ndays = 16
real, parameter :: delt = 1.0/nsteps_per_day
! integer, parameter :: nxd = 144, nyd = 90, ntsteps = 240, two_delt = 2*delt
integer, parameter :: nxd = 20, nyd = 40, ntsteps = nsteps_per_day*ndays, two_delt = 2*delt
integer :: delt_days, delt_secs
integer, parameter :: max_fields = 20 ! maximum number of fields to be interpolated

integer :: i,k,n,level
integer :: io_status
integer :: ndivs
integer :: jscomp, jecomp, iscomp, iecomp, isd,ied,jsd,jed
integer :: numfields, domain_layout(2)
integer :: num_nbrs, nbins,axes(3), interp_diagnostic_id
integer :: column_diagnostic_id1, column_diagnostic_id(max_fields)
integer :: io, ierr

real ::  missing_value = -1.e10

character(len=1) :: dest_grid
character(len=128) :: src_file, file_out, title, units, colaer
logical :: vector_field=.false., result

type(domain2d) :: domain
type(time_type) :: model_time

type(interpolate_type) :: o3, aerosol

real, dimension(:,:), allocatable :: col_data
real, dimension(:,:,:), allocatable :: model_data, p_half, p_full
real, dimension(:), allocatable :: latb_mod(:,:),lonb_mod(:,:),lon_mod(:),lat_mod(:)
real :: dx,dy
real :: dtr,tpi
real :: p_bot,p_top,lambda
character(len=64) :: names(13)
data names(:) /"so4_anthro","so4_natural","organic_carbon","black_carbon","sea_salt",&
"anthro_dust_0.2","anthro_dust_0.8","anthro_dust_2.0","anthro_dust_8.0",&
"natural_dust_0.2","natural_dust_0.8","natural_dust_2.0","natural_dust_8.0"/

integer :: out_of_bounds(1)
data out_of_bounds / CONSTANT/!, CONSTANT/!, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT, &
!ZERO, ZERO, ZERO, ZERO /

namelist /interpolator_nml/ src_file

! initialize communication modules

delt_days = INT(delt)
delt_secs = INT(delt*SECONDS_PER_DAY) - delt_days*SECONDS_PER_DAY

write(*,*) delt, delt_days,delt_secs

call mpp_init
call mpp_domains_init
call set_calendar_type(JULIAN)
call diag_manager_init
call constants_init
call time_interp_init

level = 18
tpi = 2.0*PI !4.*acos(0.)
dtr = tpi/360.

src_file = 'src_file'  ! input file containing fields to be interpolated


model_time = set_date(1979,12,1,0,0,0)

!if (numfields.ne.2.and.vector_field) call mpp_error(FATAL,'2 components of vector field not specified')
!if (numfields.gt.1.and..not.vector_field) call mpp_error(FATAL,'only 1 scalar at a time')
!if (numfields .gt. max_fields) call mpp_error(FATAL,'max num fields exceeded')

!--------------------------------------------------------------------
! namelist input
!--------------------------------------------------------------------

  read (input_nml_file, nml=interpolator_nml, iostat=io)
  ierr = check_nml_error(io, 'interpolator_nml')

! decompose model grid points
! mapping can get expensive so we distribute the task at this level

ndivs = mpp_npes()

call mpp_define_layout ((/1,nxd,1,nyd/), ndivs, domain_layout)
call mpp_define_domains((/1,nxd,1,nyd/),domain_layout, domain,xhalo=0,yhalo=0)
call mpp_get_data_domain(domain,isd,ied,jsd,jed)
call mpp_get_compute_domain (domain, iscomp, iecomp, jscomp, jecomp)

allocate(lonb_mod(nxd+1,nyd+1),lon_mod(nxd))
allocate(latb_mod(nxd+1,nyd+1),lat_mod(nyd))
allocate(col_data(isd:ied,jsd:jed)) ; col_data = 0.0
allocate(p_half(isd:ied,jsd:jed,level+1),p_full(isd:ied,jsd:jed,level))
p_top = 1.0
p_bot = 101325.0 !Model level in Pa
lambda = -1.0*log(p_top/p_bot)/(level+1)

p_half(:,:,level+1) = p_bot
do i=level,1,-1
  p_half(:,:,i)=p_half(:,:,i+1)*exp(-1.0*lambda)
enddo
do i=1,level
  p_full(:,:,i)=(p_half(:,:,i+1)+p_half(:,:,i))/2.0
enddo

allocate(model_data(isd:ied,jsd:jed,level))

dx = 360./nxd
dy = 180./nyd
do i = 1,nxd+1
  lonb_mod(i,:) = (i-1)*dx
enddo
do i = 1,nyd+1
  latb_mod(:,i) = -90. + (i-1)*dy
enddo
do i=1,nxd
  lon_mod(i)=(lonb_mod(i+1,1)+lonb_mod(i,1))/2.0
enddo
do i=1,nyd
  lat_mod(i)=(latb_mod(1,i+1)+latb_mod(1,i))/2.0
enddo

lonb_mod = lonb_mod * dtr
latb_mod = latb_mod * dtr

   axes(1) = diag_axis_init('x',lon_mod,units='degrees',cart_name='x',domain2=domain)
   axes(2) = diag_axis_init('y',lat_mod,units='degrees',cart_name='y',domain2=domain)
   axes(3) = diag_axis_init('z',p_full(isd,jsd,:),units='mb',cart_name='z')

interp_diagnostic_id =  register_diag_field('interp','ozone',axes(1:3),model_time,&
                                'interpolated_ozone_clim', 'kg/kg', missing_value)
column_diagnostic_id1 =  register_diag_field('interp','colozone',axes(1:2),model_time,&
                                'column_ozone_clim', 'kg/m2', missing_value)

do i=1,size(names(:))
colaer = 'col'//trim(names(i))
column_diagnostic_id(i) =  register_diag_field('interp',colaer,axes(1:2),model_time,&
                                'column_aerosol_clim', 'kg/m2', missing_value)
enddo

call ozone_init(o3,lonb_mod(isd:ied+1,jsd:jed+1), latb_mod(isd:ied+1,jsd:jed+1), axes, model_time, &
                data_out_of_bounds=out_of_bounds)

call init_clim_diag(o3, axes, model_time)
call sulfate_init(aerosol,lonb_mod(isd:ied+1,jsd:jed+1), latb_mod(isd:ied+1,jsd:jed+1), names, &
                data_out_of_bounds=(/CONSTANT/) )
call init_clim_diag(aerosol, axes, model_time)

do n=1,ntsteps
  if( mpp_pe() == mpp_root_pe() ) write(*,*) n

  call get_ozone(o3,model_time,p_half,model_data)

  if(interp_diagnostic_id>0) &
       result = send_data(interp_diagnostic_id,&
            model_data(iscomp:iecomp,jscomp:jecomp,:),model_time)

  if(column_diagnostic_id1>0) then

    col_data(iscomp:iecomp,jscomp:jecomp)=0.0
    do k=1,level
       col_data(iscomp:iecomp,jscomp:jecomp)= col_data(iscomp:iecomp,jscomp:jecomp)+ &
          model_data(iscomp:iecomp,jscomp:jecomp,k)* &
          (p_half(iscomp:iecomp,jscomp:jecomp,k+1)-p_half(iscomp:iecomp,jscomp:jecomp,k))/grav
    enddo
       result = send_data(column_diagnostic_id1,col_data(:,:),model_time)
  endif



  do i=1,size(names(:))

call get_anthro_sulfate(aerosol,model_time,p_half,names(i),model_data,clim_units=units)

    if(column_diagnostic_id(i)>0) then

      col_data(iscomp:iecomp,jscomp:jecomp)=0.0
      do k=1,level
        if (trim(adjustl(lowercase(units))) .eq. 'kg/m^2') then
           col_data(iscomp:iecomp,jscomp:jecomp)= col_data(iscomp:iecomp,jscomp:jecomp)+ &
              model_data(iscomp:iecomp,jscomp:jecomp,k)
        else
           col_data(iscomp:iecomp,jscomp:jecomp)= col_data(iscomp:iecomp,jscomp:jecomp)+ &
              model_data(iscomp:iecomp,jscomp:jecomp,k)* &
              (p_half(iscomp:iecomp,jscomp:jecomp,k+1)-p_half(iscomp:iecomp,jscomp:jecomp,k))/grav
        endif
      enddo
      result = send_data(column_diagnostic_id(i),&
      col_data(iscomp:iecomp,jscomp:jecomp),model_time)
    endif

  enddo

   model_time = model_time + set_time(delt_secs,delt_days)

   if (n.eq. ntsteps) call diag_manager_end(model_time)

enddo

call interpolator_end(aerosol)
call interpolator_end(o3)

deallocate(lonb_mod, lon_mod, latb_mod,lat_mod, col_data, p_half, p_full, model_data)

call mpp_exit

contains
!
!#######################################################################
!
subroutine sulfate_init(aerosol,lonb, latb, names, data_out_of_bounds, vert_interp, units)
type(interpolate_type), intent(inout)         :: aerosol
real,                   intent(in)            :: lonb(:,:),latb(:,:)
character(len=64),      intent(in)            :: names(:)
integer,                intent(in)            :: data_out_of_bounds(:)
integer,                intent(in), optional  :: vert_interp(:)
character(len=*),       intent(out),optional  :: units(:)

if (.not. file_exist("INPUT/aerosol.climatology.nc") ) return
call interpolator_init( aerosol, "aerosol.climatology.nc", lonb, latb, &
                        data_names=names, data_out_of_bounds=data_out_of_bounds, &
                        vert_interp=vert_interp, clim_units=units )

end subroutine sulfate_init
!
!#######################################################################
!
subroutine get_anthro_sulfate( sulfate, model_time, p_half, name, model_data, is, js, clim_units )
type(interpolate_type), intent(inout) :: sulfate
type(time_type), intent(in) :: model_time
real, intent(in)           :: p_half(:,:,:)
character(len=*), intent(in) :: name
character(len=*), intent(out), optional :: clim_units
real, intent(out) :: model_data(:,:,:)
integer, intent(in), optional :: is,js

call interpolator( sulfate, model_time, p_half, model_data, name, is, js, clim_units)

end subroutine get_anthro_sulfate
!
!#######################################################################
!
subroutine ozone_init( o3, lonb, latb, axes, model_time, data_out_of_bounds, vert_interp )
real,                  intent(in)           :: lonb(:,:),latb(:,:)
integer,               intent(in)           :: axes(:)
type(time_type),       intent(in)           :: model_time
type(interpolate_type),intent(inout)        :: o3
integer,               intent(in)           :: data_out_of_bounds(:)
integer,               intent(in), optional :: vert_interp(:)

if (.not. file_exist("INPUT/o3.climatology.nc") ) return
call interpolator_init( o3, "o3.climatology.nc", lonb, latb, &
                        data_out_of_bounds=data_out_of_bounds, vert_interp=vert_interp )

end subroutine ozone_init
!
!#######################################################################
!
subroutine get_ozone( o3, model_time, p_half, model_data, is, js )
type(interpolate_type),intent(inout) :: o3
type(time_type), intent(in) :: model_time
real, intent(in)           :: p_half(:,:,:)
real, intent(out) :: model_data(:,:,:)
integer, intent(in), optional :: is,js

call interpolator( o3, model_time, p_half, model_data, "ozone", is, js)

end subroutine get_ozone

end program test_interpolator
