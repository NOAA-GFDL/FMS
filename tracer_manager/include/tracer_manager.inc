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
!> @defgroup tracer_manager_mod tracer_manager_mod
!> @ingroup tracer_manager


!> @brief Subroutine to set the tracer field to the wanted profile.
!!
!> If the profile type is 'fixed' then the tracer field values are set
!! equal to the surface value.
!! If the profile type is 'profile' then the top/bottom of model and
!! surface values are read and an exponential profile is calculated,
!! with the profile being dependent on the number of levels in the
!! component model. This should be called from the part of the dynamical
!! core where tracer restarts are called in the event that a tracer
!! restart file does not exist.
!!
!!  This can be activated by adding a method to the field_table
!!  e.g.
!!  @verbose "profile_type","fixed","surface_value = 1e-12" @endverbose
!!  would return values of surf_value = 1e-12 and a multiplier of 1.0
!!  One can use these to initialize the entire field with a value of 1e-12.
!!
!!  "profile_type","profile","surface_value = 1e-12, top_value = 1e-15"
!!   In a 15 layer model this would return values of surf_value = 1e-12 and
!!   multiplier = 0.6309573 i.e 1e-15 = 1e-12*(0.6309573^15)
!!   In this case the model should be MODEL_ATMOS as you have a "top" value.
!!
!!   If you wish to initialize the ocean model, one can use bottom_value instead
!!   of top_value.
subroutine SET_TRACER_PROFILE_(model, n, tracer, err_msg)

integer, intent(in) :: model !< Parameter representing component model in use
integer, intent(in) :: n !< Tracer number
real(FMS_TM_KIND_), intent(inout), dimension(:,:,:) :: tracer !< Initialized tracer array
character(len=*), intent(out), optional :: err_msg

real(FMS_TM_KIND_)  :: surf_value, multiplier
integer :: numlevels, k, n1, flag
real(FMS_TM_KIND_) :: top_value, bottom_value
character(len=80)  :: scheme, control,profile_type
character(len=128) :: err_msg_local
character(len=11)  :: chn

integer, parameter :: lkind=FMS_TM_KIND_

if(.not.module_is_initialized) call tracer_manager_init

if (n < 1 .or. n > total_tracers(model)) then
  write(chn, '(i11)') n
  err_msg_local = ' Invalid tracer index.  Model name = '//trim(MODEL_NAMES(model))//',  Index='//trim(chn)
  if(error_handler('set_tracer_profile', err_msg_local, err_msg)) return
endif
n1 = TRACER_ARRAY(model,n)

!default values
profile_type  = 'Fixed'
surf_value = 0.0E+00_lkind
top_value  = surf_value
bottom_value = surf_value
multiplier = 1.0_lkind

tracer = surf_value

if ( query_method ( 'profile_type',model,n,scheme,control)) then
!Change the tracer_number to the tracer_manager version

  if(lowercase(trim(scheme(1:5))).eq.'fixed') then
    profile_type                   = 'Fixed'
    flag =parse(control,'surface_value',surf_value)
    multiplier = 1.0_lkind
    tracer = surf_value
  endif

  if(lowercase(trim(scheme(1:7))).eq.'profile') then
    profile_type                   = 'Profile'
    flag=parse(control,'surface_value',surf_value)
    if (surf_value .eq. 0.0_lkind) &
      call mpp_error(FATAL,'set_tracer_profile : Cannot have a zero surface value for an exponential profile. Tracer '&
                           //tracers(n1)%tracer_name//" "//control//" "//scheme)
    select case (tracers(n1)%model)
      case (MODEL_ATMOS)
        flag=parse(control,'top_value',top_value)
        if(mpp_pe()==mpp_root_pe() .and. flag == 0) &
           call mpp_error(NOTE,'set_tracer_profile : Parameter top_value needs to be defined for the tracer profile.')
      case (MODEL_OCEAN)
        flag =parse(control,'bottom_value',bottom_value)
        if(mpp_pe() == mpp_root_pe() .and. flag == 0) &
           call mpp_error(NOTE, &
                          & 'set_tracer_profile : Parameter bottom_value needs to be defined for the tracer profile.')
      case default
!   Should there be a NOTE or WARNING message here?
    end select

! If profile type is profile then set the surface value to the input
! value and calculate the vertical multiplier.
!
! Assume an exponential decay/increase from the surface to the top level
!  C = C0 exp ( -multiplier* level_number)
!  => multiplier = exp [ ln(Ctop/Csurf)/number_of_levels]
!
numlevels = size(tracer,3) -1
    select case (tracers(n1)%model)
      case (MODEL_ATMOS)
        multiplier = exp( log (top_value/surf_value) /real(numlevels,lkind))
        tracer(:,:,1) = surf_value
        do k = 2, size(tracer,3)
          tracer(:,:,k) = tracer(:,:,k-1) * multiplier
        enddo
      case (MODEL_OCEAN)
        multiplier = exp( log (bottom_value/surf_value) / real(numlevels,lkind))
        tracer(:,:,size(tracer,3)) = surf_value
        do k = size(tracer,3) - 1, 1, -1
          tracer(:,:,k) = tracer(:,:,k+1) * multiplier
        enddo
      case default
    end select
  endif !scheme.eq.profile

  if (mpp_pe() == mpp_root_pe() ) write(*,700) 'Tracer ',trim(tracers(n1)%tracer_name),    &
                            ' initialized with surface value of ',surf_value, &
                            ' and vertical multiplier of ',multiplier
  700 FORMAT (3A,E13.6,A,F13.6)

endif ! end of query scheme

end subroutine SET_TRACER_PROFILE_

!#######################################################################
!> @}
! close documentation grouping
