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
!> @addtogroup sat_vapor_pres_mod
!> @{

!#######################################################################
! <SUBROUTINE NAME="lookup_es_0d" INTERFACE="lookup_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES_0D_ ( temp, esat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es_k(temp, esat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES_0D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es_1d" INTERFACE="lookup_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES_1D_ ( temp, esat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp(:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat(:) !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 character(len=54) :: err_msg_local
 integer :: nbad !< if temperature is out of range
!-----------------------------------------------

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es_k(temp, esat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es',err_msg_local,err_msg)) return
   endif

!-----------------------------------------------

 end subroutine LOOKUP_ES_1D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es_2d" INTERFACE="lookup_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES_2D_ ( temp, esat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp(:,:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat(:,:) !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 character(len=54) :: err_msg_local
 integer :: nbad !< if temperature is out of range
!-----------------------------------------------

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es_k(temp, esat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es',err_msg_local,err_msg)) return
   endif

!-----------------------------------------------

 end subroutine LOOKUP_ES_2D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es_3d" INTERFACE="lookup_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES_3D_ ( temp, esat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp(:,:,:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat(:,:,:) !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es_k(temp, esat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es',err_msg_tmp,err_msg)) return
   endif

 end subroutine LOOKUP_ES_3D_


!#######################################################################
! <SUBROUTINE NAME="lookup_es2_0d" INTERFACE="lookup_es2">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES2_0D_ ( temp, esat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es2_k(temp, esat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es2',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES2_0D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es2_1d" INTERFACE="lookup_es2">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES2_1D_ ( temp, esat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp(:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat(:) !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 character(len=54) :: err_msg_local
 integer :: nbad !< if temperature is out of range
!-----------------------------------------------

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es2_k(temp, esat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es2',err_msg_local,err_msg)) return
   endif

!-----------------------------------------------

 end subroutine LOOKUP_ES2_1D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es2_2d" INTERFACE="lookup_es2">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES2_2D_ ( temp, esat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp(:,:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat(:,:) !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 character(len=54) :: err_msg_local
 integer :: nbad !< if temperature is out of range
!-----------------------------------------------

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es2_k(temp, esat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es2',err_msg_local,err_msg)) return
   endif

!-----------------------------------------------

 end subroutine LOOKUP_ES2_2D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es2_3d" INTERFACE="lookup_es2">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES2_3D_ ( temp, esat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp(:,:,:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat(:,:,:) !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es2_k(temp, esat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es2',err_msg_tmp,err_msg)) return
   endif

 end subroutine LOOKUP_ES2_3D_


!#######################################################################
! <SUBROUTINE NAME="lookup_es3_0d" INTERFACE="lookup_es3">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES3_0D_ ( temp, esat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es3_k(temp, esat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es3',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES3_0D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es3_1d" INTERFACE="lookup_es3">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES3_1D_ ( temp, esat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp(:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat(:) !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 character(len=54) :: err_msg_local
 integer :: nbad !< if temperature is out of range
!-----------------------------------------------

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es3_k(temp, esat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es3',err_msg_local,err_msg)) return
   endif

!-----------------------------------------------

 end subroutine LOOKUP_ES3_1D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es3_2d" INTERFACE="lookup_es3">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES3_2D_ ( temp, esat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp(:,:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat(:,:) !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 character(len=54) :: err_msg_local
 integer :: nbad !< if temperature is out of range
!-----------------------------------------------

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es3_k(temp, esat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es3',err_msg_local,err_msg)) return
   endif

!-----------------------------------------------

 end subroutine LOOKUP_ES3_2D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es3_3d" INTERFACE="lookup_es3">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES3_3D_ ( temp, esat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp(:,:,:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat(:,:,:) !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es3_k(temp, esat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es3',err_msg_tmp,err_msg)) return
   endif

 end subroutine LOOKUP_ES3_3D_


!#######################################################################
!  routines for computing derivative of es
!#######################################################################

! <SUBROUTINE NAME="lookup_des_0d" INTERFACE="lookup_des">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_DES_0D_ ( temp, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_des_k( temp, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_des',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_DES_0D_

!#######################################################################

! <SUBROUTINE NAME="lookup_des_1d" INTERFACE="lookup_des">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_DES_1D_ ( temp, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp (:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: desat(:) !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 character(len=54) :: err_msg_local
 integer :: nbad !< if temperature is out of range
!-----------------------------------------------

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   if(present(err_msg)) err_msg=''

   call lookup_des_k(temp, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es',err_msg_local,err_msg)) return
   endif
!-----------------------------------------------

 end subroutine LOOKUP_DES_1D_

!#######################################################################

! <SUBROUTINE NAME="lookup_des_2d" INTERFACE="lookup_des">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_DES_2D_ ( temp, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp (:,:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: desat(:,:) !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 character(len=54) :: err_msg_local
 integer :: nbad !< if temperature is out of range
!-----------------------------------------------

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_des_k(temp, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es',err_msg_local,err_msg)) return
   endif
!-----------------------------------------------

 end subroutine LOOKUP_DES_2D_

!#######################################################################
! <SUBROUTINE NAME="lookup_des_3d" INTERFACE="lookup_des">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_DES_3D_ ( temp, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp (:,:,:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: desat(:,:,:) !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_des_k( temp, desat, nbad )

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg=''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_des',err_msg_tmp,err_msg)) return
   endif

 end subroutine LOOKUP_DES_3D_


! <SUBROUTINE NAME="lookup_des2_0d" INTERFACE="lookup_des2">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_DES2_0D_ ( temp, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_des2_k( temp, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_des2',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_DES2_0D_

!#######################################################################

! <SUBROUTINE NAME="lookup_des2_1d" INTERFACE="lookup_des2">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_DES2_1D_ ( temp, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp (:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: desat(:) !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 character(len=54) :: err_msg_local
 integer :: nbad !< if temperature is out of range
!-----------------------------------------------

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   if(present(err_msg)) err_msg=''

   call lookup_des2_k(temp, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_des2',err_msg_local,err_msg)) return
   endif
!-----------------------------------------------

 end subroutine LOOKUP_DES2_1D_

!#######################################################################

! <SUBROUTINE NAME="lookup_des2_2d" INTERFACE="lookup_des2">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_DES2_2D_ ( temp, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp (:,:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: desat(:,:) !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 character(len=54) :: err_msg_local
 integer :: nbad !< if temperature is out of range
!-----------------------------------------------

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_des2_k(temp, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_des2',err_msg_local,err_msg)) return
   endif
!-----------------------------------------------

 end subroutine LOOKUP_DES2_2D_

!#######################################################################
! <SUBROUTINE NAME="lookup_des2_3d" INTERFACE="lookup_des2">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_DES2_3D_ ( temp, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp (:,:,:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: desat(:,:,:) !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_des2_k( temp, desat, nbad )

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg=''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_des2',err_msg_tmp,err_msg)) return
   endif

 end subroutine LOOKUP_DES2_3D_


! <SUBROUTINE NAME="lookup_des3_0d" INTERFACE="lookup_des3">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_DES3_0D_ ( temp, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_des3_k( temp, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_des3',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_DES3_0D_

!#######################################################################

! <SUBROUTINE NAME="lookup_des3_1d" INTERFACE="lookup_des3">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_DES3_1D_ ( temp, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp (:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: desat(:) !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 character(len=54) :: err_msg_local
 integer :: nbad !< if temperature is out of range
!-----------------------------------------------

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   if(present(err_msg)) err_msg=''

   call lookup_des3_k(temp, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_des3',err_msg_local,err_msg)) return
   endif
!-----------------------------------------------

 end subroutine LOOKUP_DES3_1D_

!#######################################################################

! <SUBROUTINE NAME="lookup_des3_2d" INTERFACE="lookup_des3">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_DES3_2D_ ( temp, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp (:,:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: desat(:,:) !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 character(len=54) :: err_msg_local
 integer :: nbad !< if temperature is out of range
!-----------------------------------------------

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_des3_k(temp, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_des3',err_msg_local,err_msg)) return
   endif
!-----------------------------------------------

 end subroutine LOOKUP_DES3_2D_

!#######################################################################
! <SUBROUTINE NAME="lookup_des3_3d" INTERFACE="lookup_des3">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_DES3_3D_ ( temp, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp (:,:,:) !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: desat(:,:,:) !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_des3_k( temp, desat, nbad )

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg=''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_des3',err_msg_tmp,err_msg)) return
   endif

 end subroutine LOOKUP_DES3_3D_

!========================================================================================================

!#######################################################################

! <SUBROUTINE NAME="lookup_es_des_0d" INTERFACE="lookup_es_des">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
!   <OUT NAME="desat" UNITS="pascal / degree" TYPE="real" DIM="(scalar)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES_DES_0D_ ( temp, esat, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es_des_k(temp, esat, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES_DES_0D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es_des_1d" INTERFACE="lookup_es_des">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="desat" UNITS="pascal / degree" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES_DES_1D_ ( temp, esat, desat, err_msg )

 real(kind=FMS_SVP_KIND_), dimension(:), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), dimension(:), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), dimension(:), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es_des_k(temp, esat, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES_DES_1D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es_des_2d" INTERFACE="lookup_es_des">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="desat" UNITS="pascal / degree" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES_DES_2D_ ( temp, esat, desat, err_msg )

 real(kind=FMS_SVP_KIND_), dimension(:,:), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), dimension(:,:), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), dimension(:,:), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es_des_k(temp, esat, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES_DES_2D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es_des_3d" INTERFACE="lookup_es_des">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="desat" UNITS="pascal / degree" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES_DES_3D_ ( temp, esat, desat, err_msg )

 real(kind=FMS_SVP_KIND_), dimension(:,:,:), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), dimension(:,:,:), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), dimension(:,:,:), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es_des_k(temp, esat, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES_DES_3D_

!#######################################################################
!#######################################################################

! <SUBROUTINE NAME="lookup_es2_des2_0d" INTERFACE="lookup_es2_des2">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
!   <OUT NAME="desat" UNITS="pascal / degree" TYPE="real" DIM="(scalar)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES2_DES2_0D_ ( temp, esat, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es2_des2_k(temp, esat, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es2_des2',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES2_DES2_0D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es2_des2_1d" INTERFACE="lookup_es2_des2">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="desat" UNITS="pascal / degree" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES2_DES2_1D_ ( temp, esat, desat, err_msg )

 real(kind=FMS_SVP_KIND_), dimension(:), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), dimension(:), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), dimension(:), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es2_des2_k(temp, esat, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es2_des2',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES2_DES2_1D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es2_des2_2d" INTERFACE="lookup_es2_des2">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="desat" UNITS="pascal / degree" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES2_DES2_2D_ ( temp, esat, desat, err_msg )

 real(kind=FMS_SVP_KIND_), dimension(:,:), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), dimension(:,:), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), dimension(:,:), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es2_des2_k(temp, esat, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es2_des2',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES2_DES2_2D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es2_des2_3d" INTERFACE="lookup_es2_des2">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="desat" UNITS="pascal / degree" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES2_DES2_3D_ ( temp, esat, desat, err_msg )

 real(kind=FMS_SVP_KIND_), dimension(:,:,:), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), dimension(:,:,:), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), dimension(:,:,:), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es2_des2_k(temp, esat, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es2_des2',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES2_DES2_3D_


!#######################################################################
!#######################################################################

! <SUBROUTINE NAME="lookup_es3_des3_0d" INTERFACE="lookup_es3_des3">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
!   <OUT NAME="desat" UNITS="pascal / degree" TYPE="real" DIM="(scalar)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES3_DES3_0D_ ( temp, esat, desat, err_msg )

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es3_des3_k(temp, esat, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es3_des3',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES3_DES3_0D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es3_des3_1d" INTERFACE="lookup_es3_des3">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="desat" UNITS="pascal / degree" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES3_DES3_1D_ ( temp, esat, desat, err_msg )

 real(kind=FMS_SVP_KIND_), dimension(:), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), dimension(:), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), dimension(:), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es3_des3_k(temp, esat, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es3_des3',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES3_DES3_1D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es3_des3_2d" INTERFACE="lookup_es3_des3">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="desat" UNITS="pascal / degree" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES3_DES3_2D_ ( temp, esat, desat, err_msg )

 real(kind=FMS_SVP_KIND_), dimension(:,:), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), dimension(:,:), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), dimension(:,:), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es3_des3_k(temp, esat, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es3_des3',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES3_DES3_2D_

!#######################################################################

! <SUBROUTINE NAME="lookup_es3_des3_3d" INTERFACE="lookup_es3_des3">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="desat" UNITS="pascal / degree" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine LOOKUP_ES3_DES3_3D_ ( temp, esat, desat, err_msg )

 real(kind=FMS_SVP_KIND_), dimension(:,:,:), intent(in)  :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), dimension(:,:,:), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), dimension(:,:,:), intent(out) :: desat !< derivative of saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_local

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   call lookup_es3_des3_k(temp, esat, desat, nbad)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_local,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('lookup_es3_des3',err_msg_local,err_msg)) return
   endif

 end subroutine LOOKUP_ES3_DES3_3D_

!#######################################################################

! <SUBROUTINE NAME="compute_qs_0d" INTERFACE="compute_qs">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(SCALAR)"></IN>
!   <IN NAME="press UNIT="Pascals" TYPE="real" DIM="(SCALAR)"></IN>
!   <OUT NAME="qsat" UNITS="kg(vapor)/kg(moist air)" TYPE="real" DIM="(SCALAR)"></OUT>
!   <IN NAME="q" UNIT="kg(vapor)/kg(moistair)" TYPE="real" DIM="(SCALAR)"></IN>
!   <IN NAME="hc" UNIT="fraction" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="dqsdT" UNIT="kg(vapor)/kg(moistair)/ degree Kelvin" TYPE="real" DIM="(SCALAR)"></OUT>
!   <OUT NAME="esat" UNITS="Pascals" TYPE="real" DIM="(scalar)"> </OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine COMPUTE_QS_0D_ ( temp, press, qsat, q, hc, dqsdT, esat, &
                            err_msg, es_over_liq, es_over_liq_and_ice )

 real(kind=FMS_SVP_KIND_), intent(in)    :: temp  !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(in)    :: press !< pressure
 real(kind=FMS_SVP_KIND_), intent(out)   :: qsat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: q     !< vapor relative humidity
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: hc    !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), optional :: dqsdT !< d(qsat)/dT
 real(kind=FMS_SVP_KIND_), intent(out), optional :: esat  !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg
 logical,intent(in),            optional :: es_over_liq
 logical,intent(in),            optional :: es_over_liq_and_ice

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp

 !> EPSILO and ZVIR are module level variables that are declared in r8_kind.
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(kind=FMS_SVP_KIND_), parameter :: EPSILO_loc=real(EPSILO,FMS_SVP_KIND_)
 real(kind=FMS_SVP_KIND_), parameter :: ZVIRl=real(ZVIR,FMS_SVP_KIND_)

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   if (present(es_over_liq)) then
     if (.not. (construct_table_wrt_liq)) then
       call error_mesg ('compute_qs', &
          'requesting es wrt liq, but that table not constructed', &
                                                                FATAL)
     endif
   endif
   if (present(es_over_liq_and_ice)) then
     if (.not. (construct_table_wrt_liq_and_ice)) then
       call error_mesg ('compute_qs', &
      'requesting es wrt liq and ice, but that table not constructed', &
                                                                FATAL)
     endif
   endif

   call compute_qs_k (temp, press,  EPSILO_loc, ZVIRl, qsat, nbad, q, hc, &
                       dqsdT, esat, es_over_liq, es_over_liq_and_ice)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('compute_qs',err_msg_tmp,err_msg)) return
   endif

 end subroutine COMPUTE_QS_0D_

!#######################################################################

! <SUBROUTINE NAME="compute_qs_1d" INTERFACE="compute_qs">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <IN NAME="press UNIT="Pascals" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="qsat" UNITS="kg(vapor)/kg(moist air)" TYPE="real" DIM="(:)"></OUT>
!   <IN NAME="q" UNIT="kg(vapor)/kg(moistair)" TYPE="real" DIM="(:)"></IN>
!   <IN NAME="hc" UNIT="fraction" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="dqsdT" UNIT="kg(vapor)/kg(moistair)/ degree Kelvin" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="esat" UNITS="Pascals" TYPE="real" DIM="(:)"> </OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine COMPUTE_QS_1D_ ( temp, press, qsat, q, hc, dqsdT, esat, &
                            err_msg, es_over_liq, es_over_liq_and_ice )

 real(kind=FMS_SVP_KIND_), intent(in)   :: temp(:)   !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(in)   :: press(:)  !< pressure
 real(kind=FMS_SVP_KIND_), intent(out)  :: qsat(:)   !< specific humidity
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: q(:)      !< vapor relative humidity
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: hc        !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), optional :: dqsdT(:)  !< d(qsat)/dT
 real(kind=FMS_SVP_KIND_),intent(out), optional :: esat(:) !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg
 logical,intent(in),            optional :: es_over_liq
 logical,intent(in),            optional :: es_over_liq_and_ice

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp
 !> EPSILO and ZVIR are module level variables that are declared in r8_kind.
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(kind=FMS_SVP_KIND_), parameter :: EPSILOl=real(EPSILO,FMS_SVP_KIND_)
 real(kind=FMS_SVP_KIND_), parameter :: ZVIRl=real(ZVIR,FMS_SVP_KIND_)

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   if (present(es_over_liq)) then
     if (.not. (construct_table_wrt_liq)) then
       call error_mesg ('compute_qs', &
          'requesting es wrt liq, but that table not constructed', &
                                                                FATAL)
     endif
   endif
   if (present(es_over_liq_and_ice)) then
     if (.not. (construct_table_wrt_liq_and_ice)) then
       call error_mesg ('compute_qs', &
      'requesting es wrt liq and ice, but that table not constructed', &
                                                                FATAL)
     endif
   endif

!  call compute_qs_k (temp, press,  EPSILO, ZVIR, qsat, nbad, q, dqsdT)
   call compute_qs_k (temp, press,  EPSILOl, ZVIRl, qsat, nbad, q, hc, &
                       dqsdT, esat, es_over_liq, es_over_liq_and_ice)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('compute_qs',err_msg_tmp,err_msg)) return
   endif

 end subroutine COMPUTE_QS_1D_


!#######################################################################

! <SUBROUTINE NAME="compute_qs_2d" INTERFACE="compute_qs">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <IN NAME="press UNIT="Pascals" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="qsat" UNITS="kg(vapor)/kg(moist air)" TYPE="real" DIM="(;,:)"></OUT>
!   <IN NAME="q" UNIT="kg(vapor)/kg(moistair)" TYPE="real" DIM="(:,:)"></IN>
!   <IN NAME="hc" UNIT="fraction" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="dqsdT" UNIT="kg(vapor)/kg(moistair)/ degree Kelvin" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="esat" UNITS="Pascals" TYPE="real" DIM="(:,:)"> </OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine COMPUTE_QS_2D_ ( temp, press, qsat, q, hc, dqsdT, esat, &
                            err_msg, es_over_liq, es_over_liq_and_ice )

 real(kind=FMS_SVP_KIND_), intent(in)   :: temp(:,:)   !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(in)   :: press(:,:)  !< pressure
 real(kind=FMS_SVP_KIND_), intent(out)  :: qsat(:,:)   !< specific humidity
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: q(:,:)     !< vapor relative humidity
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: hc         !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), optional :: dqsdT(:,:) !< d(qsat)/dT
 real(kind=FMS_SVP_KIND_), intent(out), optional :: esat(:,:)  !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg
 logical,intent(in),            optional :: es_over_liq
 logical,intent(in),            optional :: es_over_liq_and_ice

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp
 !> EPSILO and ZVIR are module level variables that are declared in r8_kind.
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(kind=FMS_SVP_KIND_), parameter :: EPSILOl=real(EPSILO,FMS_SVP_KIND_)
 real(kind=FMS_SVP_KIND_), parameter :: ZVIRl=real(ZVIR,FMS_SVP_KIND_)

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   if (present(es_over_liq)) then
     if (.not. (construct_table_wrt_liq)) then
       call error_mesg ('compute_qs', &
          'requesting es wrt liq, but that table not constructed', &
                                                                FATAL)
     endif
   endif
   if (present(es_over_liq_and_ice)) then
     if (.not. (construct_table_wrt_liq_and_ice)) then
       call error_mesg ('compute_qs', &
      'requesting es wrt liq and ice, but that table not constructed', &
                                                                FATAL)
     endif
   endif

!  call compute_qs_k (temp, press,  EPSILO, ZVIR, qsat, nbad, q, dqsdT)
   call compute_qs_k (temp, press,  EPSILOl, ZVIRl, qsat, nbad, q, hc, &
                       dqsdT, esat, es_over_liq, es_over_liq_and_ice)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('compute_qs',err_msg_tmp,err_msg)) return
   endif

 end subroutine COMPUTE_QS_2D_

!#######################################################################

! <SUBROUTINE NAME="compute_qs_3d" INTERFACE="compute_qs">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <IN NAME="press UNIT="Pascals" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="qsat" UNITS="kg(vapor)/kg(moist air)" TYPE="real" DIM="(;,:,:)"></OUT>
!   <IN NAME="q" UNIT="kg(vapor)/kg(moistair)" TYPE="real" DIM="(:,:,:)"></IN>
!   <IN NAME="hc" UNIT="fraction" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="dqsdT" UNIT="kg(vapor)/kg(moistair)/ degree Kelvin" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="esat" UNITS="Pascals" TYPE="real" DIM="(:,:,:)"> </OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine COMPUTE_QS_3D_ ( temp, press, qsat, q, hc, dqsdT, esat, &
                            err_msg, es_over_liq, es_over_liq_and_ice )

 real(kind=FMS_SVP_KIND_), intent(in)   :: temp(:,:,:)   !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(in)   :: press(:,:,:)  !< pressure
 real(kind=FMS_SVP_KIND_), intent(out)  :: qsat(:,:,:)   !< specific humidity
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: q(:,:,:)     !< vapor relative humidity
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: hc           !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), optional :: dqsdT(:,:,:) !< d(qsat)/dT
 real(kind=FMS_SVP_KIND_), intent(out), optional :: esat(:,:,:)  !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg
 logical,intent(in),            optional :: es_over_liq
 logical,intent(in),            optional :: es_over_liq_and_ice

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp
 !> EPSILO and ZVIR are module level variables that are declared in r8_kind.
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(kind=FMS_SVP_KIND_), parameter :: EPSILOl=real(EPSILO,FMS_SVP_KIND_)
 real(kind=FMS_SVP_KIND_), parameter :: ZVIRl=real(ZVIR,FMS_SVP_KIND_)

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   if (present(es_over_liq)) then
     if (.not. (construct_table_wrt_liq)) then
       call error_mesg ('compute_qs', &
          'requesting es wrt liq, but that table not constructed', &
                                                                FATAL)
     endif
   endif
   if (present(es_over_liq_and_ice)) then
     if (.not. (construct_table_wrt_liq_and_ice)) then
       call error_mesg ('compute_qs', &
      'requesting es wrt liq and ice, but that table not constructed', &
                                                                FATAL)
     endif
   endif

!  call compute_qs_k (temp, press,  EPSILO, ZVIR, qsat, nbad, q, dqsdT)
   call compute_qs_k (temp, press,  EPSILOl, ZVIRl, qsat, nbad, q, hc, &
                       dqsdT, esat, es_over_liq, es_over_liq_and_ice)


   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('compute_qs',err_msg_tmp,err_msg)) return
   endif

 end subroutine COMPUTE_QS_3D_

!#######################################################################
!#######################################################################

! <SUBROUTINE NAME="compute_mrs_0d" INTERFACE="compute_mrs">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(SCALAR)"></IN>
!   <IN NAME="press" UNIT="Pascals" TYPE="real" DIM="(SCALAR)"></IN>
!   <OUT NAME="mrsat" UNITS="kg(vapor)/kg(dry air)" TYPE="real" DIM="(SCALAR</OUT>
!   <IN NAME="mr" UNIT="kg(vapor)/kg(dry air)" TYPE="real" DIM="(SCALAR)"></IN>
!   <IN NAME="hc" UNIT="fraction" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="dmrsdT" UNIT="kg(vapor)/kg(dry air)/ degree Kelvin" TYPE="real" DIM="(SCALAR)"></OUT>
!   <OUT NAME="esat" UNITS="Pascals" TYPE="real" DIM="(scalar)"> </OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine COMPUTE_MRS_0D_ ( temp, press, mrsat, mr, hc, dmrsdT, esat, &
                            err_msg, es_over_liq, es_over_liq_and_ice )

 real(kind=FMS_SVP_KIND_), intent(in)   :: temp   !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(in)   :: press  !< pressure
 real(kind=FMS_SVP_KIND_), intent(out)  :: mrsat  !< mixing ratio at relative humidity
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: mr     !< vapor mixing ratio
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: hc     !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), optional :: dmrsdT !< d(mrsat)/dT
 real(kind=FMS_SVP_KIND_), intent(out), optional :: esat   !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg
 logical,intent(in),            optional :: es_over_liq
 logical,intent(in),            optional :: es_over_liq_and_ice

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp
 !> EPSILO and ZVIR are module level variables that are declared in r8_kind.
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(kind=FMS_SVP_KIND_), parameter :: EPSILOl=real(EPSILO,FMS_SVP_KIND_)
 real(kind=FMS_SVP_KIND_), parameter :: ZVIRl=real(ZVIR,FMS_SVP_KIND_)

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   if (present(es_over_liq)) then
     if (.not. (construct_table_wrt_liq)) then
       call error_mesg ('compute_mrs', &
          'requesting es wrt liq, but that table not constructed', &
                                                                FATAL)
     endif
   endif
   if (present(es_over_liq_and_ice)) then
     if (.not. (construct_table_wrt_liq_and_ice)) then
       call error_mesg ('compute_qs', &
      'requesting es wrt liq and ice, but that table not constructed', &
                                                                FATAL)
     endif
   endif

   call compute_mrs_k (temp, press, EPSILOl, ZVIRl, mrsat, nbad, mr,  &
                     hc, dmrsdT, esat, es_over_liq, es_over_liq_and_ice)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('compute_mrs',err_msg_tmp,err_msg)) return
   endif

 end subroutine COMPUTE_MRS_0D_

!#######################################################################
!#######################################################################

! <SUBROUTINE NAME="compute_mrs_1d" INTERFACE="compute_mrs">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <IN NAME="press" UNIT="Pascals" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="mrsat" UNITS="kg(vapor)/kg(dry air)" TYPE="real" DIM="(:)"></OUT>
!   <IN NAME="mr" UNIT="kg(vapor)/kg(dry air)" TYPE="real" DIM="(:)"></IN>
!   <IN NAME="hc" UNIT="fraction" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="dmrsdT" UNIT="kg(vapor)/kg(dry air)/ degree Kelvin" TYPE="real" DIM="(:)"></OUT>
!   <OUT NAME="esat" UNITS="Pascals" TYPE="real" DIM="(:)"> </OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine COMPUTE_MRS_1D_ ( temp, press, mrsat, mr, hc, dmrsdT, esat,&
                            err_msg, es_over_liq, es_over_liq_and_ice )

 real(kind=FMS_SVP_KIND_), intent(in)   :: temp(:)   !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(in)   :: press(:)  !< pressure
 real(kind=FMS_SVP_KIND_), intent(out)  :: mrsat(:)  !< mixing ratio at relative humidity
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: mr(:)     !< vapor mixing ratio
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: hc        !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), optional :: dmrsdT(:) !< d(mrsat)/dT
 real(kind=FMS_SVP_KIND_), intent(out), optional :: esat(:)   !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg
 logical,intent(in),            optional :: es_over_liq
 logical,intent(in),            optional :: es_over_liq_and_ice

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp
 !> EPSILO and ZVIR are module level variables that are declared in r8_kind.
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(kind=FMS_SVP_KIND_), parameter :: EPSILOl=real(EPSILO,FMS_SVP_KIND_)
 real(kind=FMS_SVP_KIND_), parameter :: ZVIRl=real(ZVIR,FMS_SVP_KIND_)

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   if (present(es_over_liq)) then
     if (.not. (construct_table_wrt_liq)) then
       call error_mesg ('compute_mrs', &
          'requesting es wrt liq, but that table not constructed', &
                                                                FATAL)
     endif
   endif
   if (present(es_over_liq_and_ice)) then
     if (.not. (construct_table_wrt_liq_and_ice)) then
       call error_mesg ('compute_qs', &
      'requesting es wrt liq and ice, but that table not constructed', &
                                                                FATAL)
     endif
   endif

!  call compute_mrs_k (temp, press, EPSILO, ZVIR, mrsat,  &
!                                                     nbad, mr, dmrsdT)
   call compute_mrs_k (temp, press, EPSILOl, ZVIRl, mrsat, nbad, mr,  &
                     hc, dmrsdT, esat, es_over_liq, es_over_liq_and_ice)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('compute_mrs',err_msg_tmp,err_msg)) return
   endif

 end subroutine COMPUTE_MRS_1D_

!#######################################################################

! <SUBROUTINE NAME="compute_mrs_2d" INTERFACE="compute_mrs">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <IN NAME="press" UNIT="Pascals" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="mrsat" UNITS="kg(vapor)/kg(dry air)" TYPE="real" DIM="(:,:)"></OUT>
!   <IN NAME="mr" UNIT="kg(vapor)/kg(dry air)" TYPE="real" DIM="(:,:)"></IN>
!   <IN NAME="hc" UNIT="fraction" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="dmrsdT" UNIT="kg(vapor)/kg(dry air)/ degree Kelvin" TYPE="real" DIM="(:,:)"></OUT>
!   <OUT NAME="esat" UNITS="Pascals" TYPE="real" DIM="(:,:)"> </OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine COMPUTE_MRS_2D_ ( temp, press, mrsat, mr, hc, dmrsdT, esat,&
                            err_msg, es_over_liq, es_over_liq_and_ice )

 real(kind=FMS_SVP_KIND_), intent(in)   :: temp(:,:)   !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(in)   :: press(:,:)  !< pressure
 real(kind=FMS_SVP_KIND_), intent(out)  :: mrsat(:,:)  !< mixing ratio at relative humidity
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: mr(:,:)     !< vapor mixing ratio
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: hc          !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), optional :: dmrsdT(:,:) !< d(mrsat)/dT
 real(kind=FMS_SVP_KIND_), intent(out), optional :: esat(:,:)   !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg
 logical,intent(in),            optional :: es_over_liq
 logical,intent(in),            optional :: es_over_liq_and_ice

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp
 !> EPSILO and ZVIR are module level variables that are declared in r8_kind.
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(kind=FMS_SVP_KIND_), parameter :: EPSILOl=real(EPSILO,FMS_SVP_KIND_)
 real(kind=FMS_SVP_KIND_), parameter :: ZVIRl=real(ZVIR,FMS_SVP_KIND_)

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   if (present(es_over_liq)) then
     if (.not. (construct_table_wrt_liq)) then
       call error_mesg ('compute_mrs', &
          'requesting es wrt liq, but that table not constructed', &
                                                                FATAL)
     endif
   endif
   if (present(es_over_liq_and_ice)) then
     if (.not. (construct_table_wrt_liq_and_ice)) then
       call error_mesg ('compute_qs', &
      'requesting es wrt liq and ice, but that table not constructed', &
                                                                FATAL)
     endif
   endif

!  call compute_mrs_k (temp, press, EPSILO, ZVIR, mrsat,  &
!                                                     nbad, mr, dmrsdT)
   call compute_mrs_k (temp, press, EPSILOl, ZVIRl, mrsat, nbad, mr,  &
                     hc, dmrsdT, esat, es_over_liq, es_over_liq_and_ice)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('compute_mrs',err_msg_tmp,err_msg)) return
   endif

 end subroutine COMPUTE_MRS_2D_

!#######################################################################

! <SUBROUTINE NAME="compute_mrs_3d" INTERFACE="compute_mrs">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <IN NAME="press" UNIT="Pascals" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="mrsat" UNITS="kg(vapor)/kg(dry air)" TYPE="real" DIM="(:,:,:)"></OUT>
!   <IN NAME="mr" UNIT="kg(vapor)/kg(dry air)" TYPE="real" DIM="(:,:,:)"></IN>
!   <IN NAME="hc" UNIT="fraction" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="dmrsdT" UNIT="kg(vapor)/kg(dry air)/ degree Kelvin" TYPE="real" DIM="(:,:,:)"></OUT>
!   <OUT NAME="esat" UNITS="Pascals" TYPE="real" DIM="(:,:,:)"> </OUT>
!   <OUT NAME="err_msg" TYPE="character">  </OUT>
! </SUBROUTINE>
 subroutine COMPUTE_MRS_3D_ ( temp, press, mrsat, mr, hc, dmrsdT, esat,&
                            err_msg, es_over_liq, es_over_liq_and_ice )

 real(kind=FMS_SVP_KIND_), intent(in)   :: temp(:,:,:)    !< temperature in degrees Kelvin (K)
 real(kind=FMS_SVP_KIND_), intent(in)   :: press(:,:,:)   !< pressure
 real(kind=FMS_SVP_KIND_), intent(out)  :: mrsat(:,:,:)   !< mixing ratio at relative humidity
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: mr(:,:,:)     !< vapor mixing ratio
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: hc            !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), optional :: dmrsdT(:,:,:) !< d(mrsat)/dT
 real(kind=FMS_SVP_KIND_), intent(out), optional :: esat(:,:,:)   !< saturation vapor pressure
 character(len=*), intent(out), optional :: err_msg
 logical,intent(in),            optional :: es_over_liq
 logical,intent(in),            optional :: es_over_liq_and_ice

 integer :: nbad !< if temperature is out of range
 character(len=128) :: err_msg_tmp
 !> EPSILO and ZVIR are module level variables that are declared in r8_kind.
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(kind=FMS_SVP_KIND_), parameter :: EPSILOl=real(EPSILO,FMS_SVP_KIND_)
 real(kind=FMS_SVP_KIND_), parameter :: ZVIRl=real(ZVIR,FMS_SVP_KIND_)

   if (.not.module_is_initialized) then
      if(fms_error_handler('lookup_es','sat_vapor_pres_init is not called' ,err_msg)) return
   endif

   if (present(es_over_liq)) then
     if (.not. (construct_table_wrt_liq)) then
       call error_mesg ('compute_mrs', &
          'requesting es wrt liq, but that table not constructed', &
                                                                FATAL)
     endif
   endif
   if (present(es_over_liq_and_ice)) then
     if (.not. (construct_table_wrt_liq_and_ice)) then
       call error_mesg ('compute_qs', &
      'requesting es wrt liq and ice, but that table not constructed', &
                                                                FATAL)
     endif
   endif

!  call compute_mrs_k (temp, press, EPSILO, ZVIR, mrsat,   &
!                                                    nbad, mr, dmrsdT)
   call compute_mrs_k (temp, press, EPSILOl, ZVIRl, mrsat, nbad, mr,  &
                     hc, dmrsdT, esat, es_over_liq, es_over_liq_and_ice)

   if ( nbad == 0 ) then
     if(present(err_msg)) err_msg = ''
   else
     if(show_bad_value_count_by_slice) call temp_check ( temp )
     if(show_all_bad_values) call show_all_bad ( temp )
     write(err_msg_tmp,'(a47,i7)') 'saturation vapor pressure table overflow, nbad=', nbad
     if(fms_error_handler('compute_mrs',err_msg_tmp,err_msg)) return
   endif

 end subroutine COMPUTE_MRS_3D_


!#######################################################################

 function CHECK_1D_ ( temp ) result ( nbad )
 real(kind=FMS_SVP_KIND_)   , intent(in)  :: temp(:) !< temperature in degrees Kelvin (K)
 integer :: nbad, ind, i

 !> DTINV, TMIN, TEPS are module level variables declared in r8_kind
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(FMS_SVP_KIND_) :: dtinvll !< local version of module variable dtinvl
 real(FMS_SVP_KIND_) :: tminll  !< local version of module variable tminl
 real(FMS_SVP_KIND_) :: tepsll  !< local version of module variable tepsl

 dtinvll=real(dtinv,FMS_SVP_KIND_)
 tminll=real(tmin,FMS_SVP_KIND_)
 tepsll=real(teps,FMS_SVP_KIND_)

   nbad = 0
   do i = 1, size(temp,1)
     ind = int( dtinvll*(temp(i)-tminll + tepsll) )
     if (ind < 0 .or. ind > nlim) nbad = nbad+1
   enddo

 end function CHECK_1D_

!------------------------------------------------

 function CHECK_2D_ ( temp ) result ( nbad )
 real(kind=FMS_SVP_KIND_)   , intent(in)  :: temp(:,:) !< temperature in degrees Kelvin (K)
 integer :: nbad
 integer :: j

    nbad = 0
    do j = 1, size(temp,2)
      nbad = nbad + check_1d ( temp(:,j) )
    enddo
 end function CHECK_2D_

!#######################################################################

 subroutine TEMP_CHECK_1D_ ( temp )
 real(kind=FMS_SVP_KIND_)   , intent(in) :: temp(:) !< temperature in degrees Kelvin (K)
 integer :: i, unit

   unit = stdoutunit
   write(unit,*) 'Bad temperatures (dimension 1): ', (check_1d(temp(i:i)),i=1,size(temp,1))

 end subroutine TEMP_CHECK_1D_

!--------------------------------------------------------------

 subroutine TEMP_CHECK_2D_ ( temp )
 real(kind=FMS_SVP_KIND_)   , intent(in) :: temp(:,:) !< temperature in degrees Kelvin (K)
 integer :: i, j, unit

   unit = stdoutunit
   write(unit,*) 'Bad temperatures (dimension 1): ', (check_1d(temp(i,:)),i=1,size(temp,1))
   write(unit,*) 'Bad temperatures (dimension 2): ', (check_1d(temp(:,j)),j=1,size(temp,2))

 end subroutine TEMP_CHECK_2D_

!--------------------------------------------------------------

 subroutine TEMP_CHECK_3D_ ( temp )
 real(kind=FMS_SVP_KIND_), intent(in)  :: temp(:,:,:) !< temperature in degrees Kelvin (K)
 integer :: i, j, k, unit

   unit = stdoutunit
   write(unit,*) 'Bad temperatures (dimension 1): ', (check_2d(temp(i,:,:)),i=1,size(temp,1))
   write(unit,*) 'Bad temperatures (dimension 2): ', (check_2d(temp(:,j,:)),j=1,size(temp,2))
   write(unit,*) 'Bad temperatures (dimension 3): ', (check_2d(temp(:,:,k)),k=1,size(temp,3))

 end subroutine TEMP_CHECK_3D_

!#######################################################################

 subroutine SHOW_ALL_BAD_0D_ ( temp )
 real(kind=FMS_SVP_KIND_)   , intent(in) :: temp !< temperature in degrees Kelvin (K)
 integer :: ind, unit
 !> DTINV, TMIN, TEPS are module level variables declared in r8_kind
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(FMS_SVP_KIND_) :: dtinvll !< local version of module variable dtinvl
 real(FMS_SVP_KIND_) :: tminll  !< local version of module variable tminl
 real(FMS_SVP_KIND_) :: tepsll  !< local version of module variable tepsl

 dtinvll=real(dtinv,FMS_SVP_KIND_)
 tminll=real(tmin,FMS_SVP_KIND_)
 tepsll=real(teps,FMS_SVP_KIND_)

 unit = stdoutunit
 ind = int( dtinvll*(temp-tminll+tepsll) )
 if (ind < 0 .or. ind > nlim) then
   write(unit,'(a,e10.3,a,i6)') 'Bad temperature=',temp,' pe=',mpp_pe()
 endif

 end subroutine SHOW_ALL_BAD_0D_

!--------------------------------------------------------------

 subroutine SHOW_ALL_BAD_1D_ ( temp )
 real(kind=FMS_SVP_KIND_)   , intent(in) :: temp(:) !< temperature in degrees Kelvin (K)
 integer :: i, ind, unit
 !> DTINV, TMIN, TEPS are module level variables declared in r8_kind
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(FMS_SVP_KIND_) :: dtinvll !< local version of module variable dtinvl
 real(FMS_SVP_KIND_) :: tminll  !< local version of module variable tminl
 real(FMS_SVP_KIND_) :: tepsll  !< local version of module variable tepsl

 dtinvll=real(dtinv,FMS_SVP_KIND_)
 tminll=real(tmin,FMS_SVP_KIND_)
 tepsll=real(teps,FMS_SVP_KIND_)

 unit = stdoutunit
 do i=1,size(temp)
   ind = int( dtinvll*(temp(i)-tminll+tepsll) )
   if (ind < 0 .or. ind > nlim) then
     write(unit,'(a,e10.3,a,i4,a,i6)') 'Bad temperature=',temp(i),'  at i=',i,' pe=',mpp_pe()
   endif
 enddo

 end subroutine SHOW_ALL_BAD_1D_

!--------------------------------------------------------------

 subroutine SHOW_ALL_BAD_2D_ ( temp )
 real(kind=FMS_SVP_KIND_)   , intent(in) :: temp(:,:) !< temperature in degrees Kelvin (K)
 integer :: i, j, ind, unit
 !> DTINV, TMIN, TEPS are module level variables declared in r8_kind
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(FMS_SVP_KIND_) :: dtinvll !< local version of module variable dtinvl
 real(FMS_SVP_KIND_) :: tminll  !< local version of module variable tminl
 real(FMS_SVP_KIND_) :: tepsll  !< local version of module variable tepsl

 dtinvll=real(dtinv,FMS_SVP_KIND_)
 tminll=real(tmin,FMS_SVP_KIND_)
 tepsll=real(teps,FMS_SVP_KIND_)

 unit = stdoutunit
 do j=1,size(temp,2)
 do i=1,size(temp,1)
   ind = int( dtinvll*(temp(i,j)-tminll+tepsll) )
   if (ind < 0 .or. ind > nlim) then
     write(unit,'(a,e10.3,a,i4,a,i4,a,i6)') 'Bad temperature=',temp(i,j),'  at i=',i,' j=',j,' pe=',mpp_pe()
   endif
 enddo
 enddo

 end subroutine SHOW_ALL_BAD_2D_

!--------------------------------------------------------------

 subroutine SHOW_ALL_BAD_3D_ ( temp )
 real(kind=FMS_SVP_KIND_), intent(in)  :: temp(:,:,:) !< temperature in degrees Kelvin (K)
 integer :: i, j, k, ind, unit
 !> DTINV, TMIN, TEPS are module level variables declared in r8_kind
 !! Thus they need to be converted to FMS_SVP_KIND_
 real(FMS_SVP_KIND_) :: dtinvll !< local version of module variable dtinvl
 real(FMS_SVP_KIND_) :: tminll  !< local version of module variable tminl
 real(FMS_SVP_KIND_) :: tepsll  !< local version of module variable tepsl

 dtinvll=real(dtinv,FMS_SVP_KIND_)
 tminll=real(tmin,FMS_SVP_KIND_)
 tepsll=real(teps,FMS_SVP_KIND_)

 unit = stdoutunit
 do k=1,size(temp,3)
 do j=1,size(temp,2)
 do i=1,size(temp,1)
   ind = int( dtinvll*(temp(i,j,k)-tminll+tepsll) )
   if (ind < 0 .or. ind > nlim) then
     write(unit,'(a,e10.3,a,i4,a,i4,a,i4,a,i6)') 'Bad temperature=',temp(i,j,k),'  at i=',i,' j=',j,' k=',k, &
          & ' pe=',mpp_pe()
   endif
 enddo
 enddo
 enddo

 end subroutine SHOW_ALL_BAD_3D_

!--------------------------------------------------------------

!> @}
