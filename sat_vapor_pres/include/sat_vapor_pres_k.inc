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
!> @addtogroup sat_vapor_pres_k_mod
!> @brief Kernel module to be used by @ref sat_vapor_pres_mod for
!> @{

!! table lookups and calculations

!> This routine has been generalized to return tables for any temperature range and resolution
!! The TABLEs for saturation vapor pressure are computed  with r8_kind precision since
!! these TABLES are module level variables that are decared as r8_kind.
!! The subroutines compute_es_k, compute_es_k, compute_es_liq_k, and compute_es_liq_ice_k
!! seem to be mostly used to compute the TABLE values (and thus all variables within them can be declared
!! as r8_kind).  However, these compute* subroutines have been modified to accept both r4_kind and r8_kind arguments
!! for general usage and the math can be conducted in either r4_kind and r8_kind.
!! In this initialization routine, r8_kind arguments are passed to these compute* subroutines.
!! This routine does not assume the passed in arguments are in r8_precision.
!! Thus all variables used for the computation of the TABLES (e.g. TABLE, DTABLE*, D2TABLE*) are promoted
!! to r8_kind precision.  All local variables used for computation are in r8_kind precision
!! Thus the TABLEs are constructed as accurately as possible and are promoted down to r4_kind when users
!! pass in r4_kind arguments to the LOOKUP* subroutines.
 subroutine SAT_VAPOR_PRES_INIT_K_(table_size, tcmin, tcmax, TFREEZE, HLV, RVGAS, ES0, err_msg, &
                                  use_exact_qs_input, do_simple,  &
                                  construct_table_wrt_liq, &
                                  construct_table_wrt_liq_and_ice, &
                                  teps, tmin, dtinv)

 integer, intent(in) :: table_size
 real(kind=FMS_SVP_KIND_), intent(in) :: tcmin   !< TABLE(1) = sat vapor pressure at temperature tcmin (deg C)
 real(kind=FMS_SVP_KIND_), intent(in) :: tcmax   !< TABLE(table_size) = sat vapor pressure at temperature tcmax (deg C)
 real(kind=FMS_SVP_KIND_), intent(in) :: TFREEZE !< Conversion to Kelvin
 real(kind=FMS_SVP_KIND_), intent(in) :: HLV     !< Latent heat of evaporation [J/kg]
 real(kind=FMS_SVP_KIND_), intent(in) :: RVGAS   !< Gas constant for water vapor
 real(kind=FMS_SVP_KIND_), intent(in) :: ES0     !< Humidity factor [dimensionless]
 logical, intent(in)  :: use_exact_qs_input
 logical, intent(in)  :: do_simple
 logical, intent(in)  :: construct_table_wrt_liq
 logical, intent(in)  :: construct_table_wrt_liq_and_ice
 character(len=*), intent(out) :: err_msg
 real(kind=FMS_SVP_KIND_), intent(out), optional :: teps
 real(kind=FMS_SVP_KIND_), intent(out), optional :: tmin
 real(kind=FMS_SVP_KIND_), intent(out), optional :: dtinv


!> increment used to generate derivative table
!! the following variables are used in the computation
!! of the *TABLES* (which is defined in r8_kind in sat_vapor_pres_mod)
!! Thus these variables are declared with r8_kind
  real(kind=r8_kind), dimension(3) :: tem(3), es(3)
  real(kind=r8_kind) :: hdtinv, tinrc, tfact
  integer :: i

   err_msg = ''

   if (module_is_initialized) return

   if(allocated(TABLE) .or. allocated(DTABLE) .or. allocated(D2TABLE)) then
      err_msg = 'Attempt to allocate sat vapor pressure tables when already allocated'
      return
   else
      allocate(TABLE(table_size), DTABLE(table_size), D2TABLE(table_size))
   endif

   if (construct_table_wrt_liq) then
      if(allocated(TABLE2) .or. allocated(DTABLE2) .or. allocated(D2TABLE2)) then
        err_msg = 'Attempt to allocate sat vapor pressure table2s when already allocated'
        return
      else
        allocate(TABLE2(table_size), DTABLE2(table_size), D2TABLE2(table_size))
      endif
   endif

   if (construct_table_wrt_liq_and_ice) then
      if(allocated(TABLE3) .or. allocated(DTABLE3) .or. allocated(D2TABLE3)) then
        err_msg = 'Attempt to allocate sat vapor pressure table2s when already allocated'
        return
      else
        allocate(TABLE3(table_size), DTABLE3(table_size), D2TABLE3(table_size))
      endif
   endif

      table_siz = table_size
      dtres  = (real(tcmax,r8_kind)-real(tcmin,r8_kind))/real(table_size-1,r8_kind)
      tminl  = real(tcmin,r8_kind)+real(TFREEZE,r8_kind)  ! minimum valid temp in table
      dtinvl = 1.0_r8_kind/dtres
      tepsl  = 0.5_r8_kind*dtres
      tinrc  = 0.1_r8_kind*dtres
      if(present(teps )) teps =real(tepsl,  FMS_SVP_KIND_)
      if(present(tmin )) tmin =real(tminl,  FMS_SVP_KIND_)
      if(present(dtinv)) dtinv=real(dtinvl, FMS_SVP_KIND_)

!> To be able to compute tables for any temperature range and resolution,
!! and at the same time exactly reproduce answers from memphis revision,
!! it is necessary to compute ftact differently than it is in memphis.
      tfact = 5.0_r8_kind*dtinvl

      hdtinv = 0._r8_kind*dtinvl

!> compute es tables from tcmin to tcmax
!> estimate es derivative with small +/- difference

      if (do_simple) then

        !> TABLE = 610.78ES0*exp(-HLV/RGAS[1/tem - 1.TFREEZE])
        !> DTABLE = HLV(TABLE/RVGAS)^2
        do i = 1, table_size
          tem(1) = tminl + dtres*real(i-1,r8_kind)
          TABLE(i) = real(ES0,r8_kind)*610.78_r8_kind* &
               exp( -real(HLV,r8_kind)/real(RVGAS,r8_kind) * (1.0_r8_kind/tem(1) - 1._r8_kind/real(TFREEZE,r8_kind)) )
          DTABLE(i) = real(HLV,r8_kind)*TABLE(i)/real(RVGAS,r8_kind)/tem(1)**2._r8_kind
        enddo

      else

        do i = 1, table_size
          tem(1) = tminl + dtres*real(i-1,r8_kind)
          tem(2) = tem(1)-tinrc
          tem(3) = tem(1)+tinrc
          es = compute_es_k (tem, real(TFREEZE,r8_kind))
          TABLE(i) = es(1)
          DTABLE(i) = (es(3)-es(2))*tfact
        enddo

      endif !if (do_simple)

!> compute one-half second derivative using centered differences
!! differencing des values in the table

      do i = 2, table_size-1
         D2TABLE(i) = 0.25_r8_kind*dtinvl*(DTABLE(i+1)-DTABLE(i-1))
      enddo
!> one-sided derivatives at boundaries

         D2TABLE(1) = 0.50_r8_kind*dtinvl*(DTABLE(2)-DTABLE(1))

         D2TABLE(table_size) = 0.50_r8_kind*dtinvl*(DTABLE(table_size)-DTABLE(table_size-1))

   if (construct_table_wrt_liq) then
!> compute es tables from tcmin to tcmax
!> estimate es derivative with small +/- difference

      do i = 1, table_size
        tem(1) = tminl + dtres*real(i-1,r8_kind)
        tem(2) = tem(1)-tinrc
        tem(3) = tem(1)+tinrc
!>   pass in flag to force all values to be wrt liquid
        es = compute_es_liq_k (tem, real(TFREEZE,r8_kind))
        TABLE2(i) = es(1)
        DTABLE2(i) = (es(3)-es(2))*tfact
      enddo

!> compute one-half second derivative using centered differences
!! differencing des values in the table

     do i = 2, table_size-1
       D2TABLE2(i) = 0.25_r8_kind*dtinvl*(DTABLE2(i+1)-DTABLE2(i-1))
     enddo
!> one-sided derivatives at boundaries

     D2TABLE2(1) = 0.50_r8_kind*dtinvl*(DTABLE2(2)-DTABLE2(1))

     D2TABLE2(table_size) = 0.50_r8_kind*dtinvl*(DTABLE2(table_size)-DTABLE2(table_size-1))
   endif


   if (construct_table_wrt_liq_and_ice) then
!> compute es tables from tcmin to tcmax
!> estimate es derivative with small +/- difference

      do i = 1, table_size
        tem(1) = tminl + dtres*real(i-1,r8_kind)
        tem(2) = tem(1)-tinrc
        tem(3) = tem(1)+tinrc
!>   pass in flag to force all values to be wrt liquid
        es = compute_es_liq_ice_k (tem, real(TFREEZE,r8_kind))
        TABLE3(i) = es(1)
        DTABLE3(i) = (es(3)-es(2))*tfact
      enddo

!> compute one-half second derivative using centered differences
!! differencing des values in the table

     do i = 2, table_size-1
       D2TABLE3(i) = 0.25_r8_kind*dtinvl*(DTABLE3(i+1)-DTABLE3(i-1))
     enddo
!> one-sided derivatives at boundaries

     D2TABLE3(1) = 0.50_r8_kind*dtinvl*(DTABLE3(2)-DTABLE3(1))

     D2TABLE3(table_size) = 0.50_r8_kind*dtinvl*(DTABLE3(table_size)-DTABLE3(table_size-1))
   endif

      use_exact_qs = use_exact_qs_input
      module_is_initialized = .true.

  end subroutine SAT_VAPOR_PRES_INIT_K_

!#######################################################################

 function COMPUTE_ES_K_(tem, TFREEZE) result (es)
 real(kind=FMS_SVP_KIND_), intent(in) :: tem(:)  !< temperature
 real(kind=FMS_SVP_KIND_), intent(in) :: TFREEZE !< conversion to Kelvin
 real(kind=FMS_SVP_KIND_) :: es(size(tem,1))     !< saturation vapor pressure

 real(kind=FMS_SVP_KIND_) :: x
 real(kind=FMS_SVP_KIND_) :: esice
 real(kind=FMS_SVP_KIND_) :: esh2o
 real(kind=FMS_SVP_KIND_) :: TBASW
 real(kind=FMS_SVP_KIND_) :: TBASI
 integer :: i

 integer, parameter :: lkind=FMS_SVP_KIND_ !< local kind parameter

 real(kind=FMS_SVP_KIND_), parameter :: ESBASW = 101324.60_lkind
 real(kind=FMS_SVP_KIND_), parameter :: ESBASI = 610.71_lkind

 !> one and ten are declared for code readability.  For example, one is easier to read
 !! then 1.0_lkind where lkind=FMS_SVP_KIND_ throughout the code
 real(FMS_SVP_KIND_), parameter :: one=1.0_lkind
 real(FMS_SVP_KIND_), parameter :: ten=10.0_lkind

   TBASW = TFREEZE+100.0_lkind !to Kelvin
   TBASI = TFREEZE

   do i = 1, size(tem)

!>  compute es over ice

     !> x = -9.09718(TBASI/tem-1) - 3.56654log(TBASI/tem) + 0.876793(1-tem/TBASI) + log(ESBASI)
     !! the coded equation below is the commented equation above
     if (tem(i) < TBASI) then
         x = -9.09718_lkind*(TBASI/tem(i)-one)   &
             -3.56654_lkind*log10(TBASI/tem(i)) &
             +0.876793_lkind*(one-tem(i)/TBASI) + log10(ESBASI)
         esice =ten**(x)
     else
         esice = 0.0_lkind
     endif

!>  compute es over water greater than -20 c.
!!  values over 100 c may not be valid
!!  see smithsonian meteorological tables page 350.

     !> x = -7.90298(TBASW/tem-1) + 5.02808log(TBASW/tem)
     !!     -1.3816d-07*10^[11.344(1-tem/TBASW)-1]
     !!     +8.1328d-03*10^[-3.49149(TBASW/tem-1)-1] + log(ESBASW)
     !! the coded equation below is the commented equation above
     if (tem(i) > -20.0_lkind+TBASI) then
         x = -7.90298_lkind*(TBASW/tem(i)-one)   &
             +5.02808_lkind*log10(TBASW/tem(i)) &
             -1.3816e-07_lkind*(ten**((one-tem(i)/TBASW)*11.344_lkind)-one) &
             +8.1328e-03_lkind*(ten**((TBASW/tem(i)-one)*(-3.49149_lkind))-one) &
             +log10(ESBASW)
         esh2o = ten**(x)
     else
         esh2o = 0.0_lkind
     endif

!>  derive blended es over ice and supercooled water between -20c and 0c

     !> es = 0.05*[esice*(TBASI-10)+esh2o*(tem-TBASI+20)]
     !! the coded equation below is the commented equation above
     if (tem(i) <= -20.0_lkind+TBASI) then
         es(i) = esice
     else if (tem(i) >= TBASI) then
         es(i) = esh2o
     else
         es(i) = 0.05_lkind*((TBASI-tem(i))*esice + (tem(i)-TBASI+20.0_lkind)*esh2o)
     endif

   enddo

 end function COMPUTE_ES_K_

!#######################################################################

 function COMPUTE_ES_LIQ_K_(tem, TFREEZE) result (es)
 real(kind=FMS_SVP_KIND_), intent(in) :: tem(:)  !< temperature
 real(kind=FMS_SVP_KIND_), intent(in) :: TFREEZE !< conversion to Kelvin
 real(kind=FMS_SVP_KIND_) :: es(size(tem,1))     !< saturation vapor pressure

 real(kind=FMS_SVP_KIND_) :: x
 real(kind=FMS_SVP_KIND_) :: esh2o
 real(kind=FMS_SVP_KIND_) :: TBASW
 integer :: i

 !> local kind variable
 !! one and ten are declared for code readability.  For example, one is easier to read
 !! then 1.0_lkind where lkind=FMS_SVP_KIND_ throughout the code
 integer, parameter :: lkind=FMS_SVP_KIND_
 real(kind=FMS_SVP_KIND_), parameter :: one=1.0_lkind
 real(kind=FMS_SVP_KIND_), parameter :: ten=10.0_lkind
 real(kind=FMS_SVP_KIND_), parameter :: ESBASW = 101324.60_lkind

   TBASW = TFREEZE+100.0_lkind

   do i = 1, size(tem)

!>  compute es over water for all temps.
!!  values over 100 c may not be valid
!!  see smithsonian meteorological tables page 350.

      !> x = -7.90298(TBASW/tem-1) + 5.02808log(TBASW/tem)
      !!     -1.3816d-07*10^[11.344(1-tem/TBASW)-1]
      !!     +8.1328d-03*10^[-3.49149(TBASW/tem-1)-1] + log(ESBASW)
      !! the coded equation below is the commented equation above
         x = -7.90298_lkind*(TBASW/tem(i)-one)  &
             +5.02808_lkind*log10(TBASW/tem(i)) &
             -1.3816e-07_lkind*(ten**((one-tem(i)/TBASW)*11.344_lkind)-one)   &
             +8.1328e-03_lkind*(ten**((TBASW/tem(i)-one)*(-3.49149_lkind))-one)&
             +log10(ESBASW)
         esh2o = ten**(x)

         es(i) = esh2o

   enddo

 end function COMPUTE_ES_LIQ_K_

!#######################################################################

 function COMPUTE_ES_LIQ_ICE_K_(tem, TFREEZE) result (es)
 real(kind=FMS_SVP_KIND_), intent(in) :: tem(:)  !< temperature
 real(kind=FMS_SVP_KIND_), intent(in) :: TFREEZE !< conversion to Kelvin
 real(kind=FMS_SVP_KIND_) :: es(size(tem,1))     !< saturation vapor pressure

 real(kind=FMS_SVP_KIND_) :: x
 real(kind=FMS_SVP_KIND_) :: TBASW
 real(kind=FMS_SVP_KIND_) :: TBASI
 integer :: i

 integer, parameter :: lkind=FMS_SVP_KIND_
 real(kind=FMS_SVP_KIND_), parameter :: ESBASW = 101324.60_lkind
 real(kind=FMS_SVP_KIND_), parameter :: ESBASI = 610.71_lkind
 !> one and ten are declared for code readability.  For example, one is easier to read
 !! then 1.0_lkind where lkind=FMS_SVP_KIND_ throughout the code
 real(kind=FMS_SVP_KIND_), parameter :: one=1.0_lkind
 real(kind=FMS_SVP_KIND_), parameter :: ten=10.0_lkind

   TBASW = TFREEZE+100.0_lkind
   TBASI = TFREEZE

   do i = 1, size(tem)

     if (tem(i) < TBASI) then

!>  compute es over ice
        !> x= -9.09718(TBASI/tem-1) -3.56654log(TBASI/tem) +0.87679(1-tem/TBASI)+log(EBASI)
        !! the coded equation below is the commented equation above
         x = -9.09718_lkind*(TBASI/tem(i)-one)   &
             -3.56654_lkind*log10(TBASI/tem(i)) &
             +0.876793_lkind*(one-tem(i)/TBASI) + log10(ESBASI)
         es(i) =ten**(x)
     else

!>  compute es over water
!!  values over 100 c may not be valid
!!  see smithsonian meteorological tables page 350.
        !> x = -7.90298(TBASW/tem-1) + 5.02808log(TBASW/tem)
        !!     -1.3816d-07*10^[11.344(1-tem/TBASW)-1]
        !!     +8.1328d-03*10^[-3.49149(TBASW/tem-1)-1] + log(ESBASW)
        !! the coded equation below is the commented equation above
          x = -7.90298_lkind*(TBASW/tem(i)-one) &
              +5.02808_lkind*log10(TBASW/tem(i)) &
              -1.3816e-07_lkind*(ten**((one-tem(i)/TBASW)*11.344_lkind)-one)      &
              +8.1328e-03_lkind*(ten**((TBASW/tem(i)-one)*(-3.49149_lkind))-one) &
             +log10(ESBASW)
         es(i) = ten**(x)
     endif
   enddo

 end function COMPUTE_ES_LIQ_ICE_K_

!#######################################################################

 subroutine COMPUTE_QS_K_3D_ (temp, press, eps, zvir, qs, nbad, q, hc, &
                          dqsdT, esat, es_over_liq, es_over_liq_and_ice)

 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:) :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:) :: press !< pressure
 real(kind=FMS_SVP_KIND_), intent(in)                    :: eps   !< EPSILO
 real(kind=FMS_SVP_KIND_), intent(in)                    :: zvir  !< ZVIR
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:) :: qs    !< specific humidity
 integer, intent(out)                                    :: nbad  !< if temperature is out of range
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:), optional :: q     !< vapor relative humidity
 real(kind=FMS_SVP_KIND_), intent(in),                    optional :: hc    !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:), optional :: dqsdT !< d(qs)/dT
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:), optional :: esat  !< saturation vapor pressure
 logical,intent(in),                  optional :: es_over_liq !< use es table wrt liquid only
 logical,intent(in),                  optional :: es_over_liq_and_ice

 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_, either r4_kind or r8_kind

 real(kind=FMS_SVP_KIND_), dimension(size(temp,1), size(temp,2), size(temp,3)) :: esloc
 real(kind=FMS_SVP_KIND_), dimension(size(temp,1), size(temp,2), size(temp,3)) :: desat
 real(kind=FMS_SVP_KIND_), dimension(size(temp,1), size(temp,2), size(temp,3)) :: denom
 integer :: i, j, k
 real(kind=FMS_SVP_KIND_) :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0_lkind
   endif
 if (present(es_over_liq)) then
   if (present (dqsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dqsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dqsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif
   if (nbad == 0) then
     if (present (q) .and. use_exact_qs) then
       qs = (1.0_lkind + zvir*q)*eps*esloc/press
       if (present (dqsdT)) then
         dqsdT = (1.0_lkind + zvir*q)*eps*desat/press
       endif
     else  ! (present(q))
       denom = press - (1.0_lkind - eps)*esloc
       do k=1,size(qs,3)
         do j=1,size(qs,2)
           do i=1,size(qs,1)
             if (denom(i,j,k) > 0.0_lkind) then
               qs(i,j,k) = eps*esloc(i,j,k)/denom(i,j,k)
             else
               qs(i,j,k) = eps
             endif
           end do
         end do
       end do
       if (present (dqsdT)) then
         dqsdT = eps*press*desat/denom**2
       endif
     endif ! (present(q))
   else ! (nbad = 0)
     qs = -999.0_lkind
     if (present (dqsdT)) then
       dqsdT = -999.0_lkind
     endif
     if (present (esat)) then
       esat = -999.0_lkind
     endif
   endif ! (nbad = 0)


 end subroutine COMPUTE_QS_K_3D_

!#######################################################################

 subroutine COMPUTE_QS_K_2D_ (temp, press, eps, zvir, qs, nbad, q, hc, &
                          dqsdT, esat, es_over_liq, es_over_liq_and_ice)

 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:) :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:) :: press !< pressure
 real(kind=FMS_SVP_KIND_), intent(in)                  :: eps   !< EPSILO
 real(kind=FMS_SVP_KIND_), intent(in)                  :: zvir  !< ZVIR
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:) :: qs    !< specific humidity
 integer, intent(out)                                  :: nbad  !< if temperature is out of range
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:), optional :: q     !< vapor specific humidty
 real(kind=FMS_SVP_KIND_), intent(in),                  optional :: hc    !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:), optional :: dqsdT !< d(qs)/dT
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:), optional :: esat  !<saturation vapor pressure
 logical,intent(in), optional :: es_over_liq
 logical,intent(in), optional :: es_over_liq_and_ice

 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

 real(kind=FMS_SVP_KIND_), dimension(size(temp,1), size(temp,2)) :: esloc
 real(kind=FMS_SVP_KIND_), dimension(size(temp,1), size(temp,2)) :: desat
 real(kind=FMS_SVP_KIND_), dimension(size(temp,1), size(temp,2)) :: denom
 integer :: i, j
 real(kind=FMS_SVP_KIND_)  :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0_lkind
   endif

 if (present(es_over_liq)) then
   if (present (dqsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dqsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dqsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif
   if (nbad == 0) then
     if (present (q) .and. use_exact_qs) then
       qs = (1.0_lkind + zvir*q)*eps*esloc/press
       if (present (dqsdT)) then
         dqsdT = (1.0_lkind + zvir*q)*eps*desat/press
       endif
     else  ! (present(q))
       denom = press - (1.0_lkind - eps)*esloc
      do j=1,size(qs,2)
        do i=1,size(qs,1)
          if (denom(i,j) > 0.0_lkind) then
            qs(i,j) = eps*esloc(i,j)/denom(i,j)
          else
            qs(i,j) = eps
          endif
        end do
      end do
      if (present (dqsdT)) then
        dqsdT = eps*press*desat/denom**2
      endif
    endif ! (present(q))
   else ! (nbad = 0)
     qs = -999.0_lkind
     if (present (dqsdT)) then
       dqsdT = -999.0_lkind
     endif
     if (present (esat)) then
       esat = -999.0_lkind
     endif
   endif ! (nbad = 0)


 end subroutine COMPUTE_QS_K_2D_

!#######################################################################

 subroutine COMPUTE_QS_K_1D_ (temp, press, eps, zvir, qs, nbad, q, hc, &
                          dqsdT, esat, es_over_liq, es_over_liq_and_ice)

 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:) :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:) :: press !< pressure
 real(kind=FMS_SVP_KIND_), intent(in)                :: eps   !< EPSILO
 real(kind=FMS_SVP_KIND_), intent(in)                :: zvir  !< ZVIR
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:) :: qs    !< specific humidity
 integer, intent(out)                                :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:), optional :: q     !< vapor specific humidity
 real(kind=FMS_SVP_KIND_), intent(in),                optional :: hc    !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:), optional :: dqsdT !< d(qs)/dt
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:), optional :: esat  !< saturation vapor pressure
 logical,intent(in), optional :: es_over_liq
 logical,intent(in), optional :: es_over_liq_and_ice

 integer, parameter :: lkind=FMS_SVP_KIND_

 real(kind=FMS_SVP_KIND_), dimension(size(temp,1)) :: esloc
 real(kind=FMS_SVP_KIND_), dimension(size(temp,1)) :: desat
 real(kind=FMS_SVP_KIND_), dimension(size(temp,1)) :: denom
 integer :: i
 real(kind=FMS_SVP_KIND_)    :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0_lkind
   endif

 if (present(es_over_liq)) then
   if (present (dqsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dqsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dqsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif
   if (nbad == 0) then
     if (present (q) .and. use_exact_qs) then
       qs = (1.0_lkind + zvir*q)*eps*esloc/press
       if (present (dqsdT)) then
         dqsdT = (1.0_lkind + zvir*q)*eps*desat/press
       endif
     else  ! (present(q))
       denom = press - (1.0_lkind - eps)*esloc
       do i=1,size(qs,1)
         if (denom(i) >  0.0_lkind) then
           qs(i) = eps*esloc(i)/denom(i)
         else
           qs(i) = eps
         endif
       end do
       if (present (dqsdT)) then
         dqsdT = eps*press*desat/denom**2
       endif
     endif ! (present(q))
   else ! (nbad = 0)
     qs = -999.0_lkind
     if (present (dqsdT)) then
       dqsdT = -999.0_lkind
     endif
     if (present (esat)) then
       esat = -999.0_lkind
     endif
   endif ! (nbad = 0)


 end subroutine COMPUTE_QS_K_1D_

!#######################################################################

 subroutine COMPUTE_QS_K_0D_ (temp, press, eps, zvir, qs, nbad, q, hc, &
                          dqsdT, esat, es_over_liq, es_over_liq_and_ice)

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(in)  :: press !< pressure
 real(kind=FMS_SVP_KIND_), intent(in)  :: eps   !< EPSILO
 real(kind=FMS_SVP_KIND_), intent(in)  :: zvir  !< ZVIR
 real(kind=FMS_SVP_KIND_), intent(out) :: qs    !< specific humidity
 integer, intent(out)                  :: nbad  !< if temperature is out of range
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: q     !< vapor specific humidity
 real(kind=FMS_SVP_KIND_), intent(in),  optional :: hc    !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), optional :: dqsdT !< d(qs)/dT
 real(kind=FMS_SVP_KIND_), intent(out), optional :: esat  !< saturation vapor pressure
 logical,intent(in), optional :: es_over_liq
 logical,intent(in), optional :: es_over_liq_and_ice

 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

 real(kind=FMS_SVP_KIND_) :: esloc
 real(kind=FMS_SVP_KIND_) :: desat
 real(kind=FMS_SVP_KIND_) :: denom
 real(kind=FMS_SVP_KIND_) :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0_lkind
   endif

 if (present(es_over_liq)) then
   if (present (dqsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dqsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dqsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif
   if (nbad == 0) then
     if (present (q) .and. use_exact_qs) then
       qs = (1.0_lkind + zvir*q)*eps*esloc/press
       if (present (dqsdT)) then
         dqsdT = (1.0_lkind + zvir*q)*eps*desat/press
       endif
     else  ! (present(q))
       denom = press - (1.0_lkind - eps)*esloc
       if (denom > 0.0_lkind) then
         qs = eps*esloc/denom
       else
         qs = eps
       endif
       if (present (dqsdT)) then
         dqsdT = eps*press*desat/denom**2
       endif
     endif ! (present(q))
   else ! (nbad = 0)
     qs = -999.0_lkind
     if (present (dqsdT)) then
       dqsdT = -999.0_lkind
     endif
     if (present (esat)) then
       esat = -999.0_lkind
     endif
   endif ! (nbad = 0)


 end subroutine COMPUTE_QS_K_0D_

!#######################################################################

 subroutine COMPUTE_MRS_K_3D_ (temp, press, eps, zvir, mrs, nbad,   &
                 mr, hc, dmrsdT, esat,es_over_liq, es_over_liq_and_ice)

 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:) :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:) :: press !< pressure
 real(kind=FMS_SVP_KIND_), intent(in)                    :: eps   !< EPSILO
 real(kind=FMS_SVP_KIND_), intent(in)                    :: zvir  !< ZVIR
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:) :: mrs   !< mixing ratio at relative humidity
 integer, intent(out)                                    :: nbad  !< if temperature is out of range
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:), optional :: mr !< vapor mixing ratio
 real(kind=FMS_SVP_KIND_), intent(in),                    optional :: hc !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:), optional :: dmrsdT !< d(mrs)/dT
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:), optional :: esat   !< saturation vapor pressure
 logical,intent(in), optional :: es_over_liq
 logical,intent(in), optional :: es_over_liq_and_ice

 real(FMS_SVP_KIND_), dimension(size(temp,1), size(temp,2), size(temp,3)) :: esloc
 real(FMS_SVP_KIND_), dimension(size(temp,1), size(temp,2), size(temp,3)) :: desat
 real(FMS_SVP_KIND_), dimension(size(temp,1), size(temp,2), size(temp,3)) :: denom

 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

 integer :: i, j, k
 real(FMS_SVP_KIND_) :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0_lkind
   endif

 if (present (es_over_liq)) then
   if (present (dmrsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dmrsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dmrsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif
   if (nbad == 0) then
     if (present (mr) .and. use_exact_qs) then
       mrs = (eps + mr)*esloc/press
       if (present (dmrsdT)) then
         dmrsdT =  (eps + mr)*desat/press
       endif
     else ! (present (mr))
       denom = press - esloc
       do k=1,size(mrs,3)
         do j=1,size(mrs,2)
           do i=1,size(mrs,1)
             if (denom(i,j,k) > 0.0_lkind) then
               mrs(i,j,k) = eps*esloc(i,j,k)/denom(i,j,k)
             else
               mrs(i,j,k) = eps
             endif
           end do
         end do
       end do
       if (present (dmrsdT)) then
         dmrsdT = eps*press*desat/denom**2
       endif
     endif !(present (mr))
   else
     mrs = -999.0_lkind
     if (present (dmrsdT)) then
       dmrsdT = -999.0_lkind
     endif
     if (present (esat)) then
       esat = -999.0_lkind
     endif
   endif

 end subroutine COMPUTE_MRS_K_3D_

!#######################################################################

 subroutine COMPUTE_MRS_K_2D_ (temp, press, eps, zvir, mrs, nbad,  &
                 mr, hc, dmrsdT, esat,es_over_liq, es_over_liq_and_ice)

 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:) :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:) :: press !< pressure
 real(kind=FMS_SVP_KIND_), intent(in)                  :: eps   !< EPSILO
 real(kind=FMS_SVP_KIND_), intent(in)                  :: zvir  !< ZVIR
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:) :: mrs   !< mixing ratio at relative humidity
 integer, intent(out)                                  :: nbad  !< if temperature is out of range
 real(kind=FMS_SVP_KIND_), intent(in), dimension(:,:), optional  :: mr !< vapor mixing ratio
 real(kind=FMS_SVP_KIND_), intent(in),                 optional  :: hc !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:), optional :: dmrsdT !< d(mrs)/dT
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:), optional :: esat   !< saturation vapor pressure
 logical,intent(in), optional :: es_over_liq
 logical,intent(in), optional :: es_over_liq_and_ice

 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

 real(kind=FMS_SVP_KIND_), dimension(size(temp,1), size(temp,2)) :: esloc
 real(kind=FMS_SVP_KIND_), dimension(size(temp,1), size(temp,2)) :: desat
 real(kind=FMS_SVP_KIND_), dimension(size(temp,1), size(temp,2)) :: denom
 integer :: i, j
 real(kind=FMS_SVP_KIND_)    :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0_lkind
   endif

 if (present (es_over_liq)) then
   if (present (dmrsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dmrsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dmrsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif
   if (nbad == 0) then
     if (present (mr) .and. use_exact_qs) then
       mrs = (eps + mr)*esloc/press
       if (present (dmrsdT)) then
         dmrsdT = (eps + mr)*desat/press
       endif
     else ! (present (mr))
       denom = press - esloc
       do j=1,size(mrs,2)
         do i=1,size(mrs,1)
           if (denom(i,j) > 0.0_lkind) then
             mrs(i,j) = eps*esloc(i,j)/denom(i,j)
           else
             mrs(i,j) = eps
           endif
         end do
       end do
       if (present (dmrsdT)) then
         dmrsdT = eps*press*desat/denom**2
       endif
     endif !(present (mr))
   else
     mrs = -999.0_lkind
     if (present (dmrsdT)) then
       dmrsdT = -999.0_lkind
     endif
     if (present (esat)) then
       esat = -999.0_lkind
     endif
   endif


 end subroutine COMPUTE_MRS_K_2D_

!#######################################################################

 subroutine COMPUTE_MRS_K_1D_ (temp, press, eps, zvir, mrs, nbad,  &
                 mr, hc, dmrsdT, esat,es_over_liq, es_over_liq_and_ice)

 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:)           :: temp   !< temperature
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:)           :: press  !< pressure
 real(kind=FMS_SVP_KIND_), intent(in)                          :: eps    !< EPSILO
 real(kind=FMS_SVP_KIND_), intent(in)                          :: zvir   !< ZVIR
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:)           :: mrs    !< mixing ratio at relative humidity
 integer, intent(out)                                          :: nbad   !< if temperature is out of range
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:), optional :: mr     !< vapor mixing ratio
 real(kind=FMS_SVP_KIND_), intent(in),                optional :: hc     !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:), optional :: dmrsdT !< d(mrs)/dT
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:), optional :: esat   !< saturation vapor pressure
 logical,intent(in),              optional :: es_over_liq
 logical,intent(in),              optional :: es_over_liq_and_ice

 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

 real(kind=FMS_SVP_KIND_), dimension(size(temp,1)) :: esloc
 real(kind=FMS_SVP_KIND_), dimension(size(temp,1)) :: desat
 real(kind=FMS_SVP_KIND_), dimension(size(temp,1)) :: denom
 integer :: i
 real(kind=FMS_SVP_KIND_) :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0_lkind
   endif

 if (present (es_over_liq)) then
   if (present (dmrsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dmrsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dmrsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif
   if (nbad == 0) then
     if (present (mr) .and. use_exact_qs) then
       mrs = (eps + mr)*esloc/press
       if (present (dmrsdT)) then
         dmrsdT =  (eps + mr)*desat/press
       endif
     else ! (present (mr))
       denom = press - esloc
       do i=1,size(mrs,1)
         if (denom(i) > 0.0_lkind) then
           mrs(i) = eps*esloc(i)/denom(i)
         else
           mrs(i) = eps
         endif
       end do
       if (present (dmrsdT)) then
         dmrsdT = eps*press*desat/denom**2
       endif
     endif !(present (mr))
   else
     mrs = -999.0_lkind
     if (present (dmrsdT)) then
       dmrsdT = -999.0_lkind
     endif
     if (present (esat)) then
       esat = -999.0_lkind
     endif
   endif


 end subroutine COMPUTE_MRS_K_1D_

!#######################################################################

 subroutine COMPUTE_MRS_K_0D_ (temp, press, eps, zvir, mrs, nbad,   &
                 mr, hc, dmrsdT, esat,es_over_liq, es_over_liq_and_ice)

 real(kind=FMS_SVP_KIND_), intent(in)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(in)  :: press !< pressure
 real(kind=FMS_SVP_KIND_), intent(in)  :: eps   !< EPSILO
 real(kind=FMS_SVP_KIND_), intent(in)  :: zvir  !< ZVIR
 real(kind=FMS_SVP_KIND_), intent(out) :: mrs   !< mixing ratio at relative humidity
 integer, intent(out)                  :: nbad  !< if temperature is out of range
 real(kind=FMS_SVP_KIND_), intent(in),  optional  :: mr     !< vapor mixing ratio
 real(kind=FMS_SVP_KIND_), intent(in),  optional  :: hc     !< relative humidity
 real(kind=FMS_SVP_KIND_), intent(out), optional  :: dmrsdT !< d(mrs)/dT
 real(kind=FMS_SVP_KIND_), intent(out), optional  :: esat   !<  saturation vapor pressure
 logical,intent(in), optional :: es_over_liq
 logical,intent(in), optional :: es_over_liq_and_ice

 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

 real(kind=FMS_SVP_KIND_) :: esloc
 real(kind=FMS_SVP_KIND_) :: desat
 real(kind=FMS_SVP_KIND_) :: denom
 real(kind=FMS_SVP_KIND_) :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0_lkind
   endif

 if (present (es_over_liq)) then
   if (present (dmrsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dmrsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dmrsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif
   if (nbad == 0) then
     if (present (mr) .and. use_exact_qs) then
       mrs = (eps + mr)*esloc/press
       if (present (dmrsdT)) then
         dmrsdT = (eps + mr)*desat/press
       endif
     else ! (present (mr))
       denom = press - esloc
       if (denom > 0.0_lkind) then
         mrs = eps*esloc/denom
       else
         mrs = eps
       endif
       if (present (dmrsdT)) then
         dmrsdT = eps*press*desat/denom**2
       endif
     endif !(present (mr))
   else
     mrs = -999.0_lkind
     if (present (dmrsdT)) then
       dmrsdT = -999.0_lkind
     endif
     if (present (esat)) then
       esat = -999.0_lkind
     endif
   endif


 end subroutine COMPUTE_MRS_K_0D_


!#######################################################################

 subroutine LOOKUP_ES_DES_K_3D_ (temp, esat, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:)  :: esat  !<saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:)  :: desat !< derivative of saturation vapor pressure
 integer, intent(out)                                     :: nbad  !< if temperature is out of range

 real(kind=FMS_SVP_KIND_) :: tmp !< temp-TMIN
 real(kind=FMS_SVP_KIND_) :: del !< delta T
 integer :: ind, i, j, k
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz) then
       nbad = nbad+1
     else
       del = tmp-real(dtresl,FMS_SVP_KIND_)*real(ind,FMS_SVP_KIND_)
       !> esat = TABLE + del*(TABLE + del*D2TABLE)
       !! the coded equation below is the commented equation above
       esat(i,j,k) = real(TABLE(ind+1),FMS_SVP_KIND_)              &
                         + del*( real(DTABLE(ind+1),FMS_SVP_KIND_) &
                                    + del*real(D2TABLE(ind+1),FMS_SVP_KIND_) )
       !> desat = DTABLE + 2del*D2TABLE
       !! the coded equation below is the commented equation above
       desat(i,j,k) = real(DTABLE(ind+1), FMS_SVP_KIND_) &
                        + 2.0_lkind*del*real(D2TABLE(ind+1),FMS_SVP_KIND_)
     endif
   enddo
   enddo
   enddo

 end subroutine LOOKUP_ES_DES_K_3D_

!#######################################################################

 subroutine LOOKUP_ES_DES_K_2D_ (temp, esat, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:) :: temp  !<  temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:) :: desat !< derivative of the saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)


   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
       !> esat = TABLE + del*(TABLE + del*D2TABLE)
       !! the coded equation below is the commented equation above
       esat(i,j) = real(TABLE(ind+1),FMS_SVP_KIND_) &
                 + del*( real(DTABLE(ind+1),FMS_SVP_KIND_) &
                 + del*real(D2TABLE(ind+1),FMS_SVP_KIND_) )
       !> desat = DTABLE + 2del*D2TABLE
       !! the coded equation below is the commented equation above
       desat(i,j) = real(DTABLE(ind+1),FMS_SVP_KIND_) &
                  + 2.0_lkind*del*real(D2TABLE(ind+1),FMS_SVP_KIND_)
     endif
   enddo
   enddo

 end subroutine LOOKUP_ES_DES_K_2D_

!#######################################################################

 subroutine LOOKUP_ES_DES_K_1D_ (temp, esat, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:) :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:) :: desat !< derivative of the saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range

 real(kind=FMS_SVP_KIND_) :: tmp !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del !< delta T
 integer :: ind, i
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)


   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> esat = TABLE + del*(TABLE + del*D2TABLE)
        !! the coded equation below is the commented equation above
        esat(i) = real(TABLE(ind+1),FMS_SVP_KIND_) &
               + del*( real(DTABLE(ind+1),FMS_SVP_KIND_) &
                                + del*real(D2TABLE(ind+1),FMS_SVP_KIND_) )
        !> desat = DTABLE + 2del*D2TABLE
        !! the coded equation below is the commented equation above
        desat(i) = real(DTABLE(ind+1),FMS_SVP_KIND_) &
             + 2.0_lkind*del*real(D2TABLE(ind+1),FMS_SVP_KIND_)
     endif
   enddo

 end subroutine LOOKUP_ES_DES_K_1D_

!#######################################################################

 subroutine LOOKUP_ES_DES_K_0D_ (temp, esat, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out) :: desat !< derivative of the saturation vapor pressure
 integer, intent(out) :: nbad !< if temperature is out of range

 real(kind=FMS_SVP_KIND_) :: tmp !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del !< delta T
 integer :: ind
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   tmp = temp-tminll
   ind = int(dtinvll*(tmp+tepsll))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
     del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
     !> esat = TABLE + del*(TABLE + del*D2TABLE)
     !! the coded equation below is the commented equation above
     esat = real(TABLE(ind+1),FMS_SVP_KIND_) &
            + del*( real(DTABLE(ind+1),FMS_SVP_KIND_) &
                          + del*real(D2TABLE(ind+1),FMS_SVP_KIND_) )
     !> desat = DTABLE + 2del*D2TABLE
     !! the coded equation below is the commented equation above
     desat = real(DTABLE(ind+1),FMS_SVP_KIND_) &
             + 2.0_lkind*del*real(D2TABLE(ind+1),FMS_SVP_KIND_)
   endif

 end subroutine LOOKUP_ES_DES_K_0D_

!#######################################################################

 subroutine LOOKUP_ES_K_3D_(temp, esat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:)  :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:)  :: esat !< saturavation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp - TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j, k
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
       !> esat = TABLE + del*(TABLE + del*D2TABLE)
       !! the coded equation below is the commented equation above
       esat(i,j,k) = real(TABLE(ind+1),FMS_SVP_KIND_) &
                     + del*( real(DTABLE(ind+1),FMS_SVP_KIND_) &
                                    + del*real(D2TABLE(ind+1),FMS_SVP_KIND_) )
     endif
   enddo
   enddo
   enddo

 end subroutine LOOKUP_ES_K_3D_

!#######################################################################

 subroutine LOOKUP_DES_K_3D_(temp, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:)  :: desat !< derivative of sat vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j, k
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> desat = DTABLE + 2del*D2TABLE
        !! the coded equation below is the commented equation above
       desat(i,j,k) = real(DTABLE(ind+1),FMS_SVP_KIND_) &
                      + 2.0_lkind*del*real(D2TABLE(ind+1),FMS_SVP_KIND_)
     endif
   enddo
   enddo
   enddo

 end subroutine LOOKUP_DES_K_3D_

!#######################################################################
 subroutine LOOKUP_DES_K_2D_(temp, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:)  :: desat !< derivative of sat vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j

 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> desat = DTABLE + 2del*D2TABLE
        !! the coded equation below is the commented equation above
        desat(i,j) = real(DTABLE(ind+1),FMS_SVP_KIND_) &
             + 2.0_lkind*del*real(D2TABLE(ind+1),FMS_SVP_KIND_)
     endif
   enddo
   enddo

 end subroutine LOOKUP_DES_K_2D_
!#######################################################################
 subroutine LOOKUP_ES_K_2D_(temp, esat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:)  :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:)  :: esat !< saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)


   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
       !> esat = TABLE + del*(TABLE + del*D2TABLE)
       !! the coded equation below is the commented equation above
       esat(i,j) = real(TABLE(ind+1),FMS_SVP_KIND_) &
                   + del*(real(DTABLE(ind+1),FMS_SVP_KIND_) &
                                  + del*real(D2TABLE(ind+1),FMS_SVP_KIND_)  )
     endif
   enddo
   enddo

 end subroutine LOOKUP_ES_K_2D_
!#######################################################################
 subroutine LOOKUP_DES_K_1D_(temp, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:)  :: desat !< derivative of sat vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> desat = DTABLE + 2del*D2TABLE
        !! the coded equation below is the commented equation above
        desat(i) = real(DTABLE(ind+1),FMS_SVP_KIND_) &
                  + 2.0_lkind*del*real(D2TABLE(ind+1),FMS_SVP_KIND_)
     endif
   enddo

 end subroutine LOOKUP_DES_K_1D_
!#######################################################################
 subroutine LOOKUP_ES_K_1D_(temp, esat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:)  :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:)  :: esat !< saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
       !> esat = TABLE + del*(TABLE + del*D2TABLE)
       !! the coded equation below is the commented equation above
       esat(i) = real(TABLE(ind+1),FMS_SVP_KIND_) &
                 + del*(real(DTABLE(ind+1),FMS_SVP_KIND_) &
                                + del*real(D2TABLE(ind+1),FMS_SVP_KIND_) )
     endif
   enddo

 end subroutine LOOKUP_ES_K_1D_
!#######################################################################
 subroutine LOOKUP_DES_K_0D_(temp, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out) :: desat !< derivative of sat vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp - TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   tmp = temp-tminll
   ind = int(dtinvll*(tmp+tepsll))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
      del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
      !> desat = DTABLE + 2del*D2TABLE
      !! the coded equation below is the commented equation above
      desat = real(DTABLE(ind+1),FMS_SVP_KIND_) &
           + 2.0_lkind*del*real(D2TABLE(ind+1),FMS_SVP_KIND_)
   endif

 end subroutine LOOKUP_DES_K_0D_
!#######################################################################
 subroutine LOOKUP_ES_K_0D_(temp, esat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in)     :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out)    :: esat !< saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp - TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !<  delta T
 integer :: ind
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   tmp = temp-tminll
   ind = int(dtinvll*(tmp+tepsll))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
      del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
      !> esat = TABLE + del*(TABLE + del*D2TABLE)
      !! the coded equation below is the commented equation above
      esat = real(TABLE(ind+1),FMS_SVP_KIND_) &
            + del*( real(DTABLE(ind+1),FMS_SVP_KIND_) &
                             + del*real(D2TABLE(ind+1),FMS_SVP_KIND_) )
   endif

 end subroutine LOOKUP_ES_K_0D_
!#######################################################################

 subroutine LOOKUP_ES2_DES2_K_3D_ (temp, esat, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:) :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:) :: desat !< derivative of sat vapor pressure
 integer, intent(out)                                    :: nbad  !< if temperature is out of range

 real(kind=FMS_SVP_KIND_) :: tmp !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del !< delta T
 integer :: ind, i, j, k
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz) then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> esat = TABLE + del*(TABLE + del*D2TABLE)
        !! the coded equation below is the commented equation above
        esat(i,j,k) = real(TABLE2(ind+1),FMS_SVP_KIND_) &
                     + del*( real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                                + del*real(D2TABLE2(ind+1),FMS_SVP_KIND_) )
        !> desat = DTABLE + 2del*D2TABLE
        !! the coded equation below is the commented equation above
        desat(i,j,k) = real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                      + 2.0_lkind*del*real(D2TABLE2(ind+1),FMS_SVP_KIND_)
     endif
   enddo
   enddo
   enddo

 end subroutine LOOKUP_ES2_DES2_K_3D_

!#######################################################################

 subroutine LOOKUP_ES2_DES2_K_2D_ (temp, esat, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:)  :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:)  :: desat !< derivative of sat vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range

 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> esat = TABLE + del*(TABLE + del*D2TABLE)
        !! the coded equation below is the commented equation above
        esat(i,j) = real(TABLE2(ind+1),FMS_SVP_KIND_) &
                   + del*( real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                                 + del*real(D2TABLE2(ind+1),FMS_SVP_KIND_) )
        !> desat = DTABLE + 2del*D2TABLE
        !! the coded equation below is the commented equation above
        desat(i,j) = real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                   + 2.0_lkind*del*real(D2TABLE2(ind+1),FMS_SVP_KIND_)
     endif
   enddo
   enddo

 end subroutine LOOKUP_ES2_DES2_K_2D_

!#######################################################################

 subroutine LOOKUP_ES2_DES2_K_1D_ (temp, esat, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:)  :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:)  :: desat !< derivative of sat vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range

 real(kind=FMS_SVP_KIND_) :: tmp  !< temp - TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> esat = TABLE + del*(TABLE + del*D2TABLE)
        !! the coded equation below is the commented equation above
        esat(i) = real(TABLE2(ind+1),FMS_SVP_KIND_) &
                 + del*( real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                 + del*real(D2TABLE2(ind+1),FMS_SVP_KIND_) )
        !> desat = DTABLE + 2del*D2TABLE
        desat(i) = real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                 + 2.0_lkind*del*real(D2TABLE2(ind+1),FMS_SVP_KIND_)
     endif
   enddo

 end subroutine LOOKUP_ES2_DES2_K_1D_

!#######################################################################

 subroutine LOOKUP_ES2_DES2_K_0D_ (temp, esat, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in)     :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out)    :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out)    :: desat !< derivative of sat vapor pressure
 integer, intent(out) :: nbad    !< if temperature is out of range

 real(kind=FMS_SVP_KIND_) :: tmp !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del !< delta T
 integer :: ind
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)


   nbad = 0
   tmp = temp-tminll
   ind = int(dtinvll*(tmp+tepsll))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
      del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
      !> esat = TABLE + del*(TABLE + del*D2TABLE)
      !! the coded equation below is the commented equation above
      esat = real(TABLE2(ind+1),FMS_SVP_KIND_) &
            + del*( real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                          + del*real(D2TABLE2(ind+1),FMS_SVP_KIND_) )
      !> desat = DTABLE + 2del*D2TABLE
      !! the coded equation below is the commented equation above
      desat = real(DTABLE2(ind+1),FMS_SVP_KIND_) &
             + 2.0_lkind*del*real(D2TABLE2(ind+1),FMS_SVP_KIND_)
   endif

 end subroutine LOOKUP_ES2_DES2_K_0D_

!#######################################################################

 subroutine LOOKUP_ES2_K_3D_(temp, esat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:)  :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:)  :: esat !< saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j, k
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
       !> esat = TABLE + del*(TABLE + del*D2TABLE)
       !! the coded equation below is the commented equation above
       esat(i,j,k) = real(TABLE2(ind+1),FMS_SVP_KIND_) &
                     + del*( real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                                  + del*real(D2TABLE2(ind+1),FMS_SVP_KIND_) )
     endif
   enddo
   enddo
   enddo

 end subroutine LOOKUP_ES2_K_3D_

!#######################################################################

 subroutine LOOKUP_DES2_K_3D_(temp, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:)  :: desat !< derivative of sat vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j, k
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
       !> desat = DTABLE + 2del*D2TABLE
       !! the coded equation below is the commented equation above
       desat(i,j,k) = real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                      + 2.0_lkind*del*real(D2TABLE2(ind+1),FMS_SVP_KIND_)
     endif
   enddo
   enddo
   enddo

 end subroutine LOOKUP_DES2_K_3D_

!#######################################################################
 subroutine LOOKUP_DES2_K_2D_(temp, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:)  :: desat !< derivative of sat vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
       !> desat = DTABLE + 2del*D2TABLE
       !! the coded equation below is the commented equation above
       desat(i,j) = real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                    + 2.0_lkind*del*real(D2TABLE2(ind+1),FMS_SVP_KIND_)
     endif
   enddo
   enddo

 end subroutine LOOKUP_DES2_K_2D_
!#######################################################################
 subroutine LOOKUP_ES2_K_2D_(temp, esat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:)  :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:)  :: esat !< saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,kind=FMS_SVP_KIND_)
        !> esat = TABLE + del*(TABLE + del*D2TABLE)
        !! the coded equation below is the commented equation above
        esat(i,j) = real(TABLE2(ind+1),FMS_SVP_KIND_) &
                   + del*( real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                                   + del*real(D2TABLE2(ind+1),FMS_SVP_KIND_) )
     endif
   enddo
   enddo

 end subroutine LOOKUP_ES2_K_2D_
!#######################################################################
 subroutine LOOKUP_DES2_K_1D_(temp, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:)  :: desat !< derivative of sat vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> desat = DTABLE + 2del*D2TABLE
        !! the coded equation below is the commented equation above
        desat(i) = real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                  + 2.0_lkind*del*real(D2TABLE2(ind+1),FMS_SVP_KIND_)
     endif
   enddo

 end subroutine LOOKUP_DES2_K_1D_
!#######################################################################
 subroutine LOOKUP_ES2_K_1D_(temp, esat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:)  :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:)  :: esat !< saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> esat = TABLE + del*(TABLE + del*D2TABLE)
        !! the coded equation below is the commented equation above
        esat(i) = real(TABLE2(ind+1),FMS_SVP_KIND_) &
                 + del*( real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                               + del*real(D2TABLE2(ind+1),FMS_SVP_KIND_) )
     endif
   enddo

 end subroutine LOOKUP_ES2_K_1D_
!#######################################################################
 subroutine LOOKUP_DES2_K_0D_(temp, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in)     :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out)    :: desat !< derivative of sat vapor pressure
 integer, intent(out)        :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_)    :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_)    :: del  !< delta T
 integer :: ind
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   tmp = temp-tminll
   ind = int(dtinvll*(tmp+tepsll))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
      del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
      !> desat = DTABLE + 2del*D2TABLE
      !! the coded equation below is the commented equation above
      desat = real(DTABLE2(ind+1),FMS_SVP_KIND_) &
            + 2.0_lkind*del*real(D2TABLE2(ind+1),FMS_SVP_KIND_)
   endif

 end subroutine LOOKUP_DES2_K_0D_
!#######################################################################
 subroutine LOOKUP_ES2_K_0D_(temp, esat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in)     :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out)    :: esat !< saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   tmp = temp-tminll
   ind = int(dtinvll*(tmp+tepsll))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
      del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
      !> esat = TABLE + del*(TABLE + del*D2TABLE)
      !! the coded equation below is the commented equation above
      esat = real(TABLE2(ind+1),FMS_SVP_KIND_) &
            + del*( real(DTABLE2(ind+1),FMS_SVP_KIND_) &
                       + del*real(D2TABLE2(ind+1),FMS_SVP_KIND_) )
   endif

 end subroutine LOOKUP_ES2_K_0D_
!#######################################################################

!#######################################################################

 subroutine LOOKUP_ES3_DES3_K_3D_ (temp, esat, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:)  :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:)  :: desat !< derivative of esat
 integer, intent(out)     :: nbad !< if temperature is out of range

 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !<  delta T
 integer :: ind, i, j, k
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz) then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> esat = TABLE + del*(TABLE + del*D2TABLE)
        !! the coded equation below is the commented equation above
        esat(i,j,k) = real(TABLE3(ind+1),FMS_SVP_KIND_) &
                     + del*( real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                                  + del*real(D2TABLE3(ind+1),FMS_SVP_KIND_) )
        !> desat = DTABLE + 2del*D2TABLE
        !! the coded equation below is the commented equation above
        desat(i,j,k) = real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                      + 2.0_lkind*del*real(D2TABLE3(ind+1),FMS_SVP_KIND_)
     endif
   enddo
   enddo
   enddo

 end subroutine LOOKUP_ES3_DES3_K_3D_

!#######################################################################

 subroutine LOOKUP_ES3_DES3_K_2D_ (temp, esat, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:)  :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:)  :: desat !< derivative of desat
 integer, intent(out)     :: nbad !< if temperature is out of range

 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
       !> esat = TABLE + del*(TABLE + del*D2TABLE)
       !! the coded equation below is the commented equation above
       esat(i,j) = real(TABLE3(ind+1),FMS_SVP_KIND_) &
                   + del*( real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                   + del*real(D2TABLE3(ind+1),FMS_SVP_KIND_) )
       !> desat = DTABLE + 2del*D2TABLE
       !! the coded equation below is the commented equation above
       desat(i,j) = real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                    + 2.0_lkind*del*real(D2TABLE3(ind+1),FMS_SVP_KIND_)
     endif
   enddo
   enddo

 end subroutine LOOKUP_ES3_DES3_K_2D_

!#######################################################################

 subroutine LOOKUP_ES3_DES3_K_1D_ (temp, esat, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:)  :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:)  :: desat !< derivative of esat
 integer, intent(out)     :: nbad !< if temperature is out of range

 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> esat = TABLE + del*(TABLE + del*D2TABLE)
        !! the coded equation below is the commented equation above
        esat(i) = real(TABLE3(ind+1),FMS_SVP_KIND_) &
                 + del*( real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                 + del*real(D2TABLE3(ind+1),FMS_SVP_KIND_) )
        !> desat = DTABLE + 2del*D2TABLE
        !! the coded equation below is the commented equation above
        desat(i) = real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                  + 2.0_lkind*del*real(D2TABLE3(ind+1),FMS_SVP_KIND_)
     endif
   enddo

 end subroutine LOOKUP_ES3_DES3_K_1D_

!#######################################################################

 subroutine LOOKUP_ES3_DES3_K_0D_ (temp, esat, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out) :: esat  !< saturation vapor pressure
 real(kind=FMS_SVP_KIND_), intent(out) :: desat !< derivative of esat
 integer, intent(out)     :: nbad !< if temperature is out of range

 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   tmp = temp-tminll
   ind = int(dtinvll*(tmp+tepsll))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
      del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
      !> esat = TABLE + del*(TABLE + del*D2TABLE)
      !! the coded equation below is the commented equation above
      esat = real(TABLE3(ind+1),FMS_SVP_KIND_) &
            + del*( real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                           + del*real(D2TABLE3(ind+1),FMS_SVP_KIND_) )
      !> desat = DTABLE + 2del*D2TABLE
      !! the coded equation below is the commented equation above
      desat = real(DTABLE3(ind+1),FMS_SVP_KIND_) &
             + 2.0_lkind*del*real(D2TABLE3(ind+1),FMS_SVP_KIND_)
   endif

 end subroutine LOOKUP_ES3_DES3_K_0D_

!#######################################################################

 subroutine LOOKUP_ES3_K_3D_(temp, esat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:)  :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:)  :: esat !< saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j, k
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
       !> esat = TABLE + del*(TABLE + del*D2TABLE)
       !! the coded equation below is the commented equation above
       esat(i,j,k) = real(TABLE3(ind+1),FMS_SVP_KIND_) &
                     + del*( real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                                    + del*real(D2TABLE3(ind+1),FMS_SVP_KIND_) )
     endif
   enddo
   enddo
   enddo

 end subroutine LOOKUP_ES3_K_3D_

!#######################################################################

 subroutine LOOKUP_DES3_K_3D_(temp, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:,:)  :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:,:)  :: desat!< derivatove of saturation vap pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j, k
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> desat = DTABLE + 2del*D2TABLE
        !! the coded equation below is the commented equation above
        desat(i,j,k) = real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                      + 2.0_lkind*del*real(D2TABLE3(ind+1),FMS_SVP_KIND_)
     endif
   enddo
   enddo
   enddo

 end subroutine LOOKUP_DES3_K_3D_

!#######################################################################
 subroutine LOOKUP_DES3_K_2D_(temp, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:)  :: desat !< derivative of sat vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> desat = DTABLE + 2del*D2TABLE
        !! the coded equation below is the commented equation above
        desat(i,j) = real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                    + 2.0_lkind*del*real(D2TABLE3(ind+1),FMS_SVP_KIND_)
     endif
   enddo
   enddo

 end subroutine LOOKUP_DES3_K_2D_
!#######################################################################
 subroutine LOOKUP_ES3_K_2D_(temp, esat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:,:)  :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:,:)  :: esat !< saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i, j
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> esat = TABLE + del*(TABLE + del*D2TABLE)
        !! the coded equation below is the commented equation above
        esat(i,j) = real(TABLE3(ind+1),FMS_SVP_KIND_) &
                   + del*( real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                                 + del*real(D2TABLE3(ind+1),FMS_SVP_KIND_) )
     endif
   enddo
   enddo

 end subroutine LOOKUP_ES3_K_2D_
!#######################################################################
 subroutine LOOKUP_DES3_K_1D_(temp, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:)  :: desat !< derivative of sat vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> desat = DTABLE + 2del*D2TABLE
        !! the coded equation below is the commented equation above
        desat(i) = real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                  + 2.0_lkind*del*real(D2TABLE3(ind+1),FMS_SVP_KIND_)
     endif
   enddo

 end subroutine LOOKUP_DES3_K_1D_
!#######################################################################
 subroutine LOOKUP_ES3_K_1D_(temp, esat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in),  dimension(:)  :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out), dimension(:)  :: esat !< saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind, i
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminll
     ind = int(dtinvll*(tmp+tepsll))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
        del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
        !> esat = TABLE + del*(TABLE + del*D2TABLE)
        !! the coded equation below is the commented equation above
        esat(i) = real(TABLE3(ind+1),FMS_SVP_KIND_) &
                 + del*( real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                                + del*real(D2TABLE3(ind+1),FMS_SVP_KIND_) )
     endif
   enddo

 end subroutine LOOKUP_ES3_K_1D_
!#######################################################################
 subroutine LOOKUP_DES3_K_0D_(temp, desat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in)  :: temp  !< temperature
 real(kind=FMS_SVP_KIND_), intent(out) :: desat !< derivative of sat vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll
 integer, parameter :: lkind=FMS_SVP_KIND_ !< local FMS_SVP_KIND_

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   tmp = temp-tminll
   ind = int(dtinvll*(tmp+tepsll))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
      del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
      !> desat = DTABLE + 2del*D2TABLE
      !! the coded equation below is the commented equation above
      desat = real(DTABLE3(ind+1),FMS_SVP_KIND_) &
             + 2.0_lkind*del*real(D2TABLE3(ind+1),FMS_SVP_KIND_)
   endif

 end subroutine LOOKUP_DES3_K_0D_
!#######################################################################
 subroutine LOOKUP_ES3_K_0D_(temp, esat, nbad)
 real(kind=FMS_SVP_KIND_), intent(in)     :: temp !< temperature
 real(kind=FMS_SVP_KIND_), intent(out)    :: esat !< saturation vapor pressure
 integer, intent(out)     :: nbad !< if temperature is out of range
 real(kind=FMS_SVP_KIND_) :: tmp  !< temp-TMINLL
 real(kind=FMS_SVP_KIND_) :: del  !< delta T
 integer :: ind
 !> dtres, tminl, tepsl, dtinvl are module level variables declared in r8_kind precision
 !! for precision consistency and for code readability, the *ll variables are declared
 !! and used
 real(kind=FMS_SVP_KIND_) :: dtresl
 real(kind=FMS_SVP_KIND_) :: tepsll
 real(kind=FMS_SVP_KIND_) :: tminll
 real(kind=FMS_SVP_KIND_) :: dtinvll

   dtresl=real(dtres, FMS_SVP_KIND_)
   tminll=real(tminl, FMS_SVP_KIND_)
   tepsll=real(tepsl, FMS_SVP_KIND_)
   dtinvll=real(dtinvl, FMS_SVP_KIND_)

   nbad = 0
   tmp = temp-tminll
   ind = int(dtinvll*(tmp+tepsll))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
      del = tmp-dtresl*real(ind,FMS_SVP_KIND_)
      !> esat = TABLE + del*(TABLE + del*D2TABLE)
      !! the coded equation below is the commented equation above
      esat = real(TABLE3(ind+1),FMS_SVP_KIND_) &
            + del*( real(DTABLE3(ind+1),FMS_SVP_KIND_) &
                          + del*real(D2TABLE3(ind+1),FMS_SVP_KIND_) )
   endif

 end subroutine LOOKUP_ES3_K_0D_
!#######################################################################
!> @}
! close documentation grouping
