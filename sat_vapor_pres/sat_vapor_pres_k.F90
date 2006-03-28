
!VERSION NUMBER:
!  $Id: sat_vapor_pres_k.F90,v 13.0 2006/03/28 21:42:55 fms Exp $

!module sat_vapor_pres_inter_mod

!#include "sat_vapor_pres_interfaces.h"

!end module sat_vapor_pres_inter_mod

!!RSH NOTE:  
!!       The following subroutine is applicable only to the case where
!!       the input argument temp and output argument esat are scalars.
!!       This is the way sat_vapor_pres_mod is used by donner_deep,
!!       the only module for which this kernelized routine is currently
!!       used. If this is seen to be a feasible way of addressing the
!!       use of sat_vapor_pres, then the additional interfaces currently
!!       provided in sat_vapor_pres_mod need also be converted.
!!       Note also that the functionality of temp_check is NOT included
!!       here; this should also be addressed.


!######################################################################

subroutine sat_vapor_pres_hold_tables_k   &
       (nsize_in, ermesg, table, dtable, d2table, temp, esat)         

use sat_vapor_pres_params_mod, only :  NSIZE, TCMIN, ESRES, NLIM, &
                                       TFREEZE

implicit none

integer,                   intent(in)            :: nsize_in
character(len=*),          intent(out)           :: ermesg
real, dimension(nsize_in), intent(in)            :: table, dtable,  &
                                                    d2table
real,                      intent(in)            :: temp               
real,                      intent(out)           :: esat               


      real, dimension(NSIZE), SAVE  :: table_sv, dtable_sv, d2table_sv
      logical, SAVE ::  table_present = .false.
      real    :: tmp, del, dtinv, teps, dtres, tmin
      integer :: ind


      ermesg = ' '

      if (nsize_in /= NSIZE) then
        ermesg = 'size of saved tables will differ from the input nsize'
        return
      endif

      if (.not. table_present) then
        table_sv = table
        dtable_sv = dtable
        d2table_sv = d2table
        table_present = .true.
      else
        tmin = real(tcmin) + tfreeze
        dtres = 1./real(esres)
        dtinv = real(esres)
        teps = 1./real(2*esres)
        tmp = temp-tmin
        ind = int(dtinv*(tmp+teps))
        del = tmp-dtres*real(ind)
        esat = TABLE_sv(ind+1) + del*(DTABLE_sv(ind+1) +   &
                                      del*D2TABLE_sv(ind+1))

!       if (ind < 0 .or. ind > nlim) call temp_check ( 1, temp )
          
        if (ind < 0 .or. ind > nlim) then                          
          ermesg =     &
                  ' sat_vapor_pres_hold_tables_k: ind < 0 or ind > nlim'
          return
        endif
      endif
          
end subroutine sat_vapor_pres_hold_tables_k


!######################################################################

subroutine sat_vapor_pres_lookup_es_k (temp, esat, ermesg)

use sat_vapor_pres_params_mod, only :  NSIZE

implicit none

real,             intent(in   )  :: temp
real,             intent(  out)  :: esat
character(len=*), intent(  out)  :: ermesg



      real, dimension(NSIZE) :: table_dum, dtable_dum, d2table_dum


      ermesg = ' '
      call sat_vapor_pres_hold_tables_k (NSIZE, ermesg, table_dum, &
                                 dtable_dum, d2table_dum, temp, esat) 


end subroutine sat_vapor_pres_lookup_es_k



