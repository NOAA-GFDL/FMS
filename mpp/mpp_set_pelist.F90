module mpp_set_pelist_mod


use mpp_mod, only : mpp_declare_pelist, mpp_npes, mpp_pe

implicit none
private


! public data and interface
public :: pelist1, pelist2, atm_pelist, ocn_pelist, group1, group2, &
          mpp_set_pelist_init, mpp_set_pelist_end, check_parallel,  &
          atm1_pelist, atm2_pelist, ocn1_pelist, ocn2_pelist,       &
          atm1_pe, ocn1_pe, atm2_pe, ocn2_pe, atm_pe, ocn_pe

logical :: group1=.false., group2=.false., check_parallel=.false.
logical :: atm1_pe = .false., atm2_pe = .false., atm_pe = .false.
logical :: ocn1_pe = .false., ocn2_pe = .false., ocn_pe = .false.

integer, dimension(:), allocatable :: pelist1, pelist2, atm_pelist, ocn_pelist, &
                                      atm1_pelist, atm2_pelist, ocn1_pelist, ocn2_pelist

!-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: mpp_set_pelist.F90,v 6.1 2003/04/09 21:18:11 fms Exp $'
  character(len=128) :: tagname = '$Name: inchon $'
!-----------------------------------------------------------------------

contains

!##################################################################
! set up the pelist for parallel checking
   subroutine mpp_set_pelist_init(npes1, npes2, atm1_pe_start, atm1_pe_end,  &
                                  ocn1_pe_start, ocn1_pe_end, atm2_pe_start, &
                                  atm2_pe_end, ocn2_pe_start, ocn2_pe_end )

   integer, intent(in) :: npes1, npes2 ! number of pes of the two 
                                       ! concurrent running ensembles. 
                                       ! when check is false, npes2 should be 0
   integer, intent(in) :: atm1_pe_start, atm1_pe_end, &
                          ocn1_pe_start, ocn1_pe_end, &
                          atm2_pe_start, atm2_pe_end, &
                          ocn2_pe_start, ocn2_pe_end

   ! some local variables
      integer :: npes, pe, i, npes_atm1, npes_atm2, npes_ocn1, npes_ocn2

      check_parallel = .true.
      
      npes = mpp_npes()
      pe   = mpp_pe()

      allocate(pelist1(npes1), pelist2(npes2))
      pelist1 = (/(i, i=0, npes1-1)/)
      pelist2 = (/(i, i=npes1,npes-1)/)
      group1  = (pe .GE. 0) .AND. (pe .LT. npes1)
      group2  = (pe .GE. npes1) .AND. (pe .LT. npes)
      atm1_pe = (pe .GE. atm1_pe_start) .AND. (pe .LE. atm1_pe_end)
      atm2_pe = (pe .GE. atm2_pe_start) .AND. (pe .LE. atm2_pe_end)
      ocn1_pe = (pe .GE. ocn1_pe_start) .AND. (pe .LE. ocn1_pe_end)
      ocn2_pe = (pe .GE. ocn2_pe_start) .AND. (pe .LE. ocn2_pe_end)
      atm_pe  = atm1_pe .or. atm2_pe
      ocn_pe  = ocn1_pe .or. ocn2_pe
      npes_atm1 = atm1_pe_end - atm1_pe_start + 1
      npes_atm2 = atm2_pe_end - atm2_pe_start + 1
      npes_ocn1 = ocn1_pe_end - ocn1_pe_start + 1
      npes_ocn2 = ocn2_pe_end - ocn2_pe_start + 1
      allocate(atm_pelist(npes_atm1+npes_atm2), ocn_pelist(npes_ocn1+npes_ocn2))
      allocate(atm1_pelist(npes_atm1), atm2_pelist(npes_atm2), &
               ocn1_pelist(npes_ocn1), ocn2_pelist(npes_ocn2) )
      atm1_pelist = (/(i, i = atm1_pe_start, atm1_pe_end)/)
      atm2_pelist = (/(i, i = atm2_pe_start, atm2_pe_end)/)
      ocn1_pelist = (/(i, i = ocn1_pe_start, ocn1_pe_end)/)
      ocn2_pelist = (/(i, i = ocn2_pe_start, ocn2_pe_end)/)

      atm_pelist(1:npes_atm1) = (/(i, i = atm1_pe_start, atm1_pe_end)/)
      atm_pelist(npes_atm1+1:npes_atm1+npes_atm2) = (/(i, i = atm2_pe_start, atm2_pe_end)/)   
      ocn_pelist(1:npes_ocn1) = (/(i, i = ocn1_pe_start, ocn1_pe_end)/)
      ocn_pelist(npes_ocn1+1:npes_ocn1+npes_ocn2) = (/(i, i = ocn2_pe_start, ocn2_pe_end)/)  

      call mpp_declare_pelist(pelist1)
      call mpp_declare_pelist(pelist2)
      call mpp_declare_pelist(atm_pelist)
      call mpp_declare_pelist(ocn_pelist)

      return

   end subroutine mpp_set_pelist_init

!##########################################################################
!  release the memory
   subroutine mpp_set_pelist_end

      if(check_parallel)  &
          deallocate(pelist1, pelist2, atm_pelist, ocn_pelist)

       return

   end subroutine mpp_set_pelist_end

!############################################################################

end module mpp_set_pelist_mod
