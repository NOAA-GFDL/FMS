
!>simple routine to pass non-tracer fields across components (regridding)

module gex_mod


use fms_mod,             only: lowercase, error_mesg, FATAL, NOTE
use tracer_manager_mod,  only: NO_TRACER
use field_manager_mod,   only: MODEL_LAND, MODEL_ATMOS, NUM_MODELS
use field_manager_mod,   only: fm_list_iter_type, fm_dump_list, fm_field_name_len, &
                               fm_type_name_len, fm_get_length,fm_loop_over_list, fm_init_loop, &
                               fm_string_len, fm_get_current_list, fm_path_name_len, fm_change_list
use fm_util_mod,         only: fm_util_get_real, fm_util_get_logical, fm_util_get_string
use mpp_mod,             only: mpp_root_pe, mpp_pe

implicit none ; private

public :: gex_init, gex_get_index,gex_get_n, gex_get_p, gex_name, gex_units

character(3) :: module_name = 'gex'
logical      :: initialized = .FALSE.

integer, parameter :: gex_name  = 1
integer, parameter :: gex_units = 2

type gex_type
   character(fm_field_name_len):: name
   character(fm_string_len)    :: units
end type gex_type
type gex_type_r
   type(gex_type), allocatable:: field(:)
end type gex_type_r  

integer, allocatable :: n_gex(:,:)
type(gex_type_r), allocatable :: gex_fields(:,:)

contains

!#######################################################################
!> Generic exchange between model components (initiatization)
!#######################################################################

subroutine gex_init()

   if (initialized) return

   allocate(n_gex(NUM_MODELS,NUM_MODELS))
   allocate(gex_fields(NUM_MODELS,NUM_MODELS))

   n_gex(:,:) = 0

   if (mpp_pe()==mpp_root_pe()) write(*,*) ''
   if (mpp_pe()==mpp_root_pe()) write(*,*) '####################################'
   if (mpp_pe()==mpp_root_pe()) write(*,*) '#  generic exchanged fields [gex]  #'
   if (mpp_pe()==mpp_root_pe()) write(*,*) '####################################'
   if (mpp_pe()==mpp_root_pe()) write(*,*) ''

   call gex_read_field_table('/coupler_mod/atm_lnd_ex',MODEL_ATMOS,MODEL_LAND)
   call gex_read_field_table('/coupler_mod/lnd_atm_ex',MODEL_LAND,MODEL_ATMOS)
   if (mpp_pe()==mpp_root_pe()) write(*,*) ''
   if (mpp_pe()==mpp_root_pe()) write(*,*) '####################################'
   if (mpp_pe()==mpp_root_pe()) write(*,*) ''

   initialized = .TRUE.

end subroutine gex_init

!#######################################################################
!> Generic exchange between model components - process fields for a given exchange
!#######################################################################   

subroutine gex_read_field_table(listroot,MODEL_SRC,MODEL_REC)

   integer, intent(in)                         :: MODEL_SRC,MODEL_REC
   character(len=*), intent(in)                :: listroot

   type(fm_list_iter_type) :: iter ! iterator over the list of species

   character(fm_field_name_len) :: name     = '' ! name of the species
   character(fm_type_name_len)  :: ftype ! type of the field table entry
   character(fm_path_name_len)  :: current_list ! storage for current location in the fiels manager tree
   character(fm_path_name_len)  :: listname  ! name of the field manager list

   integer :: n   

   if(fm_dump_list(listroot, recursive=.TRUE.)) then
      n_gex(MODEL_SRC,MODEL_REC) = fm_get_length(listroot)
      allocate(gex_fields(MODEL_SRC,MODEL_REC)%field(n_gex(MODEL_SRC,MODEL_REC)))
   
      call fm_init_loop(listroot,iter)
      do while (fm_loop_over_list(iter, name, ftype, n))
         gex_fields(MODEL_SRC,MODEL_REC)%field(n)%name = trim(name)
         if (mpp_pe()==mpp_root_pe()) write(*,*) listroot,n,trim(name)
   
         ! save current position in the field manager tree to restore it on exit
         current_list = fm_get_current_list()
         if (current_list .eq. ' ') call error_mesg(module_name,'Could not get the current list',FATAL)
      
         listname = trim(listroot)//'/'//trim(name)
      
         if (.not.fm_change_list(listname)) then
            call error_mesg(module_name,'Cannot change fm list to "'//trim(listname)//'"', FATAL)
         endif
    
         gex_fields(MODEL_SRC,MODEL_REC)%field(n)%units = &
                     fm_util_get_string('units',             &
                     caller = 'field_manager', default_value = '', scalar = .true.)          
                     
         if (.not.fm_change_list(current_list)) then
            call error_mesg(module_name,'Cannot change fm list back to "'//trim(current_list)//'"', FATAL)
         endif
      end do
   else
      call error_mesg('flux_exchange','Cannot dump field list "/coupler_mod/lnd_atm_ex". No additional field will be exchanged from land to atmosphere',NOTE)
   end if
   
end subroutine   


!#######################################################################
!> Generic exchange between model components - return number of fields exchanged
!#######################################################################   

function gex_get_n(MODEL_SRC,MODEL_REC)

   integer, intent(in)                         :: MODEL_SRC, MODEL_REC
   integer gex_get_n

   gex_get_n = n_gex(MODEL_SRC,MODEL_REC)

   return

end function

!#######################################################################
!> Generic exchange between model components - return name of field
!#######################################################################   

function gex_get_p(MODEL_SRC,MODEL_REC,index,property)

   integer, intent(in)   :: MODEL_SRC, MODEL_REC,index
   integer               :: property
   character(len=64)     :: gex_get_p

   if (index.le.n_gex(MODEL_SRC,MODEL_REC)) then
      if (property .eq. gex_name) then
         gex_get_p = trim(gex_fields(MODEL_SRC,MODEL_REC)%field(index)%name)
      elseif (property .eq. gex_units) then
         gex_get_p = trim(gex_fields(MODEL_SRC,MODEL_REC)%field(index)%units)
      else
         call error_mesg('flux_exchange|gex','property does not exist: '//gex_fields(MODEL_SRC,MODEL_REC)%field(index)%name,FATAL)
      end if
   else
      call error_mesg('flux_exchange|gex','requested tracer does not exist',FATAL)
   end if

   return

end function 

!#######################################################################
!> Generic exchange between model components - return index of exchange field
!#######################################################################      

function gex_get_index(MODEL_SRC,MODEL_REC,name)

   character(len=*), intent(in)                :: name !< name of the tracer                                                                                                                                          
   integer, intent(in)                         :: MODEL_SRC, MODEL_REC

   integer :: i
   integer :: gex_get_index

   gex_get_index = NO_TRACER

   do i = 1, n_gex(MODEL_SRC,MODEL_REC)
      if (lowercase(trim(name)) == trim(gex_fields(MODEL_SRC,MODEL_REC)%field(i)%name))then
         gex_get_index = i
         exit
      endif
   enddo

   return   

end function gex_get_index

end module gex_mod