!>simple routine to pass non-tracer fields across components (regrid)

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

public :: gex_init, gex_get_index,gex_get_n_ex, gex_get_property, gex_name, gex_units

character(3) :: module_name = 'gex'
logical      :: initialized = .FALSE.

integer, parameter :: gex_name  = 1
integer, parameter :: gex_units = 2

type gex_type
   character(fm_field_name_len):: name  = ''
   character(fm_string_len)    :: units = ''
   logical                     :: set   = .FALSE.
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

   if (mpp_pe()==mpp_root_pe()) then
      write(*,*) ''
      write(*,*) '####################################'
      write(*,*) '#  generic exchanged fields [gex]  #'
      write(*,*) '####################################'
      write(*,*) ''
   end if


   call gex_read_field_table('/coupler_mod/atm_to_lnd_ex',MODEL_ATMOS,MODEL_LAND)
   call gex_read_field_table('/coupler_mod/lnd_to_atm_ex',MODEL_LAND,MODEL_ATMOS)

   if (mpp_pe()==mpp_root_pe()) then
      write(*,*) ''
      write(*,*) '####################################'
      write(*,*) ''
   end if

   initialized = .TRUE.

end subroutine gex_init

!#######################################################################
!> Generic exchange between model components - process fields for a given exchange
!#######################################################################

subroutine gex_read_field_table(listroot,MODEL_SRC,MODEL_REC)

   character(len=*), intent(in) :: listroot  ! name of the field manager list
   integer,          intent(in) :: MODEL_SRC ! index of the model where the tracer comes FROM
   integer,          intent(in) :: MODEL_REC ! index of the model where the tracer goes TO

   type(fm_list_iter_type)      :: iter ! iterator over the list of tracers
   character(fm_field_name_len) :: name = '' ! name of the tracer
   character(fm_type_name_len)  :: ftype ! type of the field table entry (not used)
!    character(fm_path_name_len)  :: current_list ! storage for current location in the fiels manager tree
   character(fm_path_name_len)  :: listname  ! name of the field manager list for each tracer

   integer :: n

   if(.not.fm_dump_list(listroot, recursive=.TRUE.)) then
      call error_mesg('gex_read_field_table','Cannot dump field list "'//listroot//'". No additional field will be exchanged from land to atmosphere',NOTE)
      return
   endif

   n_gex(MODEL_SRC,MODEL_REC) = fm_get_length(listroot)
   allocate(gex_fields(MODEL_SRC,MODEL_REC)%field(n_gex(MODEL_SRC,MODEL_REC)))

   call fm_init_loop(listroot,iter)
   do while (fm_loop_over_list(iter, name, ftype, n))
      associate(fld=>gex_fields(MODEL_SRC,MODEL_REC)%field(n)) ! define a shorthand, to avoid very long expressions
         fld%name = trim(name)
         fld%set  = .FALSE.

         ! read parameters of the tracer

         ! I am not sure saving/restoring "current_list" is necessary: it seems to work
         ! without it. I just left it here commented out

         ! save current filed manager list
!             current_list = fm_get_current_list()
!             if (current_list .eq. ' ') &
!                 call error_mesg(module_name,'Could not get the current list',FATAL)

         ! switch to the list of tracer parameters
         listname = trim(listroot)//'/'//trim(name)
         if (.not.fm_change_list(listname)) then
            call error_mesg(module_name,'Cannot change fm list to "'//trim(listname)//'"', FATAL)
         endif
         ! read parameters
         fld%units = fm_util_get_string('units', caller = module_name, default_value = '', scalar = .true.)
         ! other parameters can be read here, for example:
         ! fld%molar_mass = get_spdata_real('molar_mass', caller = module_name, default_value=12.0, scalar=.true.)

         ! restore to the original list
!             if (.not.fm_change_list(current_list)) then
!                call error_mesg(module_name,'Cannot change fm list back to "'//trim(current_list)//'"', FATAL)
!             endif

         if (mpp_pe()==mpp_root_pe()) write(*,*) listroot,n,&
            ' name="'//trim(fld%name)//'"', &
            ' units="'//trim(fld%units)//'"'
      end associate
   end do
end subroutine


!#######################################################################
!> Generic exchange between model components - return number of fields exchanged
!#######################################################################

function gex_get_n_ex(MODEL_SRC,MODEL_REC)

   integer, intent(in)                         :: MODEL_SRC, MODEL_REC
   integer gex_get_n_ex

   gex_get_n_ex = n_gex(MODEL_SRC,MODEL_REC)

   return

end function

!#######################################################################
!> Generic exchange between model components - return name of field
!#######################################################################

function gex_get_property(MODEL_SRC,MODEL_REC,index,property)

   integer, intent(in)   :: MODEL_SRC, MODEL_REC,index
   integer               :: property
   character(len=64)     :: gex_get_property

   if (index.le.n_gex(MODEL_SRC,MODEL_REC)) then
      if (property .eq. gex_name) then
         gex_get_property = trim(gex_fields(MODEL_SRC,MODEL_REC)%field(index)%name)
      elseif (property .eq. gex_units) then
         gex_get_property = trim(gex_fields(MODEL_SRC,MODEL_REC)%field(index)%units)
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

function gex_get_index(MODEL_SRC,MODEL_REC,name,record)

   character(len=*), intent(in)                :: name !< name of the tracer
   integer, intent(in)                         :: MODEL_SRC, MODEL_REC
   logical, intent(in), optional               :: record    !record that this exchanged has been found and will be set

   integer :: i
   integer :: gex_get_index

   gex_get_index = NO_TRACER

   do i = 1, n_gex(MODEL_SRC,MODEL_REC)
      if (lowercase(trim(name)) == trim(gex_fields(MODEL_SRC,MODEL_REC)%field(i)%name))then
         gex_get_index = i

         if (present(record)) then
            if (record) then
               gex_fields(MODEL_SRC,MODEL_REC)%field(i)%set = .TRUE.
            end if
         else
            if (.not. gex_fields(MODEL_SRC,MODEL_REC)%field(i)%set) then
               call error_mesg('flux_exchange|gex','requested flux was never set',FATAL)
            end if
         end if

         exit
      endif
   enddo

   return

end function gex_get_index

end module gex_mod
