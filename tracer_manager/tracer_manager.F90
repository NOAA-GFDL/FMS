module tracer_manager_mod
! <CONTACT EMAIL="William.Cooke@noaa.gov">
!   William Cooke
! </CONTACT>

! <REVIEWER EMAIL="Matthew.Harrison@noaa.gov">
!   Matt Harrison
! </REVIEWER>

! <REVIEWER EMAIL="Bruce.Wyman@noaa.gov">
!   Bruce Wyman
! </REVIEWER>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   Code to manage the simple addition of tracers to the FMS code.
!     This code keeps track of the numbers and names of tracers included
!     in a tracer table.
! </OVERVIEW>

! <DESCRIPTION>
!     This code is a grouping of calls which will allow the simple
!     introduction of tracers into the FMS framework. It is designed to
!     allow users of a variety of component models interact easily with
!     the dynamical core of the model. 
!     
!     In calling the tracer manager routines the user must provide a
!     parameter identifying the model that the user is working with. This
!     parameter is defined within field_manager as MODEL_X 
!     where X is one of [ATMOS, OCEAN, LAND, ICE].
!
!     In many of these calls the argument list includes model and tracer_index. These 
!     are the parameter corresponding to the component model and the tracer_index N is 
!     the Nth tracer within the component model. Therefore a call with MODEL_ATMOS and 5 
!     is different from a call with MODEL_OCEAN and 5.
!
! </DESCRIPTION>


!----------------------------------------------------------------------

use           mpp_mod, only : mpp_error,          &
                              mpp_pe,             &
                              mpp_root_pe,        &
                              FATAL,              &
                              WARNING,            &
                              NOTE,               &
                              stdlog
use        mpp_io_mod, only : mpp_open,           &
                              mpp_close,          &
                              MPP_ASCII,          &
                              MPP_APPEND,         &
                              MPP_RDONLY
use           fms_mod, only : lowercase,          &
                              write_version_number

use field_manager_mod, only : field_manager_init, &
                              get_field_info,     &
                              get_field_methods,  &
                              MODEL_ATMOS,        &
                              MODEL_LAND,         &
                              MODEL_OCEAN,        &
                              MODEL_ICE,          &
                              NUM_MODELS,         &
                              method_type,        &
                              default_method,     &
                              parse,              &
                              fm_new_list,        &
                              fm_copy_list,       &
                              fm_change_list,     &
                              fm_modify_name,     &
                              fm_change_root,     &
                              fm_dump_list,       &
                              fm_query_method,    &
                              fm_find_methods,    &
                              fm_get_length,      &
                              fm_new_value,       &
                              fm_exists

implicit none
private

!-----------------------------------------------------------------------

public  tracer_manager_init,       &
        tracer_manager_end,        &
        check_if_prognostic,       &
        assign_tracer_field,       &
        add_members_to_family,     &
        split_family_into_members, &
        find_family_members,       &
        get_family_name,           &
        get_tracer_indices,        &
        get_tracer_index,          &
        get_tracer_names,          &
        get_tracer_field,          &
        get_tracer_tendency,       &
        get_tracer_tlevels,        &
        query_combined,            &
        query_method,              &
        set_tracer_profile,        &
        register_tracers,          &
        set_tracer_atts,           &
        get_number_tracers,        &
        tracer_requires_init,      &
        have_initialized_tracer,   &
        query_tracer_init,         &
        NO_TRACER,                 &
        MAX_TRACER_FIELDS

!-----------------------------------------------------------------------

integer            :: num_tracer_fields = 0
integer, parameter :: MAX_TRACER_FIELDS = 100
integer, parameter :: MAX_TRACER_METHOD = 20
integer, parameter :: NO_TRACER         = -1

integer :: total_tracers(NUM_MODELS), prog_tracers(NUM_MODELS), diag_tracers(NUM_MODELS), family_tracers(NUM_MODELS)
logical :: model_registered(NUM_MODELS) = .FALSE.

type, private ::  tracer_type
   character(len=32)        :: tracer_name, tracer_units
   character(len=128)       :: tracer_longname, tracer_family
   integer                  :: num_methods, model, instances
   logical                  :: is_prognostic, is_family, is_combined, instances_set
   real, pointer, dimension(:,:,:,:) :: field_tlevels => NULL()
   real, pointer, dimension(:,:,:) :: field => NULL(), field_tendency => NULL(), weight => NULL()
   logical                  :: needs_init
end type tracer_type

type, private ::  tracer_name_type
   character (len=32)    :: model_name, tracer_name, tracer_units
   character (len=128)   :: tracer_longname, tracer_family
end type tracer_name_type


type, private :: inst_type
   character (len=128)   :: name
   integer               :: instances
end type inst_type

type(tracer_type), save  :: tracers(MAX_TRACER_FIELDS)
type(inst_type)  , save  :: instantiations(MAX_TRACER_FIELDS)

character(len=128) :: version = '$Id: tracer_manager.F90,v 12.0 2005/04/14 18:02:28 fms Exp $'
character(len=128) :: tagname = '$Name: lima $'
logical            :: module_is_initialized = .false.

logical            :: verbose_local
integer            :: TRACER_ARRAY(NUM_MODELS,MAX_TRACER_FIELDS)

contains

!
!#######################################################################
!
! <SUBROUTINE NAME="tracer_manager_init">
!   <OVERVIEW>
!     Routine to initialize the tracer manager
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine writes the version and tagname to the logfile and 
!     sets the module initialization flag.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call tracer_manager_init
!   </TEMPLATE>
subroutine tracer_manager_init
! version number to logfile

  call write_version_number (version, tagname)
  module_is_initialized = .TRUE.


end subroutine tracer_manager_init
! </SUBROUTINE>
!
!#######################################################################
!
! <SUBROUTINE NAME="register_tracers">

!   <OVERVIEW>
!      A routine to register the tracers included in a component model.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine returns the total number of valid tracers, the number of
! prognostic and diagnostic tracers and the number of families of
! tracers.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call register_tracers(model, num_tracers,num_prog,num_diag,num_family)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <OUT NAME="num_tracers" TYPE="integer">
!    The total number of valid tracers within the component model.
!   </OUT>
!   <OUT NAME="num_prog" TYPE="integer">
!     The number of prognostic tracers within the component model.
!   </OUT>
!   <OUT NAME="num_diag" TYPE="integer">
!     The number of diagnostic tracers within the component model.
!   </OUT>
!   <OUT NAME="num_family" TYPE="integer">
!     The number of family tracers within the component model.
!   </OUT>

subroutine register_tracers(model, num_tracers,num_prog,num_diag,num_family)
!
! read tracer table and store tracer information associated with "model"
! in "tracers" array. 

integer,  intent(in) :: model ! model being used
integer, intent(out) :: num_tracers, num_prog, num_diag, num_family
character(len=1024)  :: record
character(len=256)    :: warnmesg

character(len=32) :: model_name, name_type, type, name
character(len=128) :: units, longname, family
integer :: iunit,n,n1, m, ntf, mod, num_tracer_methods, nfields, swop
integer :: j, log_unit, num_in_family, num_methods, ns, num_tracer_comp_model
logical :: flag,is_family_member(MAX_TRACER_FIELDS), flag_type
type(method_type), dimension(MAX_TRACER_METHOD) :: methods
integer :: instances, siz_inst,i
character(len = 32) :: digit,suffnam

character(len=128) :: list_name , control
integer            :: index_list_name
logical :: fm_success


!   <ERROR MSG="invalid model type" STATUS="FATAL">
!     The index for the model type is invalid.
!   </ERROR>
if (model .ne. MODEL_ATMOS .and. model .ne. MODEL_LAND .and. &
    model .ne. MODEL_OCEAN .and. model .ne. MODEL_ICE) call mpp_error(FATAL,'register_tracers : invalid model type')

! One should only call register_tracers once for each model type
! Therefore need to set up an array to stop the subroutine being 
! unnecssarily called multiple times.

if ( model_registered(model) ) then
! This routine has already been called for the component model.
! Fill in the values from the previous registration and return.
  num_tracers = total_tracers(model)
  num_prog    = prog_tracers(model)
  num_diag    = diag_tracers(model) 
  num_family  = family_tracers(model)
  select case(model)
    case (MODEL_ATMOS)
      name = "atmospheric"
    case (MODEL_OCEAN)
      name = "oceanic"
    case (MODEL_ICE  )
      name = "ice"
    case (MODEL_LAND )
      name = "land"
    case default
      name = "ERROR"
    end select
if (mpp_pe() == mpp_root_pe()) &
  call mpp_error(NOTE,&
  'register_tracers : This routine has already been called for the '//trim(name)//' component model.')

  return
endif

! Initialize the number of tracers to zero.
num_tracers = 0; num_prog = 0; num_diag = 0; num_family = 0

call field_manager_init(nfields)
call tracer_manager_init


!   <ERROR MSG="No tracers are available to be registered." STATUS="NOTE">
!      No tracers are available to be registered. This means that the field
!      table does not exist or is empty.
!   </ERROR>
if (nfields == 0 ) then
if (mpp_pe() == mpp_root_pe()) &
  call mpp_error(NOTE,'register_tracers : No tracers are available to be registered.')
return
endif

! search through field entries for model tracers
num_tracer_comp_model = 0

do n=1,nfields
   call get_field_info(n,type,name,mod,num_methods)

   if (mod == model .and. type == 'tracer') then
      if (get_tracer_index(model, name) == NO_TRACER) then      
         num_tracer_fields = num_tracer_fields + 1
         num_tracer_comp_model = num_tracer_comp_model + 1
         TRACER_ARRAY(model,num_tracer_comp_model)  = num_tracer_fields
!   <ERROR MSG="MAX_TRACER_FIELDS exceeded" STATUS="FATAL">
!     The maximum number of tracer fields has been exceeded.
!   </ERROR>
         if(num_tracer_fields > MAX_TRACER_FIELDS) call mpp_error(FATAL,'register_tracers: MAX_TRACER_FIELDS exceeded')
         tracers(num_tracer_fields)%model          = model
         tracers(num_tracer_fields)%tracer_name    = name
         tracers(num_tracer_fields)%tracer_units   = 'none'
         tracers(num_tracer_fields)%tracer_longname = tracers(num_tracer_fields)%tracer_name
         tracers(num_tracer_fields)%tracer_family   = 'orphan'
         tracers(num_tracer_fields)%instances_set   = .FALSE.
         num_tracer_methods     = 0
         methods = default_method ! initialize methods array
         call get_field_methods(n,methods)
         do j=1,num_methods
            select case (methods(j)%method_type) 
            case ('units')
               tracers(num_tracer_fields)%tracer_units   = methods(j)%method_name
            case ('longname')
               tracers(num_tracer_fields)%tracer_longname = methods(j)%method_name
            case ('family')
               tracers(num_tracer_fields)%tracer_family = methods(j)%method_name
            case ('instances')
!               tracers(num_tracer_fields)%instances = methods(j)%method_name
               siz_inst = parse(methods(j)%method_name,"",instances)
               tracers(num_tracer_fields)%instances = instances
               tracers(num_tracer_fields)%instances_set   = .TRUE.
            case default
               num_tracer_methods = num_tracer_methods+1
!               tracers(num_tracer_fields)%methods(num_tracer_methods) = methods(j)
            end select
         enddo
         tracers(num_tracer_fields)%num_methods = num_tracer_methods
         tracers(num_tracer_fields)%needs_init = .false.
         flag_type = query_method ('tracer_type',model,num_tracer_comp_model,name_type)
         if (flag_type .and. name_type == 'diagnostic') then
            tracers(num_tracer_fields)%is_prognostic = .false.
         else   
            tracers(num_tracer_fields)%is_prognostic = .true.
         endif   
         tracers(num_tracer_fields)%is_family = .false.
         tracers(num_tracer_fields)%is_combined = .false.
         if (tracers(num_tracer_fields)%is_prognostic) then
            num_prog = num_prog+1
         else
            num_diag = num_diag+1
         endif
      endif
   endif
enddo

! Now cycle through the tracers and add additional instances of the tracers.

ntf = num_tracer_fields
do n = 1, ntf !{
!   call get_field_info(n,type,name,mod,num_methods)

  if ( model == tracers(n)%model .and. tracers(n)%instances_set ) then !{ We have multiple instances of this tracer

    if ( ntf + tracers(n)%instances > MAX_TRACER_FIELDS ) then
      write(warnmesg, '("register_tracer : Number of tracers will exceed MAX_TRACER_FIELDS with &
                       &multiple (",I3," instances) setup of tracer ",A)') tracers(n)%instances,tracers(n)%tracer_name
      call mpp_error(FATAL, warnmesg)
    endif                        

    do i = 2, tracers(n)%instances !{
      num_tracer_fields = num_tracer_fields + 1
      num_tracer_comp_model = num_tracer_comp_model + 1
      TRACER_ARRAY(model,num_tracer_comp_model)  = num_tracer_fields
      ! Copy the original tracer type to the multiple instances.
      tracers(num_tracer_fields) = tracers(n)
      if ( query_method ('instances', model,model_tracer_number(model,n),name, control)) then !{
          
        if (i .lt. 10) then  !{
           write (suffnam,'(''suffix'',i1)') i
           siz_inst = parse(control, suffnam,digit)
           if (siz_inst == 0 ) then
             write (digit,'(''_'',i1)') i
           else
             digit = "_"//trim(digit)
           endif  
        elseif (i .lt. 100) then  !}{
           write (suffnam,'(''suffix'',i2)') i
           siz_inst = parse(control, suffnam,digit)
           if (siz_inst == 0 ) then
             write (digit,'(''_'',i2)') i
           else
             digit = "_"//trim(digit)
           endif
        else  !}{
          call mpp_error(FATAL, 'register_tracer : MULTIPLE_TRACER_SET_UP exceeds 100 for '//tracers(n)%tracer_name )
        endif  !}

        select case(model)
          case (MODEL_ATMOS)
            list_name = "/atmos_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
          case (MODEL_OCEAN)
            list_name = "/ocean_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
          case (MODEL_ICE  )
            list_name = "/ice_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
          case (MODEL_LAND )
            list_name = "/land_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
          case default
            list_name = "/default/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
        end select

        if (mpp_pe() == mpp_root_pe() ) write (*,*) "Creating list name = ",trim(list_name)//trim(digit)

        index_list_name = fm_copy_list(trim(list_name),digit, create = .true.)
        tracers(num_tracer_fields)%tracer_name = trim(tracers(num_tracer_fields)%tracer_name)//trim(digit)
      endif !}
         
      if (tracers(num_tracer_fields)%is_prognostic) then !{
         num_prog = num_prog+1
      else !}{
         num_diag = num_diag+1
      endif !}
    enddo !}
    ! Multiple instances of tracers were found so need to rename the original tracer.
    digit = "_1" 
    siz_inst = parse(control, "suffix1",digit)
    if (siz_inst > 0 ) then !{
      digit = "_"//trim(digit)
    endif !}
    fm_success = fm_modify_name(trim(list_name), trim(tracers(n)%tracer_name)//trim(digit))
    tracers(n)%tracer_name = trim(tracers(n)%tracer_name)//trim(digit)
  endif !}
enddo !}

!Now recycle through the tracers and get the family names

ntf = num_tracer_fields
do n=1,ntf

  if ( model == tracers(n)%model .and. tracers(n)%tracer_family /= 'orphan') then
      call find_family_members(tracers(n)%model,tracers(n)%tracer_family,is_family_member)
      num_in_family = 0
      do m=1,ntf
         if (is_family_member(m)) num_in_family = num_in_family+1
      end do
      if (num_in_family < 2) then ! do not set up a family tracer
         write(warnmesg,911) num_in_family,trim(tracers(n)%tracer_family)
911      format('register_tracers : There is only ',i2,' tracers for tracer family ', (a),'. Making an orphan.')
!   <ERROR MSG="There is only 1 tracer for tracer family X. Making an orphan." STATUS="NOTE">
!     A tracer has been given a family name but that family has only this member. Therefore it should be an orphan.
!   </ERROR>
         if (mpp_pe() == mpp_root_pe()) call mpp_error(NOTE,warnmesg)
         tracers(n)%tracer_family = 'orphan'
         cycle
      else
         m = get_tracer_index(tracers(n)%model,tracers(n)%tracer_family)
         if (m < 0) then 
         num_tracer_fields = num_tracer_fields + 1
         num_tracer_comp_model = num_tracer_comp_model + 1
         TRACER_ARRAY(model,num_tracer_comp_model)  = num_tracer_fields
            num_family = num_family+1
!   <ERROR MSG="MAX_TRACER_FIELDS needs to be increased" STATUS="FATAL">
!     The number of tracer fields has exceeded the maximum allowed. 
!     The parameter MAX_TRACER_FIELDS needs to be increased.
!   </ERROR>
            if (num_tracer_fields .gt. MAX_TRACER_FIELDS) &
                call mpp_error(FATAL,'register_tracers : MAX_TRACER_FIELDS needs to be increased')
            if (mpp_pe() == mpp_root_pe()) write(*,*) 'defining new tracer family: ',trim(tracers(n)%tracer_family)
            tracers(num_tracer_fields)%is_family = .true.
            tracers(num_tracer_fields)%tracer_name = trim(tracers(n)%tracer_family)
            tracers(num_tracer_fields)%model = tracers(n)%model
            tracers(num_tracer_fields)%tracer_units = tracers(n)%tracer_units
            tracers(num_tracer_fields)%tracer_longname = tracers(num_tracer_fields)%tracer_name
            tracers(num_tracer_fields)%tracer_family = 'orphan'
            num_methods = tracers(n)%num_methods
            tracers(num_tracer_fields)%num_methods = num_methods
            call check_family_parameters(model,tracers(n)%tracer_name) ! make sure family parameters are same for tracers within family
         endif
      endif
   endif
enddo


! Find any field entries with the instances keyword.
do n=1,nfields
   call get_field_info(n,type,name,mod,num_methods)

   if ( mod == model .and. type == 'instances' ) then
      call get_field_methods(n,methods)
      do j=1,num_methods

         m = get_tracer_index(mod,methods(j)%method_type)
         if (m .eq. NO_TRACER) then 
           call mpp_error(FATAL,'register_tracers : The instances keyword was found for undefined tracer '&
           //trim(methods(j)%method_type))
         else
           if ( tracers(m)%instances_set ) &
              call mpp_error(FATAL,'register_tracers : The instances keyword was found for '&
              //trim(methods(j)%method_type)//' but has previously been defined in the tracer entry')
           siz_inst = parse(methods(j)%method_name,"",instances)
           tracers(m)%instances = instances
           call mpp_error(NOTE,'register_tracers : '//trim(instantiations(j)%name)// &
                               ' will have '//trim(methods(j)%method_name)//' instances')
         endif
         if ( ntf + instances > MAX_TRACER_FIELDS ) then
           write(warnmesg, '("register_tracer : Number of tracers will exceed MAX_TRACER_FIELDS with &
                       &multiple (",I3," instances) setup of tracer ",A)') tracers(m)%instances,tracers(m)%tracer_name
           call mpp_error(FATAL, warnmesg)
         endif                        
! We have found a valid tracer that has more than one instantiation.
! We need to modify that tracer name to tracer_1 and add extra tracers for the extra instantiations.
         if (instances .eq. 1) then
           siz_inst = parse(methods(j)%method_control, 'suffix1',digit)
           if (siz_inst == 0 ) then
             digit = '_1'
           else
             digit = "_"//trim(digit)
           endif  
         endif
         do i = 2, instances
           num_tracer_fields = num_tracer_fields + 1
           num_tracer_comp_model = num_tracer_comp_model + 1
           TRACER_ARRAY(model,num_tracer_comp_model)  = num_tracer_fields
           tracers(num_tracer_fields)                =  tracers(m)
           
           
           
           if (i .lt. 10) then  !{
             write (suffnam,'(''suffix'',i1)') i
             siz_inst = parse(methods(j)%method_control, suffnam,digit)
             if (siz_inst == 0 ) then
               write (digit,'(''_'',i1)') i
             else
               digit = "_"//trim(digit)
             endif  
          elseif (i .lt. 100) then  !}{
             write (suffnam,'(''suffix'',i2)') i
             siz_inst = parse(methods(j)%method_control, suffnam,digit)
             if (siz_inst == 0 ) then
               write (digit,'(''_'',i2)') i
             else
               digit = "_"//trim(digit)
             endif
          else  !}{
            call mpp_error(FATAL, 'register_tracer : MULTIPLE_TRACER_SET_UP exceeds 100 for '&
                                  //tracers(num_tracer_fields)%tracer_name )
          endif  !}

          select case(model)
            case (MODEL_ATMOS)
              list_name = "/atmos_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
            case (MODEL_OCEAN)
              list_name = "/ocean_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
            case (MODEL_ICE  )
              list_name = "/ice_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
            case (MODEL_LAND )
              list_name = "/land_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
            case default
              list_name = "/default/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
          end select

          if (mpp_pe() == mpp_root_pe() ) write (*,*) "Creating list name = ",trim(list_name)

          index_list_name = fm_copy_list(trim(list_name),digit, create = .true.)

          tracers(num_tracer_fields)%tracer_name    =  trim(tracers(num_tracer_fields)%tracer_name)//digit
          if (tracers(num_tracer_fields)%is_prognostic) then
            num_prog = num_prog+1
          else
            num_diag = num_diag+1
          endif
        enddo
!Now rename the original tracer to tracer_1 (or if suffix1 present to tracer_'value_of_suffix1')
        siz_inst = parse(methods(j)%method_control, 'suffix1',digit)
        if (siz_inst == 0 ) then
          digit = '_1'
        else
          digit = "_"//trim(digit)
        endif  
        fm_success = fm_modify_name(trim(list_name), trim(tracers(m)%tracer_name)//trim(digit))
        tracers(m)%tracer_name    =  trim(tracers(m)%tracer_name)//trim(digit)
      enddo
   endif
enddo



num_tracers = num_prog + num_diag + num_family
! Make the number of tracers available publicly.
total_tracers(model)    = num_tracers
prog_tracers(model)     = num_prog
diag_tracers(model)     = num_diag
family_tracers(model)   = num_family
model_registered(model) = .TRUE.
!   <ERROR MSG="Families of tracers should be used with great caution." STATUS="NOTE">
!     Families of tracers were originally implemented in order to allow the advection of 
!     a family as one tracer. Unless the spatial profile of the family members is
!     similar, there will be leakage of one tracer to it's siblings. Therefore one 
!     should be very careful when utilizing families of tracers.
!   </ERROR>
if (num_family > 0 ) &
  call mpp_error(NOTE,"register_tracers : Families of tracers should be used with great caution.")

! Now sort through the tracer fields and sort them so that the 
! prognostic tracers are first. This should include the family tracers.

do n=1, num_tracers
  if (.not.check_if_prognostic(model,n) .and. n.le.(num_prog+num_family)) then 
  ! This is a diagnostic tracer so find a prognostic tracer to swop with
    do m = n, num_tracers
       if (check_if_prognostic(model,m) .and. .not.check_if_prognostic(model,n)) then
           swop = TRACER_ARRAY(model,n)
           TRACER_ARRAY(model,n) = TRACER_ARRAY(model,m)
           TRACER_ARRAY(model,m) = swop
           cycle
       endif
    enddo
  endif
enddo


do n=1, num_tracer_fields
   if(mpp_pe()==mpp_root_pe() .and. TRACER_ARRAY(model,n)> 0 ) &
      call print_tracer_info(TRACER_ARRAY(model,n))
enddo


log_unit = stdlog()
if ( mpp_pe() == mpp_root_pe() ) then
    select case (model)
      case (MODEL_ATMOS)
         model_name = "atmospheric"
      case (MODEL_OCEAN)
         model_name = "oceanic"
      case (MODEL_ICE)
         model_name = "ice"
      case (MODEL_LAND)
         model_name = "land"
      case default
    end select

   write (log_unit,15) trim(model_name),total_tracers(model)
endif

15 format ('Number of tracers in field table for ',A,' model = ',i4)


end subroutine register_tracers
!</SUBROUTINE>


function model_tracer_number(model,n)
integer, intent(in) :: model, n
integer model_tracer_number

integer :: i

model_tracer_number = NO_TRACER

do i = 1, MAX_TRACER_FIELDS
  if ( TRACER_ARRAY(model,i) == n ) &
     model_tracer_number = i

enddo

end function model_tracer_number

! <SUBROUTINE NAME="tracer_requires_init">
!   <OVERVIEW>
!      A routine to be called to alert the user that the tracer field requires 
!   initialization.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine sets a flag within the tracer_type which can be queried by
!     the user in the user tracer initialization scheme to see if initialization 
!     of the tracer field is necessary. The flag is set true if a restart file
!     for the tracer is not found.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call tracer_requires_init(model, name)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <IN NAME="name" TYPE="character">
!    The name of the tracer within the component model that requires initializing.
!   </IN>
subroutine tracer_requires_init(model, name)
integer, intent(in)             :: model
character (len = *), intent(in) :: name 

integer :: tr_init

tr_init = get_tracer_index(model, name)
tracers(tr_init)%needs_init = .true.

end subroutine tracer_requires_init
!</SUBROUTINE>


! <SUBROUTINE NAME="have_initialized_tracer">
!   <OVERVIEW>
!      A routine to show that the user has initialized the tracer.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine sets a flag within the tracer_type to show that the user has
!     initialized the tracer field.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call have_initialized_tracer(model, name)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <IN NAME="name" TYPE="character">
!    The name of the tracer within the component model that requires initializing.
!   </IN>
subroutine have_initialized_tracer(model, name)
integer, intent(in)             :: model
character (len = *), intent(in) :: name 

integer :: tr_init

tr_init = get_tracer_index(model, name)
tracers(tr_init)%needs_init = .false.

end subroutine have_initialized_tracer
!</SUBROUTINE>

! <FUNCTION NAME="query_tracer_init">
!   <OVERVIEW>
!      A function to be called to alert the user that the tracer field requires 
!   initialization.
!   </OVERVIEW>
!   <DESCRIPTION>
!     A function to query the initializtion status of a tracer.
!     The return value is true if the tracer requires initialization 
!     i.e. a restart file for the tracer is not found.
!   </DESCRIPTION>
!   <TEMPLATE>
!     if (query_tracer_init(model, name) ) then
!        tracer initialization code
!     endif
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <IN NAME="name" TYPE="character">
!    The name of the tracer within the component model that requires initializing.
!   </IN>
function query_tracer_init(model, name)
integer, intent(in)             :: model
character (len = *), intent(in) :: name 
logical                         :: query_tracer_init

integer :: tr_init

tr_init = get_tracer_index(model, name)
query_tracer_init = tracers(tr_init)%needs_init

end function query_tracer_init
!</FUNCTION>

! <SUBROUTINE NAME="get_number_tracers">
!   <OVERVIEW>
!      A routine to return the number of tracers included in a component model.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine returns the total number of valid tracers, the number of
! prognostic and diagnostic tracers and the number of families of
! tracers.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_number_tracers(model, num_tracers,num_prog,num_diag,num_family)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <OUT NAME="num_tracers" TYPE="integer, optional">
!    The total number of valid tracers within the component model.
!   </OUT>
!   <OUT NAME="num_prog" TYPE="integer, optional">
!     The number of prognostic tracers within the component model.
!   </OUT>
!   <OUT NAME="num_diag" TYPE="integer, optional">
!     The number of diagnostic tracers within the component model.
!   </OUT>
!   <OUT NAME="num_family" TYPE="integer, optional">
!     The number of family tracers within the component model.
!   </OUT>
subroutine get_number_tracers(model, num_tracers , num_prog, num_diag, num_family)

integer,  intent(in) :: model
integer, intent(out), optional :: num_tracers, num_prog, num_diag, num_family

!   <ERROR MSG="Model number is invalid." STATUS="FATAL">
!     The index of the component model is invalid.
!   </ERROR>
if (model .ne. MODEL_ATMOS .and. model .ne. MODEL_LAND .and. &
    model .ne. MODEL_OCEAN .and. model .ne. MODEL_ICE)  &
    call mpp_error(FATAL,"get_number_tracers : Model number is invalid.")

if (present(num_tracers)) num_tracers = total_tracers(model)
if (present(num_prog))    num_prog    = prog_tracers(model)
if (present(num_diag))    num_diag    = diag_tracers(model)
if (present(num_family))  num_family  = family_tracers(model)

end subroutine get_number_tracers
!</SUBROUTINE>


! <SUBROUTINE NAME="get_tracer_indices">

!   <OVERVIEW>
!     Routine to return the component model tracer indices as defined within
!     the tracer manager.
!   </OVERVIEW>
!   <DESCRIPTION>
!     If several models are being used or redundant tracers have been written to
! the tracer_table, then the indices in the component model and the tracer
! manager may not have a one to one correspondence. Therefore the component
! model needs to know what index to pass to calls to tracer_manager routines in
! order that the correct tracer information be accessed.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_tracer_indices(model, ind, prog_ind, diag_ind, fam_ind)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <OUT NAME="ind" TYPE="integer, optional" DIM="(:)" >
! An array containing the tracer manager defined indices for
!             all the tracers within the component model.
!   </OUT>
!   <OUT NAME="prog_ind" TYPE="integer, optional" DIM="(:)" >
! An array containing the tracer manager defined indices for
!             the prognostic tracers within the component model.
!   </OUT>
!   <OUT NAME="diag_ind" TYPE="integer, optional" DIM="(:)" >
! An array containing the tracer manager defined indices for
!             the diagnostic tracers within the component model.
!   </OUT>
!   <OUT NAME="fam_ind" TYPE="integer, optional" DIM="(:)" >
! An array containing the tracer manager defined indices for
!             the family tracers within the component model.
!   </OUT>
subroutine get_tracer_indices(model, ind, prog_ind, diag_ind, fam_ind)

integer, intent(in) :: model
integer, intent(out), dimension(:), optional :: ind, prog_ind, diag_ind, fam_ind

integer :: i, j, nf, np, nd, n

nf=0;nd=0;np=0;n=0

! Initialize arrays with dummy values
if (PRESENT(ind))      ind      = NO_TRACER
if (PRESENT(prog_ind)) prog_ind = NO_TRACER
if (PRESENT(diag_ind)) diag_ind = NO_TRACER
if (PRESENT(fam_ind))  fam_ind  = NO_TRACER

do i = 1, MAX_TRACER_FIELDS
j = TRACER_ARRAY(model,i)
 if ( j .gt. 0) then
   if ( model == tracers(j)%model) then
      if (PRESENT(ind)) then
         n=n+1
!   <ERROR MSG="index array size too small in get_tracer_indices" STATUS="FATAL">
!     The global index array is too small and cannot contain all the tracer numbers.
!   </ERROR>
         if (n > size(ind(:))) call mpp_error(FATAL,'get_tracer_indices : index array size too small in get_tracer_indices')
         ind(n) = i
      endif

      if (tracers(j)%is_family.and.PRESENT(fam_ind)) then
         nf=nf+1
!   <ERROR MSG="family array size too small in get_tracer_indices" STATUS="FATAL">
!     The family index array is too small and cannot contain all the tracer numbers.
!   </ERROR>
         if (nf > size(fam_ind(:))) call mpp_error(FATAL,'get_tracer_indices : family array size too small in get_tracer_indices')
         fam_ind(nf) = i
         cycle
      endif

      if (tracers(j)%is_prognostic.and.PRESENT(prog_ind)) then
         np=np+1
!   <ERROR MSG="prognostic array size too small in get_tracer_indices" STATUS="FATAL">
!     The prognostic index array is too small and cannot contain all the tracer numbers.
!   </ERROR>
         if ( np > size( prog_ind(:)))call mpp_error(FATAL,&
                                          'get_tracer_indices : prognostic array size too small in get_tracer_indices')
         prog_ind(np) = i
      else if (.not.tracers(j)%is_prognostic .and. .not. tracers(j)%is_family &
             .and.PRESENT(diag_ind)) then
         nd = nd+1
!   <ERROR MSG="diagnostic array size too small in get_tracer_indices" STATUS="FATAL">
!     The diagnostic index array is too small and cannot contain all the tracer numbers.
!   </ERROR>
         if (nd > size(diag_ind(:))) call mpp_error(FATAL,&
                                         'get_tracer_indices : diagnostic array size too small in get_tracer_indices')
         diag_ind(nd) = i
      endif
   endif
 endif
enddo

return
end subroutine get_tracer_indices
!</SUBROUTINE>

!<FUNCTION NAME= "get_tracer_index">
!   <OVERVIEW>
!     Function which returns the number assigned to the tracer name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This is a function which returns the index, as implied within the component model.
!   </DESCRIPTION>
!   <TEMPLATE>
!     value=get_tracer_index(model, name, indices, verbose)
!   </TEMPLATE>
!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <IN NAME="name" TYPE="character">
!     The name of the tracer (as assigned in the field table).
!   </IN>
!   <IN NAME="indices" TYPE="integer, optional" DIM="(:)">
!     An array of the component model indices. This array can be found by 
!            calling get_tracer_indices.
!   </IN>
!   <IN NAME="verbose" TYPE="logical, optional">
!     A flag to allow the message saying that a tracer with this name has not 
!     been found. This should only be used for debugging purposes.
!   </IN>
!   <OUT NAME="get_tracer_index" TYPE="integer">
!     The index of the tracer named "name". 
!     If indices is passed then the result is the array index which 
! corresponds to tracer named "name".
!   </OUT>
function get_tracer_index(model, name, indices, verbose)

integer, intent(in)                         :: model
character(len=*), intent(in)                :: name
integer, intent(in), dimension(:), optional :: indices
logical, intent(in), optional               :: verbose
integer :: get_tracer_index

integer :: i

get_tracer_index = NO_TRACER

if (PRESENT(indices)) then
    do i = 1, size(indices(:))
       if (model == tracers(indices(i))%model .and. lowercase(trim(name)) == trim(tracers(indices(i))%tracer_name)) then
           get_tracer_index = i
           exit
       endif
    enddo
else
    do i=1, num_tracer_fields
       if (lowercase(trim(name)) == trim(tracers(TRACER_ARRAY(model,i))%tracer_name)) then
           get_tracer_index = i!TRACER_ARRAY(model,i)
           exit
       endif
    enddo
end if



verbose_local=.FALSE.
if (present(verbose)) verbose_local=verbose

if (verbose_local) then
! <ERROR MSG="tracer with this name not found: X" STATUS="NOTE">
  if (get_tracer_index == NO_TRACER ) call mpp_error(NOTE,'get_tracer_index : tracer with this name not found: '//trim(name))
! </ERROR>
endif
   
return

end function get_tracer_index
!</FUNCTION>


! <SUBROUTINE NAME="assign_tracer_field" >
!   <OVERVIEW>
!     Routine to point the appropriate field within the tracer_type to the 
! appropriate field within the component model.
!   </OVERVIEW>
!   <DESCRIPTION>
!     The generality provided here is that one can point the three
! dimensional tracer field at either a two time level scheme [data and
! tendency] or a three time level scheme [data_tlevels]. The tracer manager
! points the appropriate tracer_type field at the data supplied from the 
! component model.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call assign_tracer_field(model,index, data, data_tlevels, tendency)
!   </TEMPLATE>
!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="index" TYPE="integer">
!     The tracer number that you wish to assign a tracer
!                  field for.
!   </IN>
!   <IN NAME="data" TYPE="real, target, optional" DIM="(:,:,:)" >
!     The 3D field that is associated with the present time 
!                  step in the component model.
!   </IN>
!   <IN NAME="tendency" TYPE="real, target, optional" DIM="(:,:,:)" >
!     The 3D field that is associated with the tendency time
!                  step in the component model.
!   </IN>
!   <IN NAME="data_tlevels" TYPE="real, target, optional" DIM="(:,:,:,:)" >
!     The 4D field that is associated with the tracer field 
!                  in the component model.
!   </IN>
subroutine assign_tracer_field(model,index, data, data_tlevels, tendency)

integer, intent(in)                        :: model, index
real, intent(in), dimension(:,:,:), target, optional   :: data, tendency
real, intent(in), dimension(:,:,:,:), target, optional :: data_tlevels

integer :: check

!   <ERROR MSG="invalid index" STATUS="FATAL">
!     The index that has been passed to this routine is invalid.
!   </ERROR>
if (index < 0 .or. index > num_tracer_fields) call mpp_error(FATAL,'assign_tracer_field : invalid index')

if (PRESENT(data)) tracers(TRACER_ARRAY(model,index))%field => data
if (PRESENT(data_tlevels)) tracers(TRACER_ARRAY(model,index))%field_tlevels => data_tlevels
if (PRESENT(tendency)) tracers(TRACER_ARRAY(model,index))%field_tendency => tendency


check =0
if (PRESENT(data)) check = check + 1
if (PRESENT(data_tlevels)) check = check + 1
if (PRESENT(tendency)) check = check + 1

!   <ERROR MSG="At least one of data, data_tlevels or tendency must be passed in here." STATUS="FATAL">
!     At least one of data, data_tlevels or tendency must be passed to assign_tracer_field
!     Otherwise there is not much point in calling this routine.
!   </ERROR>
if (check == 0) call mpp_error(FATAL,'assign_tracer_field : At least one of data, data_tlevels or tendency must be passed in here.')

return
end subroutine assign_tracer_field
!</SUBROUTINE>

!
!#######################################################################
!

! <SUBROUTINE NAME="tracer_manager_end" >
!   <OVERVIEW>
!     Routine to write to the log file that the tracer manager is ending.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Routine to write to the log file that the tracer manager is ending.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call tracer_manager_end
!   </TEMPLATE>

subroutine tracer_manager_end

integer :: log_unit

log_unit = stdlog()
if ( mpp_pe() == mpp_root_pe() ) then
   write (log_unit,'(/,(a))') 'Exiting tracer_manager, have a nice day ...'
endif

module_is_initialized = .FALSE.

end subroutine tracer_manager_end
!</SUBROUTINE>


!
!#######################################################################
!
subroutine print_tracer_info(i)
!
! Routine to print out the components of the tracer.
! This is useful for informational purposes.
! Used in register_tracers.
!
! Arguments:
! INTENT IN
!  i            : index of the tracer that is being printed.
!
integer, intent(in) :: i
integer :: j,log_unit

character(len=80) :: name,control
logical   :: flag


log_unit = stdlog()
write(log_unit, *)'----------------------------------------------------'
write(log_unit, *) 'Contents of tracer entry ', i
write(log_unit, *) 'Model type and field name'
write(log_unit, *) 'Model                : ', tracers(i)%model
write(log_unit, *) 'Field name           : ', trim(tracers(i)%tracer_name)
write(log_unit, *) 'Tracer family        : ', trim(tracers(i)%tracer_family)
write(log_unit, *) 'Tracer units         : ', trim(tracers(i)%tracer_units)
write(log_unit, *) 'Tracer longname      : ', trim(tracers(i)%tracer_longname)
write(log_unit, *) 'Tracer is_prognostic : ', tracers(i)%is_prognostic
write(log_unit, *) 'Tracer is_family     : ', tracers(i)%is_family
write(log_unit, *) 'Tracer is_combined   : ', tracers(i)%is_combined

if (associated(tracers(i)%field)) then
   write(log_unit, '(a,3i5)') 'Size of tracer field : ', size(tracers(i)%field,1), size(tracers(i)%field,2),&
                                                         size(tracers(i)%field,3)
else
   write(log_unit,*) 'Tracer field         : not associated'
endif

if (associated(tracers(i)%field_tendency)) then
   write(log_unit, '(a,3i5)') 'Size of tracer tendency: ', size(tracers(i)%field_tendency,1), size(tracers(i)%field_tendency,2),&
                                                           size(tracers(i)%field_tendency,3)
else
   write(log_unit,*) 'Tracer tendency      : not associated'
endif

if (associated(tracers(i)%field_tlevels)) then
   write(log_unit, '(a,3i5)') 'Size of tracer tlevels : ', size(tracers(i)%field_tlevels,1), size(tracers(i)%field_tlevels,2),&
                                                           size(tracers(i)%field_tlevels,3), size(tracers(i)%field_tlevels,4)
else
   write(log_unit,*) 'Tracer tlevels       : not associated'
endif


write(log_unit, *)'----------------------------------------------------'

900 FORMAT(A,2(1x,E12.6))
901 FORMAT(E12.6,1x,E12.6)


end subroutine print_tracer_info

!
!#######################################################################
!
!<FUNCTION NAME= "get_tracer_field">
!   <OVERVIEW>
!     A function to retrieve the present timestep data.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Function to point to the 3D field associated with a tracer.
!   </DESCRIPTION>
!   <TEMPLATE>
!     array=get_tracer_field(model, tracer_index)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="tracer_index" TYPE="integer">
!     The tracer number within the component model.
!   </IN>
!   <OUT NAME="data"  TYPE="real, pointer" DIM="(:,:,:)">
!     The tracer field is returned in this array.
!   </OUT>
function get_tracer_field(model, tracer_index) result (data)

integer              :: model, tracer_index
real, pointer        :: data(:,:,:)

integer :: n

!Convert local model index to tracer_manager index
!   <ERROR MSG="invalid index" STATUS="FATAL">
!     The index that has been passed to this routine is invalid.
!          Check the index that is being passed corresponds to a valid
!          tracer name.
!   </ERROR>
if (TRACER_ARRAY(model,tracer_index) < 1 .or. TRACER_ARRAY(model,tracer_index) > num_tracer_fields) &
    call mpp_error(FATAL,'get_tracer_field : invalid index')
!   <ERROR MSG="tracer field array not allocated" STATUS="FATAL">
!      The tracer array has not been allocated. This means that a
!          call to assign_tracer_field is absent in the code.
!   </ERROR>
if (.not. associated(tracers(TRACER_ARRAY(model,tracer_index))%field)) &
    call mpp_error(FATAL,'get_tracer_field : tracer field array not allocated')
data =>  tracers(TRACER_ARRAY(model,tracer_index))%field

end function get_tracer_field
!</FUNCTION>
!
!#######################################################################
!
!<FUNCTION NAME= "get_tracer_tlevels">
!   <OVERVIEW>
!     A function to retrieve the three time levels  data.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Function to point to the 4D field associated with a tracer.
!   </DESCRIPTION>
!   <TEMPLATE>
!     array=get_tracer_tlevels(model, tracer_index)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="tracer_index" TYPE="integer">
!     The tracer number within the component model.
!   </IN>
!   <OUT NAME="data"  TYPE="real, pointer" DIM="(:,:,:,:)">
!     The tracer field is returned in this array.
!   </OUT>
function get_tracer_tlevels(model, tracer_index) result (data)

integer  :: model, tracer_index
real, pointer        :: data(:,:,:,:)

integer :: n

!Convert local model index to tracer_manager index
!   <ERROR MSG="invalid index" STATUS="FATAL">
!     The index that has been passed to this routine is invalid.
!          Check the index that is being passed corresponds to a valid
!          tracer name.
!   </ERROR>
if (TRACER_ARRAY(model,tracer_index) < 1 .or. TRACER_ARRAY(model,tracer_index) > num_tracer_fields) &
    call mpp_error(FATAL,'get_tracer_tlevels : invalid index')
!   <ERROR MSG="tracer field array not allocated" STATUS="FATAL">
!      The tracer array has not been allocated. This means that a
!          call to assign_tracer_field is absent in the code.
!   </ERROR>
if (.not. associated(tracers(TRACER_ARRAY(model,tracer_index))%field_tlevels)) &
    call mpp_error(FATAL,'get_tracer_tlevels : tracer field array not allocated')
data =>  tracers(TRACER_ARRAY(model,tracer_index))%field_tlevels

end function get_tracer_tlevels
!</FUNCTION>
!
!#######################################################################
!
!<FUNCTION NAME= "get_tracer_tendency">
!   <OVERVIEW>
!     A function to retrieve the tendency data.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Function to point to the 3D field associated with a tracer.
!   </DESCRIPTION>
!   <TEMPLATE>
!     array=get_tracer_tendency(model, tracer_index)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="tracer_index" TYPE="integer">
!     The tracer number within the component model.
!   </IN>
!   <OUT NAME="data"  TYPE="real, pointer" DIM="(:,:,:)">
!     The tracer tendency field is returned in this array.
!   </OUT>
function get_tracer_tendency(model, tracer_index) result (data)

integer, intent(in)  :: model
integer :: tracer_index
real, pointer        :: data(:,:,:)

integer :: n

!Convert local model index to tracer_manager index
!   <ERROR MSG="invalid index" STATUS="FATAL">
!     The index that has been passed to this routine is invalid.
!          Check the index that is being passed corresponds to a valid
!          tracer name.
!   </ERROR>
if (TRACER_ARRAY(model,tracer_index) < 1 .or. TRACER_ARRAY(model,tracer_index) > num_tracer_fields) &
    call mpp_error(FATAL,'get_tracer_tendency : invalid index')
!   <ERROR MSG="tracer tendency field array not allocated" STATUS="FATAL">
!      The tracer array has not been allocated. This means that a
!          call to assign_tracer_field is absent in the code.
!   </ERROR>
if (.not. associated(tracers(TRACER_ARRAY(model,tracer_index))%field_tendency)) &
    call mpp_error(FATAL,'get_tracer_tendency : tracer tendency field array not allocated')
data =>  tracers(TRACER_ARRAY(model,tracer_index))%field_tendency

end function get_tracer_tendency
!</FUNCTION>

!#######################################################################
!
! <SUBROUTINE NAME="get_tracer_names" >
!   <OVERVIEW>
!     Routine to find the names associated with a tracer number.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine can return the name, long name and units associated
!     with a tracer.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_tracer_names(model,n,name,longname, units)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number.
!   </IN>
!   <OUT NAME="name" TYPE="character" >
!     Field name associated with tracer number.
!   </OUT>
!   <OUT NAME="longname" TYPE="character, optional" >
!     The long name associated with tracer number.
!   </OUT>
!   <OUT NAME="units" TYPE="character, optional" >
!     The units associated with tracer number.
!   </OUT>

subroutine get_tracer_names(model,n,name,longname, units)

integer,          intent(in)  :: model, n
character (len=*),intent(out) :: name
character (len=*), intent(out), optional :: longname, units

if ( n == NO_TRACER ) then
    call mpp_error(NOTE,'get_tracer_names : Trying to get names for a non-existent tracer ')
endif  

!Convert local model index to tracer_manager index
!if (n < 1 .or. n > num_tracer_fields) &
if (TRACER_ARRAY(model,n) < 1 .or. TRACER_ARRAY(model,n) > num_tracer_fields) &
    call mpp_error(FATAL,'get_tracer_names : invalid tracer index for '//trim(name))

name=tracers(TRACER_ARRAY(model,n))%tracer_name
if (PRESENT(longname)) longname =tracers(TRACER_ARRAY(model,n))%tracer_longname
if (PRESENT(units)) units =tracers(TRACER_ARRAY(model,n))%tracer_units
!name=tracers(n)%tracer_name
!if (PRESENT(longname)) longname =tracers(n)%tracer_longname
!if (PRESENT(units)) units =tracers(n)%tracer_units

end subroutine get_tracer_names
!</SUBROUTINE>
!
!#######################################################################
!
! <SUBROUTINE NAME="get_family_name" >
!   <OVERVIEW>
!     Routine to return the family name for tracer n.
!   </OVERVIEW>
!   <DESCRIPTION>
!     You may wish to use this routine to retrieve the name of the family
!     that a tracer belongs to.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_family_name(model,n,name)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number that you want the family name for.
!   </IN>
!   <OUT NAME="name" TYPE="character" >
!     The family name.
!   </OUT>
subroutine get_family_name(model,n,name)

integer,          intent(in)  :: model, n
character (len=*),intent(out) :: name

!Convert local model index to tracer_manager index
name=tracers(TRACER_ARRAY(model,n))%tracer_family

end subroutine get_family_name
!</SUBROUTINE>
!
!#######################################################################
!
!<FUNCTION NAME= "check_if_prognostic">
!   <OVERVIEW>
!    Function to see if a tracer is prognostic or diagnostic.
!   </OVERVIEW>
!   <DESCRIPTION>
!    All tracers are assumed to be prognostic when read in from the field_table
!    However a tracer can be changed to a diagnostic tracer by adding the line
!    "tracer_type","diagnostic"
!    to the tracer description in field_table.
!   </DESCRIPTION>
!   <TEMPLATE>
!     logical =check_if_prognostic(model, n)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number that you want the family name for.
!   </IN>
!   <OUT NAME="check_if_prognostic" TYPE="logical">
!     A logical flag set TRUE if the tracer is 
!                        prognostic.
!   </OUT>
function check_if_prognostic(model, n)

integer, intent(in) :: model, n
logical             :: check_if_prognostic

!Convert local model index to tracer_manager index

check_if_prognostic = tracers(TRACER_ARRAY(model,n))%is_prognostic

end function check_if_prognostic
!</FUNCTION>
!
!#######################################################################
!
! <SUBROUTINE NAME="find_family_members" >
!   <OVERVIEW>
!     Subroutine to find which tracers are members of family family_name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Subroutine to find which tracers are members of family family_name.
! This will return a logical array where the array positions 
! corresponding to the tracer numbers for family members are set .TRUE.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call find_family_members(model, family_name,is_family_member)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="family_name" TYPE="character">
!      The family name of the members one is seeking.
!   </IN>
!   <OUT NAME="is_family_member" TYPE="logical" DIM ="(:)">
!      A logical array where the tracer number is used as 
!                     the index to signify which tracer is part of the family.
!                     i.e. If tracers 1, 3, and 7 are part of the same family
!                     then is_family_member(1), is_family_member(3), and 
!                     is_family_member(7) are set TRUE.
!   </OUT>
subroutine find_family_members(model, family_name,is_family_member)

integer, intent(in)         :: model
character(len=*),intent(in) :: family_name
logical , intent(out)       :: is_family_member(:)
integer ::n

is_family_member = .false.

if (size(is_family_member(:)) < num_tracer_fields) call mpp_error(FATAL,'find_family_members : array too short')

do n=1,num_tracer_fields
   if(trim(tracers(TRACER_ARRAY(model,n))%tracer_family) == trim(family_name) .and. &
      tracers(TRACER_ARRAY(model,n))%model == model ) then
      is_family_member(n)= .true.
   endif
enddo

end subroutine find_family_members
!</SUBROUTINE>

!
!#######################################################################
!
subroutine check_family_parameters(model, family_name)
!
! Subroutine to check that the advection and diffusion schemes and 
! controls of each tracer in a family is the same.
! INTENT IN
!  model       : The model that you are calling this subroutine from.
!  family_name : The name that has been given to the family of tracers of interest
!
character(len=*),intent(in) :: family_name
integer, intent(in)         :: model
logical            :: is_family_member(num_tracer_fields), flag
integer            :: nfam,n
character(len=128) :: name,advectvert, advecthoriz, diffusvert, diffushoriz

nfam = get_tracer_index(model, family_name)
flag=.true.

if (nfam > 0) then 
   if( query_method ('advection_scheme_vert',model,nfam,name)) advectvert = name
   if( query_method ('advection_scheme_horiz',model,nfam,name)) advecthoriz = name
   if( query_method ('diffusion_scheme_vert',model,nfam,name)) diffusvert = name
   if( query_method ('diffusion_scheme_horiz',model,nfam,name)) diffushoriz = name
   call find_family_members(model,family_name,is_family_member)
   do n=1,num_tracer_fields
      if(is_family_member(n)) then
         if( query_method ('advection_scheme_vert',model,n,name))  then
            if (name .ne. advectvert)  flag=.FALSE.
         endif
         if( query_method ('advection_scheme_horiz',model,n,name)) then
            if (name .ne. advecthoriz)  flag=.FALSE.
         endif
         if( query_method ('diffusion_scheme_vert',model,n,name))  then
            if (name .ne. diffusvert)  flag=.FALSE.
         endif
         if( query_method ('diffusion_scheme_horiz',model,n,name)) then
            if (name .ne. diffushoriz)  flag=.FALSE.
         endif
      endif
      if(.not. flag) call mpp_error     &
           (FATAL,'check_family_parameters : Family members do not have the same parameters for advection and diffusion')
   enddo
endif

end subroutine check_family_parameters
!
!#######################################################################
!
! <SUBROUTINE NAME="add_members_to_family" >
!   <OVERVIEW>
!     Routine to sum up the members of a family of tracers so that they may
! be advected and diffused as one tracer.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Routine to sum up the members of a family of tracers so that they may
! be advected and diffused as one tracer. This should only be used in
! conjunction with split_family_into_members and should be placed before
! the advection scheme is called. 
!   </DESCRIPTION>
!   <TEMPLATE>
!     call add_members_to_family(model,family_name, cur, prev, next)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number.
!   </IN>
!   <IN NAME="cur" TYPE="integer, optional">
!     Array index for the current time step. This is only of use 
!                with a three timestep model.
!   </IN>
!   <IN NAME="prev" TYPE="integer, optional">
!     Array index for the previous time step. This is only of use 
!                with a three timestep model.
!   </IN>
!   <IN NAME="next" TYPE="integer, optional">
!     Array index for the next time step. This is only of use 
!                with a three timestep model.
!   </IN>

!   <NOTE>
! This should be used with extreme caution. 
! Unless the family member distributions are similar to each other spatially, 
! advection as one tracer and subsequent splitting will result in a different
! result to advecting each tracer separately. The user should understand the 
! possible repercussions of this before using it.
!</NOTE>
subroutine add_members_to_family(model,family_name, cur, prev, next)

integer, intent(in)                     :: model
character(len=*), intent(in)            :: family_name
integer, intent(in), optional           :: cur, prev, next

!
logical :: is_family_member(num_tracer_fields)
integer :: nfam,n, min_indx, t0, tm1, tp1, siz(num_tracer_fields,3)

nfam=0
nfam = get_tracer_index(model,family_name)
!This is the family number within whichever model
!is calling this routine. So need to convert it to 
! the tracer_manager tracer number.
nfam = TRACER_ARRAY(model,nfam)
! If the tracer is not a family it doesn't need to be summed.
if (.not.tracers(nfam)%is_family) return


if (nfam > 0) then
   call find_family_members(model, family_name,is_family_member)    
   if (associated(tracers(nfam)%field).and.associated(tracers(nfam)%field_tendency)) then 
       tracers(nfam)%field = 0.0 ! initialize
       do n=1,num_tracer_fields
          if (is_family_member(n)) then
              if (.not.associated(tracers(n)%field)) &
                  call mpp_error(FATAL,'add_members_to_family : current tracer field not associated')
              siz(n,1) = size(tracers(n)%field,1)
              siz(n,2) = size(tracers(n)%field,2)
              siz(n,3) = size(tracers(n)%field,3)
              tracers(nfam)%field=tracers(nfam)%field+tracers(n)%field 
          endif
       enddo
! Now divide the members by the family total
! The members should then add up to 1.
       do n=1,num_tracer_fields
          if (is_family_member(n)) then
             if (.not.associated(tracers(n)%weight)) allocate(tracers(n)%weight(siz(n,1),siz(n,2),siz(n,3)))
              where(tracers(nfam)%field /= 0.0)
                  tracers(n)%weight = tracers(n)%field/tracers(nfam)%field
              elsewhere
                  tracers(n)%weight = 0.0
              end where
              tracers(n)%is_combined = .true.
          endif
       enddo
    else if (associated(tracers(nfam)%field_tlevels)) then
       if (.not.PRESENT(cur).or..not.PRESENT(prev).or..not.PRESENT(next)) then
           call mpp_error(FATAL,'add_members_to_family : need to specify time level indices to add family members')
       endif
       ! pointer array indices start at 1
       min_indx = min(cur,prev,next)
       t0 = cur-min_indx+1
       tm1 = prev-min_indx+1
       tp1 = next-min_indx+1
       tracers(nfam)%field_tlevels(:,:,:,:)=0.0 ! initialize
       do n=1,num_tracer_fields
          if (is_family_member(n)) then
             if (.not.associated(tracers(n)%field_tlevels)) &
                 call mpp_error(FATAL,'add_members_to_family : tracer time levels not associated')
             tracers(nfam)%field_tlevels(:,:,:,t0)=tracers(nfam)%field_tlevels(:,:,:,t0)+&
                  tracers(n)%field_tlevels(:,:,:,t0)
             tracers(nfam)%field_tlevels(:,:,:,tm1)=tracers(nfam)%field_tlevels(:,:,:,tm1)+&
                  tracers(n)%field_tlevels(:,:,:,tm1)
             siz(n,1) = size(tracers(n)%field_tlevels,1)
             siz(n,2) = size(tracers(n)%field_tlevels,2)
             siz(n,3) = size(tracers(n)%field_tlevels,3)
          endif
       enddo
! Now divide the members by the family total
! The members should then add up to 1.
       do n=1,num_tracer_fields
          if (is_family_member(n)) then
              if (.not.associated(tracers(n)%weight)) allocate(tracers(n)%weight(siz(n,1),siz(n,2),siz(n,3)))
              where(tracers(nfam)%field_tlevels(:,:,:,t0) /= 0.0)  
                  tracers(n)%weight(:,:,:) = tracers(n)%field_tlevels(:,:,:,t0)/&
                                                tracers(nfam)%field_tlevels(:,:,:,t0)
              elsewhere
                  tracers(n)%weight = 0.0
              end where
              tracers(n)%is_combined = .true.
          endif
       enddo
    else
       call mpp_error(FATAL,'add_members_to_family : need to associate tracer pointers to use families')
    endif
 endif

 return

end subroutine add_members_to_family
!</SUBROUTINE>

!
!#######################################################################
!
! Subroutine that sets the present value of the member of a tracer 
! family according to the fraction of the family that it was in the 
! previous step.
!
! <SUBROUTINE NAME="split_family_into_members" >
!   <OVERVIEW>
!      Subroutine that sets the present value of the member of a tracer 
! family according to the fraction of the family that it was in the 
! previous step.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Subroutine that sets the present value of the member of a tracer 
! family according to the fraction of the family that it was in the 
! previous step.
!
! This splits the transported family into the constituent members. This
! should only be used in conjunction with <I>add_members_to_family</I> and
! should be placed after the advection scheme is called.
!
!   </DESCRIPTION>
!   <TEMPLATE>
!     call split_family_into_members(model,family_name,cur,prev,next)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="family_name" TYPE="character">
!     The name of the family of tracers that you would 
!                like to split up.
!   </IN>
!   <IN NAME="cur" TYPE="integer, optional">
!     Array index for the current time step. This is only of use 
!                with a three timestep model.
!   </IN>
!   <IN NAME="prev" TYPE="integer, optional">
!     Array index for the previous time step. This is only of use 
!                with a three timestep model.
!   </IN>
!   <IN NAME="next" TYPE="integer, optional">
!     Array index for the next time step. This is only of use 
!                with a three timestep model.
!   </IN>

!   <NOTE>
! This should be used with extreme caution. 
! Unless the family member distributions are similar to each other spatially, 
! advection as one tracer and subsequent splitting will result in a different
! result to advecting each tracer separately. The user should understand the 
! possible repercussions of this before using it.
!</NOTE>
subroutine split_family_into_members(model,family_name,cur,prev,next)

integer, intent(in) :: model
character(len=*),intent(in) :: family_name
integer, intent(in), optional :: cur, prev, next

logical :: is_family_member(num_tracer_fields)
integer :: nfam,n, min_indx, t0, tm1, tp1

! need to make sure tracers have been combined already

nfam=0
nfam = get_tracer_index(model,family_name)
! This is the family tracer number within whichever model
! is calling this routine. So need to convert it to 
! the tracer_manager tracer number.
nfam = TRACER_ARRAY(model,nfam)

! If the family name is not the same as the tracer name then the
! tracer is not a family tracer (It may be a member of a family)
if (tracers(nfam)%tracer_name .ne. tracers(nfam)%tracer_family) return

if (mpp_pe() ==mpp_root_pe() ) write(*,*) 'split',model, trim(family_name),nfam,TRACER_ARRAY(model,nfam)

if (nfam > 0) then 
   call find_family_members(model,family_name,is_family_member)
   if (associated(tracers(nfam)%field).and.associated(tracers(nfam)%field_tendency)) then 
      do n=1,num_tracer_fields
         if (is_family_member(n)) then           
            if (.not. tracers(n)%is_combined) &
                call mpp_error(FATAL,'split_family_into_members : call to split family into members when fields are not combined')
            tracers(n)%field = tracers(n)%weight*tracers(nfam)%field 
            tracers(n)%is_combined = .false.
         endif
      enddo
   else if (associated(tracers(nfam)%field_tlevels)) then
      if (.not.PRESENT(cur).or..not.PRESENT(prev).or..not.PRESENT(next)) then
         call mpp_error(FATAL,'split_family_into_members : need to specify time level indices to split family members')
      endif
      ! pointer array indices start at 1
      min_indx = min(cur,prev,next)
      t0 = cur-min_indx+1
      tm1 = prev-min_indx+1
      tp1 = next-min_indx+1
      do n=1,num_tracer_fields
         if (is_family_member(n)) then
            if (.not. tracers(n)%is_combined) &
                call mpp_error(FATAL,'split_family_into_members : call to split family into members when fields are not combined')
            tracers(n)%field_tlevels(:,:,:,tp1) = tracers(nfam)%field_tlevels(:,:,:,tp1)*tracers(n)%weight(:,:,:)
            tracers(n)%is_combined = .false.
         endif
      enddo
   else
      call mpp_error(FATAL,'split_family_into_members : need to associate tracer pointers to use families')
   endif
endif


end subroutine split_family_into_members
!</SUBROUTINE>

!
!#######################################################################
!
! <SUBROUTINE NAME="set_tracer_profile" >
!   <OVERVIEW>
!     Subroutine to set the tracer field to the wanted profile.
!   </OVERVIEW>
!   <DESCRIPTION>
!     If the profile type is 'fixed' then the tracer field values are set 
! equal to the surface value.
! If the profile type is 'profile' then the top/bottom of model and
! surface values are read and an exponential profile is calculated,
! with the profile being dependent on the number of levels in the
! component model. This should be called from the part of the dynamical
! core where tracer restarts are called in the event that a tracer
! restart file does not exist.
!
!  This can be activated by adding a method to the field_table
! e.g.
!  "profile_type","fixed","surface_value = 1e-12"
!  would return values of surf_value = 1e-12 and a multiplier of 1.0
!  One can use these to initialize the entire field with a value of 1e-12.
!
!  "profile_type","profile","surface_value = 1e-12, top_value = 1e-15"
!   In a 15 layer model this would return values of surf_value = 1e-12 and 
!   multiplier = 0.6309573 i.e 1e-15 = 1e-12*(0.6309573^15)
!   In this case the model should be MODEL_ATMOS as you have a "top" value.
!
!   If you wish to initialize the ocean model, one can use bottom_value instead
!   of top_value.

!   </DESCRIPTION>
!   <TEMPLATE>
!     call set_tracer_profile(model, n, tracer)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number.
!   </IN>
!   <INOUT NAME="tracer_array" TYPE="real">
!     The initialized tracer array.
!   </INOUT>
subroutine set_tracer_profile(model, n, tracer)

integer,  intent(in)  :: model, n
   real, intent(inout), dimension(:,:,:) :: tracer

real    :: surf_value, multiplier
integer :: numlevels, k, m
real    :: top_value, bottom_value
character(len=80) :: scheme, control,profile_type, name
integer :: flag

!default values
profile_type  = 'Fixed'
surf_value = 0.0E+00
top_value  = surf_value
bottom_value = surf_value
multiplier = 1.0

tracer = surf_value

if ( query_method ( 'profile_type',model,n,scheme,control)) then
!Change the tracer_number to the tracer_manager version

  if(lowercase(trim(scheme(1:5))).eq.'fixed') then
    profile_type                   = 'Fixed'
    flag =parse(control,'surface_value',surf_value)
    multiplier = 1.0
    tracer = surf_value
    call get_tracer_names(model, n, name)
    call have_initialized_tracer(model, name)
  endif

  if(lowercase(trim(scheme(1:7))).eq.'profile') then
    profile_type                   = 'Profile'
    flag=parse(control,'surface_value',surf_value)
    if (surf_value .eq. 0.0) &
      call mpp_error(FATAL,'set_tracer_profile : Cannot have a zero surface value for an exponential profile. Tracer '&
                           //tracers(TRACER_ARRAY(model,n))%tracer_name//" "//control//" "//scheme)
    select case (tracers(TRACER_ARRAY(model,n))%model)
      case (MODEL_ATMOS)
        flag=parse(control,'top_value',top_value)
        if(mpp_pe()==mpp_root_pe() .and. flag == 0) &
           call mpp_error(NOTE,'set_tracer_profile : Parameter top_value needs to be defined for the tracer profile.')
      case (MODEL_OCEAN)
        flag =parse(control,'bottom_value',bottom_value)
        if(mpp_pe() == mpp_root_pe() .and. flag == 0) &
           call mpp_error(NOTE,'set_tracer_profile : Parameter bottom_value needs to be defined for the tracer profile.')
      case default
    end select

! If profile type is profile then set the surface value to the input
! value and calculate the vertical multiplier.
! 
! Assume an exponential decay/increase from the surface to the top level
!  C = C0 exp ( -multiplier* level_number)
!  => multiplier = exp [ ln(Ctop/Csurf)/number_of_levels]
!
numlevels = size(tracer,3) -1
    if (associated(tracers(TRACER_ARRAY(model,n))%field)) numlevels = size(tracers(TRACER_ARRAY(model,n))%field,3) -1
    if (associated(tracers(TRACER_ARRAY(model,n))%field_tlevels)) numlevels = size(tracers(TRACER_ARRAY(model,n))%field_tlevels,3)-1
    select case (tracers(TRACER_ARRAY(model,n))%model)
      case (MODEL_ATMOS)
        multiplier = exp( log (top_value/surf_value) /numlevels)
        tracer(:,:,1) = surf_value
        do k = 2, size(tracer,3)
          tracer(:,:,k) = tracer(:,:,k-1) * multiplier
        enddo
        call get_tracer_names(model, n, name)
        call have_initialized_tracer(model, name)
      case (MODEL_OCEAN)
        multiplier = exp( log (bottom_value/surf_value) /numlevels)
        tracer(:,:,size(tracer,3)) = surf_value
        do k = size(tracer,3) - 1, 1, -1
          tracer(:,:,k) = tracer(:,:,k+1) * multiplier
        enddo
        call get_tracer_names(model, n, name)
        call have_initialized_tracer(model, name)
      case default
    end select
  endif !scheme.eq.profile

  if (mpp_pe() == mpp_root_pe() ) write(*,700) 'Tracer ',trim(tracers(TRACER_ARRAY(model,n))%tracer_name),    &
                            ' initialized with surface value of ',surf_value, &
                            ' and vertical multiplier of ',multiplier
  700 FORMAT (3A,E12.6,A,F10.6)

endif ! end of query scheme

end subroutine set_tracer_profile
!</SUBROUTINE>

!
!#######################################################################
!
! <FUNCTION NAME="query_method" >
!   <OVERVIEW>
!     A function to query the "methods" associated with each tracer.
!   </OVERVIEW>
!   <DESCRIPTION>
!     A function to query the "methods" associated with each tracer. The
!  "methods" are the parameters of the component model that can be
!  adjusted by user by placing formatted strings, associated with a
!  particular tracer, within the field table.
!  These methods can control the advection, wet deposition, dry
!  deposition or initial profile of the tracer in question. Any
!  parametrization can use this function as long as a routine for parsing
!  the name and control strings are provided by that routine.
!   </DESCRIPTION>
!   <TEMPLATE>
!     logical =query_method  (method_type, model, n, name, control)
!   </TEMPLATE>

!   <IN NAME="method_type" TYPE="character">
!     The method that is being requested.
!   </IN>
!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number that you want the family name for.
!   </IN>
!   <OUT NAME="name" TYPE="character">
!     A string containing the modified name to be used with
!     method_type. i.e. "2nd_order" might be the default for 
!     advection. One could use "4th_order" here to modify 
!     that behaviour.
!   </OUT>
!   <OUT NAME="control" TYPE="character, optional">
!     A string containing the modified parameters that are 
!     associated with the method_type and name.
!   </OUT>
!   <OUT NAME="query_method" TYPE="logical">
!      A flag to show whether method_type exists with regard to
!      tracer n. If method_type is not present then one must
!      have default values.
!   </OUT>

!<NOTE>
!  At present the tracer manager module allows the initialization of a tracer
!  profile if a restart does not exist for that tracer. 
!  Options for this routine are as follows
!
!  Tracer profile setup
!  ==================================================================
!  |method_type  |method_name  |method_control                      |
!  ==================================================================
!  |profile_type |fixed        |surface_value = X                   |
!  |profile_type |profile      |surface_value = X, top_value = Y    |(atmosphere)
!  |profile_type |profile      |surface_value = X, bottom_value = Y |(ocean)
!  ==================================================================
!
!</NOTE>
function query_method  (method_type, model, n, name, control)
!
!  A function to query the schemes associated with each tracer. 
!  
!  INTENT IN
!   method_type  : The method that is being requested.
!   model        : The model that you are calling this function from.
!   n            : The tracer number.
!  INTENT OUT
!   name         : A string containing the modified name to be used with
!                  method_type. i.e. "2nd_order" might be the default for 
!                  advection. One could use "4th_order" here to modify 
!                  that behaviour.
!   control      : A string containing the modified parameters that are 
!                  associated with the method_type and name.
!   query_method : A flag to show whether method_type exists with regard 
!                  to tracer n. If method_type is not present then one
!                  must have default values.

character(len=*), intent(in)            :: method_type
integer         , intent(in)            :: model, n
character(len=*), intent(out)           :: name
character(len=*), intent(out), optional :: control
logical                                 :: query_method

integer :: m, n1, nmethods
character(len=256) :: list_name, control_tr
character(len=20)  :: nindex
!Convert the local model tracer number to the tracer_manager version.

if ( n == NO_TRACER ) then
    call mpp_error(NOTE,'query_method : Trying to query method for a non-existent tracer ')
    call mpp_error(FATAL,'query_method : Query method is : '//method_type )
    
endif  
if (n < 1 .or. n > num_tracer_fields) then 
   write(nindex, '(I20)') n
   if (mpp_pe() .eq. mpp_root_pe()) &
     call mpp_error(NOTE,'query_method : Invalid tracer index '//trim(nindex))
   n1 = NO_TRACER
   return
endif

name=" "
query_method = .false.
n1 = TRACER_ARRAY(model,n)

select case(model)
 case (MODEL_ATMOS)
  list_name = "/atmos_mod/tracer/"//trim(tracers(n1)%tracer_name)//"/"//trim(method_type)
 case (MODEL_OCEAN)
  list_name = "/ocean_mod/tracer/"//trim(tracers(n1)%tracer_name)//"/"//trim(method_type)
 case (MODEL_ICE  )
  list_name = "/ice_mod/tracer/"//trim(tracers(n1)%tracer_name)//"/"//trim(method_type)
 case (MODEL_LAND )
  list_name = "/land_mod/tracer/"//trim(tracers(n1)%tracer_name)//"/"//trim(method_type)
 case default
  list_name = "/default/tracer/"//trim(tracers(n1)%tracer_name)//"/"//trim(method_type)
end select

  name = ''
  control_tr = ''
  query_method = fm_query_method(list_name, name, control_tr)

if ( present(control)) control = trim(control_tr)

end function query_method
!</FUNCTION>

! <FUNCTION NAME="query_combined" >
!   <OVERVIEW>
!     A function to query whether families of tracers have been combined already.
!   </OVERVIEW>
!   <DESCRIPTION>
!     A function to query whether families of tracers have been combined already.
!  This function should only be used in conjunction with add_members_to_family 
!  and split_family_into_members.
!   </DESCRIPTION>
!   <TEMPLATE>
!     logical =query_combined (model, index)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="index" TYPE="integer">
!     Tracer number.
!   </IN>
!   <OUT NAME="query_combined" TYPE="logical" >
!     A flag to show whether the tracer family has been combined.
!   </OUT>
function query_combined (model, index) 

integer,intent(in) :: model, index
logical :: query_combined

if (index < 0 .or. index > num_tracer_fields) call mpp_error(FATAL,'query_combined : invalid tracer index')

query_combined = .false.
if (tracers(TRACER_ARRAY(model,index))%is_combined) query_combined = .true.

return
end function query_combined
!</FUNCTION>

!<SUBROUTINE NAME="set_tracer_atts">
!   <OVERVIEW>
!     A subroutine to allow the user set the tracer longname and units from the 
!     tracer initialization routine.
!   </OVERVIEW>
!   <DESCRIPTION>
!     A function to allow the user set the tracer longname and units from the 
!     tracer initialization routine. It seems sensible that the user who is 
!     coding the tracer code will know what units they are working in and it 
!     is probably safer to set the value in the tracer code rather than in 
!     the field table.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call set_tracer_atts(model, name, longname, units)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="name" TYPE="character">
!     Tracer name.
!   </IN>
!   <OUT NAME="longname" TYPE="character, optional">
!     A string describing the longname of the tracer for output to NetCDF files
!   </OUT>
!   <OUT NAME="units" TYPE="character, optional">
!     A string describing the units of the tracer for output to NetCDF files
!   </OUT>
subroutine set_tracer_atts(model, name, longname, units)

integer, intent(in)                    :: model
character(len=*), intent(in)           :: name
character(len=*), intent(in), optional :: longname, units

integer :: n, index
logical :: success
character(len=128) :: list_name


n = get_tracer_index(model,name)

if ( n .ne. NO_TRACER ) then
    tracers(TRACER_ARRAY(model,n))%tracer_units   = units
    tracers(TRACER_ARRAY(model,n))%tracer_longname = longname
  select case(model)
    case(MODEL_ATMOS) 
      list_name = "/atmos_mod/tracer/"//trim(name)
    case(MODEL_OCEAN) 
      list_name = "/ocean_mod/tracer/"//trim(name)
    case(MODEL_LAND) 
      list_name = "/land_mod/tracer/"//trim(name)
    case(MODEL_ICE) 
      list_name = "/ice_mod/tracer/"//trim(name)
    case DEFAULT 
      list_name = "/"//trim(name)
  end select      

! Method_type is a list, method_name is a name of a parameter and method_control has the value.
!    list_name = trim(list_name)//"/longname"
  if ( fm_exists(list_name)) then
    success = fm_change_list(list_name)
    if ( present(longname) ) then
      if ( longname .ne. "" ) index = fm_new_value('longname',longname)
    endif
    if ( present(units) ) then
      if (units .ne. "" ) index = fm_new_value('units',units)
    endif
  endif  

    
else
    call mpp_error(NOTE,'set_tracer_atts : Trying to set longname and/or units for non-existent tracer : '//trim(name))
endif

end subroutine set_tracer_atts
!</SUBROUTINE>


!<SUBROUTINE NAME="set_tracer_method">
!   <OVERVIEW> 
!      A subroutine to allow the user to set some tracer specific methods.
!   </OVERVIEW>
!   <DESCRIPTION>
!      A subroutine to allow the user to set methods for a specific tracer. 
!   </DESCRIPTION>
!   <TEMPLATE>
!     call set_tracer_method(model, name, method_type, method_name, method_control)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="name" TYPE="character">
!     Tracer name.
!   </IN>
!   <IN NAME="method_type" TYPE="character">
!     The type of the method to be set.
!   </IN>
!   <IN NAME="method_name" TYPE="character">
!     The name of the method to be set.
!   </IN>
!   <IN NAME="method_control" TYPE="character">
!     The control parameters of the method to be set.
!   </IN>
     
subroutine set_tracer_method(model, name, method_type, method_name, method_control)


integer, intent(in)                    :: model
character(len=*), intent(in)           :: name
character(len=*), intent(in)           :: method_type
character(len=*), intent(in)           :: method_name
character(len=*), intent(in)           :: method_control

integer :: n, num_method, index
logical :: success
character(len=128) :: list_name

n = get_tracer_index(model,name)

if ( n .ne. NO_TRACER ) then
  tracers(n)%num_methods = tracers(n)%num_methods + 1
  num_method = tracers(n)%num_methods

  select case(model)
    case(MODEL_ATMOS)
      list_name = "/atmos_mod/tracer/"//trim(name)
    case(MODEL_OCEAN)
      list_name = "/ocean_mod/tracer/"//trim(name)
    case(MODEL_LAND)
      list_name = "/land_mod/tracer/"//trim(name)
    case(MODEL_ICE)
      list_name = "/ice_mod/tracer/"//trim(name)
    case DEFAULT
      list_name = "/"//trim(name)
  end select      

  if ( method_control .ne. "" ) then
! Method_type is a list, method_name is a name of a parameter and method_control has the value.
    list_name = trim(list_name)//"/"//trim(method_type)
    if ( fm_exists(list_name)) then
      success = fm_change_list(list_name)
      index = fm_new_value(method_type,method_control)
    endif
  else
    call mpp_error(NOTE,'set_tracer_method : Trying to set a method for non-existent tracer : '//trim(name))
  endif
endif

end subroutine set_tracer_method
!</SUBROUTINE>

end module tracer_manager_mod

#ifdef test_tracer_manager

program test

use field_manager_mod, only : MODEL_ATMOS, &
                              MODEL_OCEAN, &
                              field_manager_init, &
                              fm_dump_list, &
                              fm_query_method, &
                              fm_exists, &
                              fm_new_value, &
                              fm_get_value, &
                              fm_change_list
use tracer_manager_mod
use mpp_mod, only : mpp_npes, mpp_pe, mpp_root_pe, FATAL, mpp_error, mpp_exit
use mpp_domains_mod, only : domain2d, mpp_define_domains, cyclic_global_domain, &
                            mpp_get_data_domain, mpp_get_compute_domain, mpp_update_domains, &
                            mpp_get_global_domain
use mpp_io_mod, only : mpp_io_init, mpp_io_exit
use diag_manager_mod, only : register_diag_field, send_data, diag_axis_init, diag_manager_init, diag_manager_end
use time_manager_mod

implicit none

integer, parameter :: nx = 100, ny = 100, nz = 10, ntsteps = 20, delt = 1, two_delt = 2*delt
real, parameter :: advect_speed = .23
integer :: num_tracer_prog_atmos, num_tracer_diag_atmos, num_tracer_fam_atmos, num_tracers_atmos
integer :: num_tracer_prog_ocean, num_tracer_diag_ocean, num_tracer_fam_ocean, num_tracers_ocean
integer :: i, axes(3)
integer, allocatable, dimension(:) :: atmos_prog_ind, atmos_diag_ind, atmos_fam_ind, atmos_prog_tracer_diagnostic_id
integer, allocatable, dimension(:) :: atmos_diag_tracer_diagnostic_id
integer, allocatable, dimension(:) :: ocean_prog_ind, ocean_diag_ind, ocean_fam_ind

real, allocatable, dimension(:,:,:,:,:), target :: atmos_tracers, ocean_tracers, atmos_fam_tracers
real, allocatable, dimension(:,:,:,:), target :: atmos_diag_tracers, ocean_diag_tracers
real, allocatable, dimension(:,:,:), target :: atmos_tendency, ocean_tendency  
real, allocatable, dimension(:,:,:,:) :: atmos_vel, ocean_vel

real ::  missing_value = -1.e10
real, dimension(max(nx,ny,nz)) :: data
logical :: result, flag, success

type(time_type) :: model_time
type(domain2d) :: domain ! just using a single domain type for all tracers here since
                         ! allocating storage in a single array.  This is inefficient in the
                         ! case where different tracers use different order advection schemes
                         ! need to consider extending mpp_update_domains to only update a portion
                         ! of the halo OR using a derived type for storing tracer arrays with 
                         ! a domain2d component (this seems unworkable because diag_axis_init would
                         ! need different calls for different domains ...)
   
character(len=128) :: name, longname, units, scheme_name, control, family_name, list_name
character(len=512) :: meth_name, meth_control
integer :: halo, ndivs, isc, iec, jsc, jec, isd, ied, jsd, jed, tau, taum1, taup1, tmp, n, m, mm
integer, dimension(MAX_TRACER_FIELDS) :: atmos_order, ocean_order
integer :: co2_index, color_index, temp_index,k,nfields, index, jsphum
real :: surf_value,multiplier, param

call mpp_io_init
call set_calendar_type(JULIAN)
call diag_manager_init

taum1=1
tau=2
taup1=3

model_time = set_date(1982,1,1,0,0,0)

!Register the tracers for the atmospheric model.
call register_tracers(MODEL_ATMOS, num_tracers_atmos, num_tracer_prog_atmos, num_tracer_diag_atmos, &
                     num_tracer_fam_atmos)

write(*,*) "The number of tracers in the atmosphere AFTER THE CALL TO register_tracers are ", &
            num_tracers_atmos, num_tracer_prog_atmos, num_tracer_diag_atmos, num_tracer_fam_atmos                     
!Try to register the tracers for the atmospheric model a second time.
call register_tracers(MODEL_ATMOS, num_tracers_atmos, num_tracer_prog_atmos, num_tracer_diag_atmos, &
                     num_tracer_fam_atmos)
! This should result in a note that register_tracers has already been called for the atmospheric model.
! Silly exmaple of how to use get_number_tracers. 
! As the ocean tracers have not been registered by the tracer manager these will all be zero.
call get_number_tracers(MODEL_OCEAN, num_tracers_ocean, num_tracer_prog_ocean, num_tracer_diag_ocean, &
                     num_tracer_fam_ocean)
write(*,*) "The number of tracers in the ocean BEFORE THE CALL TO register_tracers are ", &
            num_tracers_ocean, num_tracer_prog_ocean, num_tracer_diag_ocean, num_tracer_fam_ocean                     

!Register the tracers for the oceanic model.
call register_tracers(MODEL_OCEAN, num_tracers_ocean, num_tracer_prog_ocean, num_tracer_diag_ocean, &
                     num_tracer_fam_ocean)

! The ocean tracers have now been registered by the tracer manager so these will be non-zero if there
! are ocean tracers in the field table.
write(*,*) "The number of tracers in the ocean AFTER THE CALL TO register_tracers are ", &
            num_tracers_ocean, num_tracer_prog_ocean, num_tracer_diag_ocean, num_tracer_fam_ocean                     

write(*,*) 'Using fm_dump_list("/",.true.) to provide full listing of fields'
result = fm_dump_list("/",.true.)

!Set the longname and units for the relative humidity tracer

write(*,*) 'Changing longname and units of the radon_contr tracer to "control_radon" and "g/g" '
write(*,*) 'Using call set_tracer_atts(MODEL_ATMOS,"radon_contr","control_radon","g/g")'
call set_tracer_atts(MODEL_ATMOS,"radon_contr","control_radon","g/g")

write(*,*) 'Using fm_dump_list("/atmos_mod/tracer/radon_contr",.true.) to provide listing of radon_contr field'
result = fm_dump_list("/atmos_mod/tracer/radon_contr",.true.)

!write(*,*) 'Using fm_query_method("/atmos_mod/tracer/radon",meth_name,meth_control) to query &
!           &for methods of the radon field.'
result = fm_query_method("/atmos_mod/TRACER/tsurf",meth_name,meth_control)
write(*,*) "meth_name = ",trim(meth_name)
write(*,*) "meth_control = ",trim(meth_control)

     jsphum = get_tracer_index ( model_atmos, 'sphum' )
write(*,*) "tracer index for sphum is  ",jsphum
     call get_tracer_names ( model_atmos, jsphum, name, longname, units )

write(*,*) 'name associated with ',jsphum,' is ',trim(name)

do i=1,max(nx,ny,nz)
   data(i) = float(i)
enddo

ndivs = mpp_npes()


! Allocate space for the various tracers
call alloc_tracers

! atmos initialization

do m=1,num_tracer_prog_atmos
       call get_tracer_names(MODEL_ATMOS,m,name,longname,units)
index= get_tracer_index(MODEL_ATMOS,name)
 write(*,*) trim(name)," ",trim(longname)," ", trim(units),index,atmos_prog_ind(m)

if(query_method ( 'emissions',MODEL_ATMOS,index,longname,units)) then
 write(*,*) trim(longname), " ", trim(units)

endif
enddo
do m=1,num_tracer_prog_atmos
! Find the name, longname and units of a tracer given a model and index.
   call get_tracer_names(MODEL_ATMOS, atmos_prog_ind(m), name, longname, units)
! Provide a flag that the tracer requires initialization.
   call tracer_requires_init(MODEL_ATMOS, name)
! Now query to see whether the tracer has been initialized.
! A .true. response means the tracer requires initialization.
! In this example we have not initialized the tracer yet.
write(*,*) 'Using query_tracer_init(MODEL_ATMOS, '//trim(name)//') to find out if tracer needs initialization'

   if ( query_tracer_init(MODEL_ATMOS, name) ) then
     write(*,*) 'Tracer '//trim(name)//' needs initialization'
! Set a gaussian peak for each tracer as an initial condition.
     call init_tracer(atmos_tracers(isc:iec,jsc:jec,:,m,tau),float(m), domain)
! Now tell the tracer manager that the tracer has been initialized.
     call have_initialized_tracer(MODEL_ATMOS, name)
! Query again. The answer here should always be false as the routine above will
! have told the tracer_manager that it has been initialized.
     if ( query_tracer_init(MODEL_ATMOS, name) ) &
       write(*,*) 'Tracer '//trim(name)//' STILL needs initialization'
   endif

! If a 'profile_type' method is defined in the field table then use it here.
   if(query_method ( 'profile_type',MODEL_ATMOS,atmos_prog_ind(m),longname,units)) then
      call set_tracer_profile(MODEL_ATMOS,atmos_prog_ind(m),atmos_tracers(:,:,:,m,tau))
   endif

! Special case for radon
!   if ( name == 'radon' ) call radon_init( atmos_tracers(:,:,:,m,tau) )

   call mpp_update_domains(atmos_tracers(:,:,:,m,tau), domain)
!Set the previous (tau - 1) timestep value equal to the present timestep value 
   atmos_tracers(:,:,:,m,taum1) = atmos_tracers(:,:,:,m,tau)
enddo

do m=1,num_tracer_diag_atmos
! Set a gaussian peak for each tracer as an initial condition.
   call init_tracer(atmos_diag_tracers(isc:iec,jsc:jec,:,m),float(m), domain)
! If a 'profile_type' method is defined in the field table then use it here.
   if(query_method ( 'profile_type',MODEL_ATMOS,atmos_diag_ind(m),longname,units)) then
      call set_tracer_profile(MODEL_ATMOS,atmos_diag_ind(m),atmos_diag_tracers(:,:,:,m))
   endif
   call mpp_update_domains(atmos_diag_tracers(:,:,:,m), domain)
enddo




! Ocean initialization
do m=1,num_tracer_prog_ocean
  call get_tracer_names(MODEL_OCEAN,ocean_prog_ind(m),name)
  list_name = "/ocean_mod/tracer/"//trim(name)
  select case(name)
    case( 'age_exp1')
      if ( fm_exists(list_name)) then
        success = fm_change_list(list_name)
        success =  fm_get_value(trim(list_name)//'/nlat',param)
        write(*,*) 'Changing nlat limits from ',param,' to 45.0 for ',list_name
        index = fm_new_value('nlat',45.0 )
        success =  fm_get_value(trim(list_name)//'/slat',param)
        write(*,*) 'Changing slat limits from ',param,' to -45.0 for ',list_name
        index = fm_new_value('slat',-45.0)
      endif  
     
    case( 'age_conv')
      write(*,*) 'Changing latitude limits for ',list_name
      if ( fm_exists(list_name)) then
        success = fm_change_list(list_name)
        write(*,*) "Using index = fm_new_value('nlat',30.0) to change value of nlat (still real)"
        index = fm_new_value('nlat',30.0)

! slat was originally real, but it can be changed to integer by using an integer argument AND the optional append = .true. argument 
        write(*,*) "Using index = fm_new_value('slat',-30, append = .true.) to change value of slat (now integer)"
        index = fm_new_value('slat',-30, append = .true.)
      endif  
     
    case( 'age_5')
!If the age_5 tracer exists then change the scale_in factor
      if ( fm_exists(list_name)) then
        success = fm_change_list(list_name)
        write(*,*) 'Changing scale_in factor for ',list_name
        index = fm_new_value('scale_in',0.5)
      endif  

    case( 'biotic1')
! If the biotic1 tracer exists, we will add a parameter for po4 uptake
      if ( fm_exists(list_name)) then
        success = fm_change_list(list_name)
        write(*,*) 'Adding a parameter to the list for ',list_name
        index = fm_new_value('po4_uptake',1e-5, create=.true., append = .true.)
      endif  

    case default
    
  end select
  
enddo

success = fm_dump_list('/ocean_mod/tracer', .true.)

! atmos time loop

!atmos_tendency = 0.0

do n=1,ntsteps
!write(*,*) 'Step ',n
   do m=1,num_tracer_prog_atmos

!This next line is more for convenience to check by name which tracer it is you are about to advect 
     call get_tracer_names(MODEL_ATMOS,atmos_prog_ind(m),name)
     call lateral_advection(atmos_order(m), atmos_tendency,tracer_now = atmos_tracers(:,:,:,m,tau),&
          vel_now = atmos_vel)
!          vel_now = atmos_vel(:,:,:,:))
     atmos_tracers(isc:iec,jsc:jec,:,m,taup1) = atmos_tracers(isc:iec,jsc:jec,:,m,taum1) - &
          atmos_tendency(isc:iec,jsc:jec,:)*two_delt
     call mpp_update_domains(atmos_tracers(:,:,:,m,taup1), domain)
     if (atmos_prog_tracer_diagnostic_id(m) > 0) result =  send_data(atmos_prog_tracer_diagnostic_id(m),&
                                   atmos_tracers(isc:iec,jsc:jec,:,m,tau),model_time)
   enddo
   do m=1,num_tracer_diag_atmos
      if (m.eq.color_index .and. temp_index .ne. NO_TRACER) &
         call calculate_color(atmos_tracers(:,:,:,temp_index,tau),atmos_diag_tracers(:,:,:,m))
      if (atmos_diag_tracer_diagnostic_id(m) > 0) &
           result = send_data(atmos_diag_tracer_diagnostic_id(m),&
          atmos_diag_tracers(isc:iec,jsc:jec,:,m),model_time)
   enddo

   do m=1,num_tracer_prog_ocean

   enddo

   tmp=taup1
   taup1=taum1
   taum1 = tau
   tau = tmp
   if (n.eq. ntsteps) call diag_manager_end(model_time)
   model_time = model_time + set_time(delt,0)      

enddo


call mpp_io_exit
call mpp_exit


contains 


subroutine alloc_tracers

if (num_tracer_prog_atmos > 0) then
! Allocate space for the atmospheric prognostic tracers
   allocate(atmos_prog_ind(num_tracer_prog_atmos))
   allocate(atmos_prog_tracer_diagnostic_id(num_tracer_prog_atmos))
! Get an array of indices for the atmospheric prognostic (prog_ind) tracers
   call get_tracer_indices(MODEL_ATMOS,prog_ind=atmos_prog_ind)

! Using the indices array retrieved above find the tracer index for co2 and temp
   co2_index = get_tracer_index(MODEL_ATMOS,'co2',atmos_prog_ind)
   temp_index = get_tracer_index(MODEL_ATMOS,'temp',atmos_prog_ind)
   
! NO_TRACER is returned if co2 does not exist so check versus NO_TRACER to test whether it exists.
   if (co2_index /= NO_TRACER .and. mpp_pe() == mpp_root_pe()) write(*,'(a)') 'co2 exists in tracer table '

   halo=1 ! default halo size for 2nd order advection
   do n=1,num_tracer_prog_atmos
     atmos_order(n)=2! Default advection scheme
     flag = query_method ('advection_scheme_horiz',MODEL_ATMOS,atmos_prog_ind(n),scheme_name,control)
     if(flag) then ! There is an advection scheme in the table so use that instead.
       select case (trim(scheme_name))
         case ('2nd_order') 
           atmos_order(n) = 2
         case ('4th_order')
           atmos_order(n) = 4
           halo = 2
         case default
           call mpp_error(FATAL,'invalid tracer advection scheme')
       end select
     endif
   enddo
   call mpp_define_domains( (/1,nx,1,ny/), (/1,ndivs/), domain, xflags = CYCLIC_GLOBAL_DOMAIN, &
                            yflags = CYCLIC_GLOBAL_DOMAIN, xhalo = halo, yhalo = halo)
   axes(1) = diag_axis_init('x',data(1:nx),units='meters',cart_name='x',domain2=domain)
   axes(2) = diag_axis_init('y',data(1:ny),units='meters',cart_name='y',domain2=domain)
   axes(3) = diag_axis_init('z',data(1:nz),units='mb',cart_name='z')
   call mpp_get_data_domain(domain,isd,ied,jsd,jed)
   call mpp_get_compute_domain(domain,isc,iec,jsc,jec)
   allocate(atmos_tracers(isd:ied,jsd:jed,nz,num_tracer_prog_atmos,1:3))
   allocate(atmos_tendency(isd:ied,jsd:jed,nz))
   allocate(atmos_vel(isd:ied,jsd:jed,nz,2))
   atmos_vel=advect_speed
   do i=1, num_tracer_prog_atmos
      if (mpp_pe() == mpp_root_pe()) write(*,'(a,i3,a)') 'Assigning tracer index ',atmos_prog_ind(i),' field to target array' 
      call get_tracer_names(MODEL_ATMOS,atmos_prog_ind(i),name,longname,units)
      atmos_prog_tracer_diagnostic_id(i) =  register_diag_field('tracers', trim(name),axes(1:3),model_time,&
                                trim(longname), trim(units), missing_value)      
      call assign_tracer_field(MODEL_ATMOS, atmos_prog_ind(i), data_tlevels=atmos_tracers(:,:,:,i,:))
   enddo
endif

! Allocate space for the oceanic prognostic tracers
if (num_tracer_prog_ocean > 0) then
   allocate(ocean_prog_ind(num_tracer_prog_ocean))
!   allocate(ocean_prog_tracer_diagnostic_id(num_tracer_prog_ocean))
   call get_tracer_indices(MODEL_OCEAN,prog_ind=ocean_prog_ind)

   halo=1 ! default halo size for 2nd order advection
   do n=1,num_tracer_prog_atmos
     ocean_order(n)=2! Default advection scheme
     flag = query_method ('horizontal-advection-scheme',MODEL_OCEAN,ocean_prog_ind(n),scheme_name,control)
     if(flag) then ! There is an advection scheme in the table so use that instead.
       select case (trim(scheme_name))
         case ('2nd_order') 
           ocean_order(n) = 2
         case ('4th_order')
           ocean_order(n) = 4
         case ('quicker') ! Have it be 4th order for simplicity.
           ocean_order(n) = 4
           halo = 2
         case default
           call mpp_error(FATAL,'invalid tracer advection scheme')
       end select
     endif
   enddo
endif

if (num_tracer_diag_atmos > 0) then
! Allocate space for the atmospheric diagnostic tracers
   allocate(atmos_diag_ind(num_tracer_diag_atmos))
   allocate(atmos_diag_tracer_diagnostic_id(num_tracer_diag_atmos))
   allocate(atmos_diag_tracers(isd:ied,jsd:jed,nz,num_tracer_diag_atmos))
   call get_tracer_indices(MODEL_ATMOS,diag_ind=atmos_diag_ind)
   color_index = get_tracer_index(MODEL_ATMOS,'color',atmos_diag_ind)
   do m=1,num_tracer_diag_atmos
      call get_tracer_names(MODEL_ATMOS,atmos_diag_ind(m),name,longname,units)
      atmos_diag_tracer_diagnostic_id(m) = register_diag_field('tracers',trim(name),axes(1:3),model_time,&
                         trim(longname),trim(units),missing_value)
   enddo
endif


end subroutine alloc_tracers

  subroutine lateral_advection(order, tendency, tracer_now, vel_now, tracer_prev, vel_prev, tracer_next, vel_next)

    integer, intent(in) :: order
    real, intent(inout), dimension(:,:,:) :: tendency, tracer_now
    real, intent(in), dimension(:,:,:,:) :: vel_now
    real, intent(in), dimension(:,:,:), optional :: tracer_prev, tracer_next
    real, intent(in), dimension(:,:,:,:), optional :: vel_prev, vel_next

    integer :: ni, nj, nk, i, j, k
    real :: a,b
    real, allocatable, dimension(:) :: fe, fn, fw, fs

    a=7.0/12.0;b=-1.0/12.0
    ni=size(tendency,1);nj=size(tendency,2);nk=size(tendency,3)

    allocate(fe(nk), fw(nk), fn(nk), fs(nk))

    select case (order)
    case (2)
       do j=2,nj-1
          do i=2,ni-1
             fe = 0.5*(tracer_now(i,j,:)+tracer_now(i+1,j,:))*vel_now(i+1,j,:,1)
             fn = 0.5*(tracer_now(i,j,:)+tracer_now(i,j+1,:))*vel_now(i,j+1,:,2)
             fw = 0.5*(tracer_now(i,j,:)+tracer_now(i-1,j,:))*vel_now(i,j,:,1) 
             fs = 0.5*(tracer_now(i,j,:)+tracer_now(i,j-1,:))*vel_now(i,j,:,2)
             tendency(i,j,:) = fe - fw + fn - fs
          enddo
       enddo
    case (4)
       do j=3,nj-2
          do i=3,ni-2
             fe = (a*(tracer_now(i,j,:)+tracer_now(i+1,j,:))+b*(tracer_now(i+2,j,:)+tracer_now(i-1,j,:)))*vel_now(i+1,j,:,1)
             fn = (a*(tracer_now(i,j,:)+tracer_now(i,j+1,:))+b*(tracer_now(i,j+2,:)+tracer_now(i,j-1,:)))*vel_now(i,j+1,:,2)
             fw = (a*(tracer_now(i,j,:)+tracer_now(i-1,j,:))+b*(tracer_now(i+1,j,:)+tracer_now(i-2,j,:)))*vel_now(i,j,:,1)
             fs = (a*(tracer_now(i,j,:)+tracer_now(i,j-1,:))+b*(tracer_now(i,j+1,:)+tracer_now(i,j-2,:)))*vel_now(i,j,:,2)
             tendency(i,j,:) = fe - fw + fn - fs
          enddo
       enddo
    case default
       call mpp_error(FATAL,'invalid advection scheme')
    end select

    deallocate(fe, fw, fn, fs)

    return
  end subroutine lateral_advection

  subroutine init_tracer(tracer, param, domain)

    real, intent(out), dimension(:,:,:) :: tracer
    real, intent(in)                    :: param
    type(domain2d), intent(in)          :: domain

    integer :: ni, nj, nk, i, j, isg, ieg, jsg, jeg, isc, iec, jsc, jec
    real :: xscale, yscale
    
    call mpp_get_global_domain(domain, xbegin=isg, xend=ieg, ybegin=jsg, yend=jeg, xsize=ni, ysize=nj)
    call mpp_get_compute_domain(domain,isc, iec, jsc, jec)

    xscale = ni/10
    yscale = nj/10

    do j=1,size(tracer,2)
       do i=1,size(tracer,1)
          tracer(i,j,:) = param*exp(-1*(isc+i-1-ni/2)**2/xscale**2)*exp(-1*(jsc+j-1-nj/2)**2/yscale**2)
       enddo
    enddo

    return
  end subroutine init_tracer

  subroutine calculate_color(temp,color)

    real, dimension(:,:,:), intent(in) :: temp
    real, dimension(:,:,:), intent(out) :: color

    color(:,:,:) = 0.5*temp(:,:,:)

    return

  end subroutine calculate_color

end program test
#endif


