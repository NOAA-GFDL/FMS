
               module column_diagnostics_mod



use mpp_io_mod,             only:  mpp_io_init, mpp_open, MPP_ASCII, &
                                   MPP_OVERWR, MPP_SEQUENTIAL,   &
                                   MPP_MULTI, mpp_close
use fms_mod,                only:  fms_init, mpp_pe, mpp_root_pe, &
                                   file_exist, check_nml_error, &
                                   error_mesg, FATAL, NOTE, WARNING, &
                                   close_file, open_namelist_file, &
                                   stdlog, write_version_number
use time_manager_mod,       only:  time_manager_init, month_name, &
                                   get_date, time_type
use constants_mod,          only:  constants_init, PI, RADIAN

!-------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!       module to locate and mark desired diagnostic columns         
!
!
!--------------------------------------------------------------------
  



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------


character(len=128)  :: version =  '$Id: column_diagnostics.F90,v 10.0 2003/10/24 22:01:26 fms Exp $'
character(len=128)  :: tag     =  '$Name: jakarta $'



!---------------------------------------------------------------------
!-------  interfaces --------

public    column_diagnostics_init,  &
          initialize_diagnostic_columns,  &
          column_diagnostics_header,   &
          close_column_diagnostics_units


!private 


!--------------------------------------------------------------------
!----    namelist -----

integer      :: dummy

namelist / column_diagnostics_nml /              &
                                      dummy

!--------------------------------------------------------------------
!-------- public data  -----


!--------------------------------------------------------------------
!------ private data ------


real, dimension(:), allocatable :: latb_deg, lonb_deg
logical    :: column_diagnostics_initialized = .false.
real       :: dellat, dellon

!-------------------------------------------------------------------
!-------------------------------------------------------------------



                        contains



!####################################################################

subroutine column_diagnostics_init (lonb_in, latb_in)

!--------------------------------------------------------------------
real, dimension(:), intent(in)   :: lonb_in, latb_in
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    column_diagnostics_init is the constructor for 
!    column_diagnostics_mod.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    local variables:
!
      integer    :: unit, ierr, io

!--------------------------------------------------------------------
!   local variables:
!
!       unit       unit number for nml file
!       ierr       error return flag
!       io         error return code
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    if routine has already been executed, return.
!--------------------------------------------------------------------
      if (column_diagnostics_initialized) return

!---------------------------------------------------------------------
!    verify that all modules used by this module have been initialized.
!----------------------------------------------------------------------
      call mpp_io_init
      call fms_init
      call time_manager_init
      call constants_init
 
!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read (unit, nml=column_diagnostics_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'column_diagnostics_nml')
        enddo
10      call close_file (unit)
      endif
 
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tag)
      if (mpp_pe() == mpp_root_pe())    &
                    write (stdlog(), nml=column_diagnostics_nml)

!--------------------------------------------------------------------
!    save the input lat and lon fields. define the delta of latitude
!    and longitude.
!--------------------------------------------------------------------
      allocate ( latb_deg (size(latb_in)) )
      allocate ( lonb_deg (size(lonb_in)) )
      latb_deg = latb_in*RADIAN
      lonb_deg = lonb_in*RADIAN
      dellat = latb_in(2) - latb_in(1)
      dellon = lonb_in(2) - lonb_in(1)


!--------------------------------------------------------------------
      column_diagnostics_initialized = .true.


end subroutine column_diagnostics_init 



!####################################################################


subroutine initialize_diagnostic_columns     &
                   (module, num_diag_pts_latlon, num_diag_pts_ij,  &
                    global_i , global_j , global_lat_latlon,   &
                    global_lon_latlon,  do_column_diagnostics,  &
                    diag_lon, diag_lat, diag_i, diag_j, diag_units)

!---------------------------------------------------------------------
!    initialize_diagnostic_columns returns the (i, j, lat, lon) coord-
!    inates of any diagnostic columns that are located on the current
!    processor.
!----------------------------------------------------------------------

!---------------------------------------------------------------------
character(len=*),      intent(in)    :: module
integer,               intent(in)    :: num_diag_pts_latlon,  &
                                        num_diag_pts_ij
integer, dimension(:), intent(in)    :: global_i, global_j   
real   , dimension(:), intent(in)    :: global_lat_latlon,    &
                                        global_lon_latlon 
logical, dimension(:), intent(out)   :: do_column_diagnostics
integer, dimension(:), intent(out)   :: diag_i, diag_j        
real   , dimension(:), intent(out)   :: diag_lat, diag_lon
integer, dimension(:), intent(out)   :: diag_units
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    intent(in) variables:
!
!       module                module calling this subroutine
!       num_diag_pts_latlon   number of diagnostic columns specified
!                             by lat-lon  coordinates
!       num_diag_pts_ij       number of diagnostic columns specified
!                             by global (i,j) coordinates
!       global_i              specified global i coordinates
!       global_j              specified global j coordinates
!       global_lat_latlon     specified global lat coordinates
!       global_lon_latlon     specified global lon coordinates
!
!    intent(out) variables:
!
!       do_column_diagnostics is a diagnostic column in this jrow ?
!       diag_i                processor i indices of diagnstic columns
!       diag_j                processor j indices of diagnstic columns
!       diag_lat              latitudes of diagnostic columns 
!                             [ degrees ]
!       diag_lon              longitudes of diagnostic columns 
!                             [ degrees ]
!       diag_units            unit number for each diagnostic column
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    local variables:

      real, dimension(size(diag_i)) :: global_lat, global_lon

      integer            ::  num_diag_pts
      integer            ::  i, j, nn
      character(len=8)   ::  char
      character(len=32)  ::   filename

!--------------------------------------------------------------------
!    local variables:
!
!       global_lat      latitudes for all diagnostic columns [ degrees ]
!       global_lon      longitudes for all diagnostic columns 
!                       [ degrees ]
!       num_diag_pts    total number of diagnostic columns
!       i, j, nn        do loop indices
!       char            character string for diaganostic column index
!       filename        filename for output file for diagnostic column 
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    initialize column_diagnostics flag and diag unit numbers. define 
!    total number of diagnostic columns.
!----------------------------------------------------------------------
      do_column_diagnostics = .false.
      diag_units(:) = -1
      diag_i(:) = -99
      diag_j(:) = -99
      diag_lat(:) = -999.
      diag_lon(:) = -999.
      num_diag_pts = size(diag_i)

!--------------------------------------------------------------------
!    define an array of lat-lon values for all diagnostic columns.
!--------------------------------------------------------------------
      do nn = 1, num_diag_pts_latlon
        global_lat(nn) = global_lat_latlon(nn)
        global_lon(nn) = global_lon_latlon(nn)
      end do

      do nn = 1, num_diag_pts_ij
        global_lat(nn+num_diag_pts_latlon) =    &
                         ((-0.5*acos(-1.0) + 0.5*dellat) + &
                         (global_j (nn)-1) *dellat)*RADIAN
        global_lon(nn+num_diag_pts_latlon) = (0.5*dellon +     &
                          (global_i (nn)-1)*dellon)*RADIAN
      end do   

!----------------------------------------------------------------------
!    loop over all diagnostic points to check for their presence on 
!    this processor.
!----------------------------------------------------------------------
      do nn=1,num_diag_pts
!----------------------------------------------------------------------
!    verify that the values of lat and lon are valid.
!----------------------------------------------------------------------
        if (global_lon(nn) >= 0. .and. global_lon(nn) <= 360.0) then
        else
          call error_mesg ('column_diagnostics_mod', &
               ' invalid longitude', FATAL)
        endif
        if (global_lat(nn) >= -90.0 .and. global_lat(nn) <= 90.0) then 
        else
          call error_mesg ('column_diagnostics_mod', &
               ' invalid latitude', FATAL)
        endif

!--------------------------------------------------------------------
!    determine if the diagnostics column is within the current 
!    processor's domain. if so, set a logical flag indicating the
!    presence of a diagnostic column on the particular row, define the 
!    i and j processor-coordinates and the latitude and longitude of 
!    the diagnostics column.
!--------------------------------------------------------------------
        do j=1,size(latb_deg) - 1
          if (global_lat(nn) .ge. latb_deg(j) .and.  &
              global_lat(nn) .lt. latb_deg(j+1)) then
            do i=1,size(lonb_deg) - 1
              if (global_lon(nn) .ge. lonb_deg(i)     .and.&
                  global_lon(nn) .lt. lonb_deg(i+1))  then
                do_column_diagnostics(j) = .true.
                diag_j(nn) = j
                diag_i(nn) = i
                diag_lon(nn) = 0.5*(lonb_deg(i) + lonb_deg(i+1))
                diag_lat(nn) = 0.5*(latb_deg(j) + latb_deg(j+1))
                write (char, '(i2)') nn
                filename = trim(module) // '_point' //    &
                           trim(adjustl(char)) // '.out'
                call mpp_open (diag_units(nn), filename, &
                               form=MPP_ASCII, &
                               action=MPP_OVERWR,  &
                               access=MPP_SEQUENTIAL,  &
                               threading=MPP_MULTI, nohdrs=.true.)
                exit
              endif
            end do
            exit
          endif
        end do
        
      end do

!---------------------------------------------------------------------


end subroutine initialize_diagnostic_columns




!####################################################################

subroutine column_diagnostics_header     &
                              (module, diag_unit, Time, nn, diag_lon, &
                               diag_lat, diag_i, diag_j)

!--------------------------------------------------------------------
!    column_diagnostics_header writes out information concerning
!    time and location of following data into the column_diagnostics
!    output file.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*),      intent(in)  :: module
type(time_type),       intent(in)  :: Time 
integer,               intent(in)  :: diag_unit
integer,               intent(in)  :: nn
real,    dimension(:), intent(in)  :: diag_lon, diag_lat
integer, dimension(:), intent(in)  :: diag_i, diag_j         

!--------------------------------------------------------------------
!    intent(in) variables
!
!       module     module name calling this subroutine
!       Time       current model time [ time_type ]
!       diag_unit  unit number for column_diagnostics output
!       nn         index of diagnostic column currently active
!       diag_lon   longitude of current diagnostic column [ degrees ]
!       diag_lat   latitude of current diagnostic column [ degrees ]
!       diag_i     i coordinate of current diagnostic column
!       diag_j     j coordinate of current diagnostic column
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!     local variables:

      integer           :: year, month, day, hour, minute, second 
      character(len=8)  :: mon
      character(len=64) :: header

!--------------------------------------------------------------------
!     local variables:
!    
!       year, month, day, hour, minute, seconds   
!                      integers defining the current time
!       mon            character string for the current month
!       header         title for the output 
!        
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    convert the time type to a date and time for printing. convert 
!    month to a character string.
!--------------------------------------------------------------------
      call get_date (Time, year, month, day, hour, minute, second)
      mon = month_name(month)

!---------------------------------------------------------------------
!    write timestamp and column location information to the diagnostic
!    columns output unit.
!---------------------------------------------------------------------
      write (diag_unit,'(a)')  ' '
      write (diag_unit,'(a)')  ' '
      write (diag_unit,'(a)')   &
              '======================================================'
      write (diag_unit,'(a)')  ' '
      header = '               PRINTING ' // module // '  DIAGNOSTICS' 
      write (diag_unit,'(a)')  header                          
      write (diag_unit,'(a)')  ' '
      write (diag_unit,'(a, i6,2x, a,i4,i4,i4,i4)')  ' time stamp:',  &
                                           year, trim(mon), day, &
                                           hour, minute, second
      write (diag_unit,'(a, i4)')      &
            ' DIAGNOSTIC POINT COORDINATES, point #', nn
      write (diag_unit,'(a)')  ' '
      write (diag_unit,'(a,f8.3,a,f8.3)') ' longitude = ',    &
                   diag_lon(nn), ' latitude  = ', diag_lat(nn)
      write (diag_unit,'(a, i4, a,i6,a,i6)')    &
                               ' on processor # ', mpp_pe(),   &
                               ' :   processor i =', diag_i(nn),     &
                               ' ,   processor j =', diag_j(nn)
      write (diag_unit,'(a)')  ' '

!---------------------------------------------------------------------



end subroutine column_diagnostics_header



!######################################################################

subroutine close_column_diagnostics_units (diag_units)

!---------------------------------------------------------------------
!    close_column_diagnostics_units closes any open column_diagnostics
!    files associated with the calling module.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
integer, dimension(:), intent(in)  :: diag_units
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!    intent(in) variable:
!
!      diag_units    array of column diagnostic unit numbers
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    local variable

      integer   :: nn    ! do loop index

!--------------------------------------------------------------------
!    close the unit associated with each diagnostic column.
!--------------------------------------------------------------------
      do nn=1, size(diag_units)
        if (diag_units(nn) /= -1) then
          call mpp_close (diag_units(nn))
        endif
      end do

!---------------------------------------------------------------------


end subroutine close_column_diagnostics_units


!#####################################################################




               end module column_diagnostics_mod
