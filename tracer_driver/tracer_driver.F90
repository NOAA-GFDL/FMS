
module tracer_driver_mod

!-----------------------------------------------------------------------

use utilities_mod, only:  file_exist, open_file, error_mesg,  &
                          FATAL, WARNING, NOTE, get_my_pe, close_file
use constants_mod, only:  grav

implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  tracer_driver, tracer_driver_init, tracer_driver_end

!-----------------------------------------------------------------------
!----------- namelist -------------------

integer :: nradon=0  ! tracer number for radon
integer :: ntest =0  ! starting tracer number for test tracers
                     ! test tracers are indexed from ntest to size(r,4)

namelist /tracer_driver_nml/ nradon

logical :: do_init=.true.

integer :: numtrace
!-----------------------------------------------------------------------
!---- version number -----
character(len=128) :: version = '$Id: tracer_driver.F90,v 1.3 2001/04/13 15:24:34 fms Exp $'
character(len=128) :: tag = '$Name: fez $'
!-----------------------------------------------------------------------
      integer, parameter :: max_tracers = 30
      character(len=16) :: field(max_tracers)

contains

!#######################################################################

 subroutine tracer_driver (lon, lat, land, phalf, r, rdt, kbot)

!-----------------------------------------------------------------------
   real, intent(in),    dimension(:,:)     :: lon, lat
logical, intent(in),    dimension(:,:)     :: land
   real, intent(in),    dimension(:,:,:)   :: phalf
   real, intent(in),    dimension(:,:,:,:) :: r
   real, intent(inout), dimension(:,:,:,:) :: rdt
integer, intent(in),    dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
real, dimension(size(r,1),size(r,2),size(r,3)) :: pwt, rtnd
integer :: k, kd, nt, ntr
!-----------------------------------------------------------------------

      if (do_init)  call error_mesg ('tracer_driver',  &
                    'tracer_driver_init must be called first.', FATAL)

      if (numtrace == 0) return

!-----------------------------------------------------------------------
      kd=size(r,3); nt=size(r,4)

      do k=1,kd
         pwt(:,:,k)=phalf(:,:,k+1)-phalf(:,:,k)
      enddo

!--------------- compute radon source-sink tendency --------------------

   if (nradon > 0 .and. nradon <= nt) then
      if (present(kbot)) then
         call radon_sourcesink (lon,lat,land,pwt,r(:,:,:,nradon),  &
                                rtnd,kbot)
      else
         call radon_sourcesink (lon,lat,land,pwt,r(:,:,:,nradon),rtnd)
      endif
      rdt(:,:,:,nradon)=rdt(:,:,:,nradon)+rtnd(:,:,:)
   endif

!-----------------------------------------------------------------------

 end subroutine tracer_driver

!#######################################################################

 subroutine radon_sourcesink (lon, lat, land, pwt, radon, radon_dt,  &
                              kbot)

!-----------------------------------------------------------------------
   real, intent(in),  dimension(:,:)   :: lon, lat
logical, intent(in),  dimension(:,:)   :: land
   real, intent(in),  dimension(:,:,:) :: pwt, radon
   real, intent(out), dimension(:,:,:) :: radon_dt
integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
   real, dimension(size(radon,1),size(radon,2),size(radon,3)) ::  &
         source, sink
logical, dimension(size(radon,1),size(radon,2)) ::  maskeq,masknh
   real  radon_flux, dtr, deg60, deg70, deg300, deg336
integer  i,j,kb,id,jd,kd
!-----------------------------------------------------------------------

      id=size(radon,1); jd=size(radon,2); kd=size(radon,3)

      dtr=acos(0.0)/90.
      deg60=60.*dtr; deg70=70.*dtr; deg300=300.*dtr; deg336=336.*dtr

!----------- compute radon source ------------
!
!  rn222 flux is 3.69e-21 kg/m*m/sec over land for latitudes lt 60n
!   between 60n and 70n the source  = source * .5
!
!  molecular wt. of air is 28.9644 gm/mole
!  molecular wt. of radon is 222 gm/mole
!  scaling facter to get reasonable mixing ratio is 1.e+21
!
!  source = 3.69e-21 * g * 28.9644 * 1.e+21/(pwt * 222.) or
!
!  source = g * .4814353 / pwt
!
!  must initialize all rn to .001
!

      radon_flux = 3.69e-21 * grav * 28.9644 * 1.e+21 / 222.
      source = 0.0
      maskeq = land .and. lat > -deg60 .and. lat < deg60
      masknh = land .and. lat >= deg60 .and. lat < deg70

      if (present(kbot)) then
          do j=1,jd
          do i=1,id
             kb=kbot(i,j)
             if (maskeq(i,j)) source(i,j,kb)=radon_flux/pwt(i,j,kb)
             if (masknh(i,j)) source(i,j,kb)=0.5*radon_flux/pwt(i,j,kb)
          enddo
          enddo
      else
          where (maskeq) source(:,:,kd)=radon_flux/pwt(:,:,kd)
          where (masknh) source(:,:,kd)=0.5*radon_flux/pwt(:,:,kd)
          where (masknh .and. lon > deg300 .and. lon < deg336)  &
               source(:,:,kd)=0.0
      endif

!------- compute radon sink --------------
!
!  rn222 has a half-life time of 3.83days 
!   (corresponds to an e-folding time of 5.52 days)
!
!  sink = 1./(86400.*5.52) = 2.09675e-6
!

!!    where (radon(:,:,:) >= 0.0)
         sink(:,:,:) = -2.09675e-6*radon(:,:,:)
!!    elsewhere
!!       sink(:,:,:) = 0.0
!!    endwhere

!------- tendency ------------------

      radon_dt=source+sink

!-----------------------------------------------------------------------

 end subroutine radon_sourcesink

!#######################################################################

 subroutine tracer_driver_init (g, r, mask)

!-----------------------------------------------------------------------
!
!   g    = gravitational constant (in m2/s2) NOT USED WILL BE REMOVED
!   r    = tracer fields dimensioned as (nlon,nlat,nlev,ntrace)
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
      real, intent(in) :: g
      real, intent(inout), dimension(:,:,:,:) :: r
      real, intent(in),    dimension(:,:,:), optional :: mask
!-----------------------------------------------------------------------
      integer  log_unit,unit,io,index
      character(len=16) ::  fld

!---- read namelist ----

      if ( file_exist('input.nml')) then
         unit = open_file ('input.nml', action='read')
         io=1; do while (io .ne. 0)
            read  (unit, nml=tracer_driver_nml, iostat=io, end=10)
         enddo
  10     call close_file (unit)
      endif

!---- write namelist ------------------

      log_unit = open_file ('logfile.out', action='append')
      if ( get_my_pe() == 0 ) then
           write (log_unit,'(/,80("="),/(a))') trim(version), trim(tag)
           write (log_unit,nml=tracer_driver_nml)
      endif

!----- read restart file ----------

      field = '                '

      if (file_exist('INPUT/tracer_driver.res')) then
         unit = open_file ('INPUT/tracer_driver.res', action='read')
         do
            read (unit,200,end=20) fld,index
            if (index <= max_tracers) field(index) = fld
         enddo
 200     format (a16,i3)
  20     call close_file (unit)
      endif

!---- dummy checks --------

      if (nradon > max_tracers) call error_mesg ('tracer_driver_init', &
                         'radon tracer number too large.', FATAL)

!----- count number of tracers ------

      numtrace = 0
      if (nradon > 0) numtrace = numtrace + 1

      if (numtrace == 0) then
          do_init = .false.
          go to 99   !return
      else
          if (size(r,4) == 0) call error_mesg ('tracer_driver_init',  &
                               'number of tracers equals zero.', FATAL)
      endif

!----- set initial value of radon ------------

      if (nradon > 0 ) then
         if (trim(field(nradon)) /= 'radon') then
            if (present(mask)) then
                r(:,:,:,nradon) = 0.001*mask(:,:,:)
            else
                r(:,:,:,nradon) = 0.001
            endif
            field(nradon) = 'radon'
            if (get_my_pe() == 0) write (log_unit,30) nradon
  30        format ('radon was initialized as tracer number ',i2)
         endif
      endif


  99  call close_file (log_unit)
      do_init = .false.

!-----------------------------------------------------------------------

 end subroutine tracer_driver_init

!#######################################################################

 subroutine tracer_driver_end

!-----------------------------------------------------------------------
   integer :: unit
   integer :: n

      if (numtrace    == 0) return
      if (get_my_pe() /= 0) return

!---------- open file and write restart file ---------------------------
!---------- file defines which tracers are used ------------------------

      unit = open_file ('RESTART/tracer_driver.res', action='write')
      do n = 1, numtrace 
        write (unit,100) field(n), n
      enddo
      call close_file (unit)

 100  format (a16,i3)

!-----------------------------------------------------------------------

 end subroutine tracer_driver_end

!#######################################################################

end module tracer_driver_mod

