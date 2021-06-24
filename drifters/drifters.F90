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
!FDOC_TAG_GFDL fdoc.pl generated xml skeleton

#include "fms_switches.h"
#define _FLATTEN(A) reshape((A), (/size((A))/) )

!> @defgroup drifters_mod drifters_mod
!> @ingroup drifters
!! @brief <TT>Drifters_mod</TT>is a module designed to advect a set of particles, in parallel or
!!   sequentially, given an prescribed velocity field.
!! @author Alexander Pletzer
!!
!> Drifters are idealized point particles with positions that evolve in time according
!! to a prescribed velocity field, starting from some initial conditions. Drifters have
!! no mass, no energy, no size, and no friction and therefore have no impact on the
!! dynamics of the underlying system. The only feature that distinguishes a drifter
!! from another is its trajectory. This makes drifters ideal for tracking pollution
!! clouds and probing fields (e.g. temperature, salinity) along ocean currents, to name
!! a few applications.
!! Drifters can mimic real experiments such as the Argo floats
!! http://www.metoffice.com/research/ocean/argo/ukfloats.html.
!!
!! When run in parallel, on a 2d decomposed domain, <TT>drifters_mod</TT> will handle all the
!! bookkeeping and communication transparently for the user. This involves adding/removing
!! drifters as they enter/leave a processor element (PE) domain. Note that the number of drifters
!! can vary greatly both between PE domains and within a PE domain in the course of a simulation; the drifters'
!! module will also manage dynamically the memory for the user.
!!
!! There are a number of basic assumptions which could make the drifters' module
!! ill-suited for some tasks. First and foremost, it is assumed that the motion of
!! drifters is not erratic but follows deterministic trajectories. Furthermore,
!! drifters should not cross both compute and data domain boundaries within less
!! than a time step. This limitation is imposed by the Runge-Kutta integration
!! scheme, which must be able to complete, within a time step, a trajectory
!! calculation that starts inside the compute domain and ends inside the data domain. Therefore, the drifters,
!! as they are presently modelled, are unlikely to work for very fast objects.
!! This constraint also puts a upper limit to the domain decomposition, although
!! it can often be remedied by increasing the number of ghost nodes.
!!
!! Another fundamental assumption is that the (e.g. velocity) fields are structured,
!! on a per PE domain basis. There is no support for locally nested or unstrucured
!! meshes. Meshes need not be smooth and continuous across PE domains, however.

!> @file
!> @brief File for @ref drifters_mod

!> @addtogroup drifters_mod
!> @{
module drifters_mod

#ifdef _SERIAL

! serial code
#define _MPP_PE 0
#define _MPP_ROOT 0
#define _MPP_NPES 1
#define _TYPE_DOMAIN2D integer

#else

! parallel code
  use mpp_mod        , only : mpp_pe, mpp_npes
  use mpp_domains_mod, only : domain2d
#define _MPP_PE mpp_pe()
#define _MPP_ROOT mpp_root_pe()
#define _MPP_NPES mpp_npes()
#define _TYPE_DOMAIN2D type(domain2d)

#endif

  use drifters_core_mod,  only: drifters_core_type, drifters_core_new, drifters_core_del, assignment(=)

  use drifters_input_mod, only: drifters_input_type, drifters_input_new, drifters_input_del, assignment(=)

  use drifters_io_mod,    only: drifters_io_type, drifters_io_new, drifters_io_del, drifters_io_set_time_units, &
                                drifters_io_set_position_names, drifters_io_set_position_units, &
                                drifters_io_set_field_names, drifters_io_set_field_units, drifters_io_write

  use drifters_comm_mod,  only: drifters_comm_type, drifters_comm_new, drifters_comm_del, drifters_comm_set_pe_neighbors, &
                                drifters_comm_set_domain, drifters_comm_gather, drifters_comm_update

  use cloud_interpolator_mod, only: cld_ntrp_linear_cell_interp, cld_ntrp_locate_cell, cld_ntrp_get_cell_values
  implicit none
  private

  public :: drifters_type, assignment(=), drifters_push, drifters_compute_k, drifters_set_field
  public :: drifters_new, drifters_del, drifters_set_domain, drifters_set_pe_neighbors
  public :: drifters_set_v_axes, drifters_set_domain_bounds, drifters_positions2lonlat
  public :: drifters_print_checksums, drifters_save, drifters_write_restart, drifters_distribute

  integer, parameter, private :: MAX_STR_LEN = 128
! Include variable "version" to be written to log file.
#include<file_version.h>
  real :: DRFT_EMPTY_ARRAY(0)
  !> @}

  !> @brief Holds all data needed for drifters communication, io, and input.
  !> @ingroup drifters_mod
  type drifters_type
     ! Be sure to update drifters_new, drifters_del and drifters_copy_new
     ! when adding members
     type(drifters_core_type)  :: core
     type(drifters_input_type) :: input
     type(drifters_io_type)    :: io
     type(drifters_comm_type)  :: comm
     real    :: dt             !< total dt, over a complete step
     real    :: time
     ! fields
     real, allocatable :: fields(:,:)
     ! velocity field axes
     real, allocatable :: xu(:) !< velocity field axes
     real, allocatable :: yu(:) !< velocity field axes
     real, allocatable :: zu(:) !< velocity field axes
     real, allocatable :: xv(:) !< velocity field axes
     real, allocatable :: yv(:) !< velocity field axes
     real, allocatable :: zv(:) !< velocity field axes
     real, allocatable :: xw(:) !< velocity field axes
     real, allocatable :: yw(:) !< velocity field axes
     real, allocatable :: zw(:) !< velocity field axes
     ! Runge Kutta coefficients holding intermediate results (positions)
     real, allocatable :: temp_pos(:,:) !< Runge Kutta coefficients holding
                                        !! intermediate results (positions)
     real, allocatable :: rk4_k1(:,:) !< Runge Kutta coefficients holding
                                      !! intermediate results (positions)
     real, allocatable :: rk4_k2(:,:) !< Runge Kutta coefficients holding
                                      !! intermediate results (positions)
     real, allocatable :: rk4_k3(:,:) !< Runge Kutta coefficients holding
                                      !! intermediate results (positions)
     real, allocatable :: rk4_k4(:,:) !< Runge Kutta coefficients holding
                                      !! intermediate results (positions)
     ! store filenames for convenience
     character(len=MAX_STR_LEN) :: input_file !< store filenames for convenience
     character(len=MAX_STR_LEN) :: output_file !< store filenames for convenience
     ! Runge Kutta stuff
     integer :: rk4_step !< Runge Kutta stuff
     logical :: rk4_completed !< Runge Kutta stuff
     integer :: nx, ny
     logical, allocatable   :: remove(:)
  end type drifters_type

  !> @brief Assignment override for @ref drifters_type
  !> @ingroup drifters_mod
  interface assignment(=)
     module procedure drifters_copy_new
  end interface

  !> @brief "Push" a given drifter at a given velocity for either 2D or 3D data
  !> @ingroup drifters_mod
  interface drifters_push
    module procedure drifters_push_2
    module procedure drifters_push_3
  end interface

  !> @ingroup drifters_mod
  interface drifters_compute_k
     module procedure drifters_computek2d
     module procedure drifters_computek3d
  end interface

  !> @brief Set the value of a given drifter field
  !> @ingroup drifters_mod
  interface drifters_set_field
    module procedure drifters_set_field_2d
    module procedure drifters_set_field_3d
  end interface

!> @addtogroup drifters_mod
!> @{

contains

  !> @brief Will read positions stored in the netCDF file <TT>input_file</TT>.
  !!
  !>   The trajectories will be saved in files <TT>output_file.PE</TT>,
  !!   one file per PE domain.
  subroutine drifters_new(self, input_file, output_file, ermesg)

    type(drifters_type) :: self !< Opaque data structure.
    character(len=*), intent(in)  :: input_file !< NetCDF input file name containing initial positions.
    character(len=*), intent(in)  :: output_file !< NetCDF output file. Will contain trajectory
                                                 !! positions and interpolated fields.
    character(len=*), intent(out) :: ermesg !< Error message (if any).

    integer nd, nf, npdim, i
    character(len=6) :: pe_str

    ermesg = ''

    self%input_file  = input_file
    self%output_file = output_file

    call drifters_input_new(self%input, input_file, ermesg)
    if(ermesg/='') return

    ! number of dimensions
    nd = size(self%input%velocity_names)
    ! estimate for the max number of particles (will resize if exceeded)
    npdim = int(1.3*size(self%input%positions, 2))
    call drifters_core_new(self%core, nd=nd, npdim=npdim, ermesg=ermesg)
    if(ermesg/='') return

    ! number of fields
    nf = size(self%input%field_names)

    ! one output file per PE
    pe_str = '    '
    write(pe_str, '(i6)') _MPP_PE
    pe_str = adjustr(pe_str)
    do i = 1, 5
       if(pe_str(i:i)==' ') pe_str(i:i)='0'
    enddo
    call drifters_io_new(self%io, output_file//'.'//pe_str, nd, nf, ermesg)
    if(ermesg/='') return

    call drifters_comm_new(self%comm)
    if(ermesg/='') return

    ! Set meta data
    call drifters_io_set_time_units(self%io, name=self%input%time_units, &
         & ermesg=ermesg)

    call drifters_io_set_position_names(self%io, names=self%input%position_names, &
         & ermesg=ermesg)
    if(ermesg/='') return
    call drifters_io_set_position_units(self%io, names=self%input%position_units, &
         & ermesg=ermesg)
    if(ermesg/='') return

    call drifters_io_set_field_names(self%io, names=self%input%field_names, &
         & ermesg=ermesg)
    if(ermesg/='') return
    call drifters_io_set_field_units(self%io, names=self%input%field_units, &
         & ermesg=ermesg)
    if(ermesg/='') return

    self%dt   = -1
    self%time = -1
    self%rk4_step = 0
    self%nx       = 0
    self%ny       = 0
    self%rk4_completed = .FALSE.

    allocate(self%rk4_k1(self%core%nd, self%core%npdim))
    self%rk4_k1 = -huge(1.)
    allocate(self%rk4_k2(self%core%nd, self%core%npdim))
    self%rk4_k2 = -huge(1.)
    allocate(self%rk4_k3(self%core%nd, self%core%npdim))
    self%rk4_k3 = -huge(1.)
    allocate(self%rk4_k4(self%core%nd, self%core%npdim))
    self%rk4_k4 = -huge(1.)
    allocate(self%remove(self%core%npdim))
    self%remove = .FALSE.
    allocate(self%temp_pos(nd, self%core%npdim))
    self%temp_pos = -huge(1.)

    allocate(self%fields(nf, self%core%npdim))
    self%fields = -huge(1.)

  end subroutine drifters_new

  !============================================================================

  !> @brief Destructor, call this to reclaim memory from data used for drifters.
  subroutine drifters_del(self, ermesg)
    type(drifters_type) :: self !< Opaque data structure.
    character(len=*), intent(out) :: ermesg !< Error message (if any).

    integer flag
    ermesg = ''
    deallocate(self%fields, stat=flag)
    deallocate(self%xu, stat=flag)
    deallocate(self%yu, stat=flag)
    deallocate(self%zu, stat=flag)
    deallocate(self%xv, stat=flag)
    deallocate(self%yv, stat=flag)
    deallocate(self%zv, stat=flag)
    deallocate(self%xw, stat=flag)
    deallocate(self%yw, stat=flag)
    deallocate(self%zw, stat=flag)
    deallocate(self%temp_pos, stat=flag)
    deallocate(self%rk4_k1, stat=flag)
    deallocate(self%rk4_k2, stat=flag)
    deallocate(self%rk4_k3, stat=flag)
    deallocate(self%rk4_k4, stat=flag)
    deallocate(self%remove, stat=flag)

    call drifters_core_del(self%core, ermesg)
    if(ermesg/='') return
    call drifters_input_del(self%input, ermesg)
    if(ermesg/='') return
    call drifters_io_del(self%io, ermesg)
    if(ermesg/='') return
    call drifters_comm_del(self%comm)
    if(ermesg/='') return

  end subroutine drifters_del

  !============================================================================
  !> @brief Copy a drifter state into a new state. Note: this will not open new files; this will
  !!   copy all members into a new container.
  subroutine drifters_copy_new(new_instance, old_instance)

    type(drifters_type), intent(in)    :: old_instance !< Old data structure.
    type(drifters_type), intent(inout) :: new_instance !< New data structure.

    character(len=MAX_STR_LEN) :: ermesg

    ermesg = ''

    ! make sure new_instance is empty
    call drifters_del(new_instance, ermesg)
    if(ermesg/='') return

    new_instance%core  = old_instance%core
    new_instance%input = old_instance%input
    new_instance%io    = old_instance%io
    new_instance%comm  = old_instance%comm

     new_instance%dt     = old_instance%dt
     new_instance%time   = old_instance%time

     allocate(new_instance%fields( size(old_instance%fields, 1), &
          &                        size(old_instance%fields, 2) ))
     new_instance%fields = old_instance%fields

     allocate(new_instance%xu( size(old_instance%xu) ))
     allocate(new_instance%yu( size(old_instance%yu) ))
     allocate(new_instance%zu( size(old_instance%zu) ))
     new_instance%xu = old_instance%xu
     new_instance%yu = old_instance%yu
     new_instance%zu = old_instance%zu
     allocate(new_instance%xv( size(old_instance%xv) ))
     allocate(new_instance%yv( size(old_instance%yv) ))
     allocate(new_instance%zv( size(old_instance%zv) ))
     new_instance%xv = old_instance%xv
     new_instance%yv = old_instance%yv
     new_instance%zv = old_instance%zv
     allocate(new_instance%xw( size(old_instance%xw) ))
     allocate(new_instance%yw( size(old_instance%yw) ))
     allocate(new_instance%zw( size(old_instance%zw) ))
     new_instance%xw = old_instance%xw
     new_instance%yw = old_instance%yw
     new_instance%zw = old_instance%zw

     allocate(new_instance%temp_pos( size(old_instance%temp_pos,1), &
          &                          size(old_instance%temp_pos,2) ))
     new_instance%temp_pos = old_instance%temp_pos
     allocate(new_instance%rk4_k1( size(old_instance%rk4_k1,1), &
          &                        size(old_instance%rk4_k1,2) ))
     allocate(new_instance%rk4_k2( size(old_instance%rk4_k2,1), &
          &                        size(old_instance%rk4_k2,2) ))
     allocate(new_instance%rk4_k3( size(old_instance%rk4_k3,1), &
          &                        size(old_instance%rk4_k3,2) ))
     allocate(new_instance%rk4_k4( size(old_instance%rk4_k4,1), &
          &                        size(old_instance%rk4_k4,2) ))
     new_instance%rk4_k1 = old_instance%rk4_k1
     new_instance%rk4_k2 = old_instance%rk4_k2
     new_instance%rk4_k3 = old_instance%rk4_k3
     new_instance%rk4_k4 = old_instance%rk4_k4

     new_instance%rk4_step = old_instance%rk4_step
     new_instance%rk4_completed = old_instance%rk4_completed
     new_instance%nx = old_instance%nx
     new_instance%ny = old_instance%ny

     allocate(new_instance%remove(size(old_instance%remove)))
     new_instance%remove = old_instance%remove


  end subroutine drifters_copy_new

  !============================================================================
  !> @brief Set the compute, data, and global domain boundaries.
  !! @details The data domain extends beyond the compute domain and is shared between
  !!   two or more PE domains. A particle crossing the compute domain boundary
  !!   will trigger a communication with one or more neighboring domains. A particle
  !!   leaving the data domain will be removed from the list of particles.
  !!
  !! <br>Example usage:
  !! @code{.F90}
  !!   call   drifters_set_domain(self, &
  !!     & xmin_comp, xmax_comp, ymin_comp, ymax_comp, &
  !!     & xmin_data, xmax_data, ymin_data, ymax_data, &
  !!     & xmin_glob, xmax_glob, ymin_glob, ymax_glob, &
  !!     & ermesg)
  !! @endcode
  subroutine drifters_set_domain(self, &
       & xmin_comp, xmax_comp, ymin_comp, ymax_comp, &
       & xmin_data, xmax_data, ymin_data, ymax_data, &
       & xmin_glob, xmax_glob, ymin_glob, ymax_glob, &
       & ermesg)
    type(drifters_type) :: self !< Opaque data structure.
    ! compute domain boundaries
    real, optional, intent(in) :: xmin_comp !< Min of longitude-like axis on compute domain.
    real, optional, intent(in) :: xmax_comp !< Max of longitude-like axis on compute domain.
    real, optional, intent(in) :: ymin_comp !< Min of latitude-like axis on compute domain.
    real, optional, intent(in) :: ymax_comp !< Max of latitude-like axis on compute domain.
    ! data domain boundaries
    real, optional, intent(in) :: xmin_data !< Min of longitude-like axis on data domain.
    real, optional, intent(in) :: xmax_data !< Max of longitude-like axis on data domain.
    real, optional, intent(in) :: ymin_data !< Min of latitude-like axis on data domain.
    real, optional, intent(in) :: ymax_data !< Max of latitude-like axis on data domain.
    ! global boundaries (only specify those if domain is periodic)
    real, optional, intent(in) :: xmin_glob !< Min of longitude-like axis on global domain.
    real, optional, intent(in) :: xmax_glob !< Max of longitude-like axis on global domain.
    real, optional, intent(in) :: ymin_glob !< Min of latitude-like axis on global domain.
    real, optional, intent(in) :: ymax_glob !< Max of latitude-like axis on global domain.
    character(len=*), intent(out) :: ermesg !< Error message (if any).

    ermesg = ''
    if(present(xmin_comp)) self%comm%xcmin = xmin_comp
    if(present(xmax_comp)) self%comm%xcmax = xmax_comp
    if(present(ymin_comp)) self%comm%ycmin = ymin_comp
    if(present(ymax_comp)) self%comm%ycmax = ymax_comp

    if(present(xmin_data)) self%comm%xdmin = xmin_data
    if(present(xmax_data)) self%comm%xdmax = xmax_data
    if(present(ymin_data)) self%comm%ydmin = ymin_data
    if(present(ymax_data)) self%comm%ydmax = ymax_data

    if(present(xmin_glob)) self%comm%xgmin = xmin_glob
    if(present(xmax_glob)) self%comm%xgmax = xmax_glob
    if(present(ymin_glob)) self%comm%ygmin = ymin_glob
    if(present(ymax_glob)) self%comm%ygmax = ymax_glob

    ! Note: the presence of both xgmin/xgmax will automatically set the
    ! periodicity flag
    if(present(xmin_glob) .and. present(xmax_glob)) self%comm%xperiodic = .TRUE.
    if(present(ymin_glob) .and. present(ymax_glob)) self%comm%yperiodic = .TRUE.

  end subroutine drifters_set_domain

  !============================================================================
  !> @brief Given an MPP based deomposition, set the PE numbers that are adjacent to this
  !!   processor.
  !!
  !> This will allow several PEs to track the trajectories of particles in the buffer regions.
  subroutine drifters_set_pe_neighbors(self, domain, ermesg)

    type(drifters_type) :: self !< Opaque data structure.
    _TYPE_DOMAIN2D      :: domain !< MPP domain.
    character(len=*), intent(out) :: ermesg !< Error message (if any).

    ermesg = ''

    call drifters_comm_set_pe_neighbors(self%comm, domain)

  end subroutine drifters_set_pe_neighbors

  !============================================================================
#define _DIMS 2
#define drifters_push_XXX drifters_push_2
#include "drifters_push.h"
#undef _DIMS
#undef drifters_push_XXX

  !============================================================================
#define _DIMS 3
#define drifters_push_XXX drifters_push_3
#include "drifters_push.h"
#undef _DIMS
#undef drifters_push_XXX

  !============================================================================
  subroutine drifters_modulo(self, positions, ermesg)
    type(drifters_type) :: self
    real, intent(inout) :: positions(:,:)
    character(len=*), intent(out) :: ermesg

    integer ip, np
    real x, y

    ermesg = ''
    np = self%core%np

    if(self%comm%xperiodic) then
       do ip = 1, np
          x = positions(1, ip)
          positions(1, ip) = self%comm%xgmin + &
               & modulo(x - self%comm%xgmin, self%comm%xgmax-self%comm%xgmin)
       enddo
    endif

    if(self%comm%yperiodic) then
       do ip = 1, np
          y = positions(2, ip)
          positions(2, ip) = self%comm%ygmin + &
               & modulo(y - self%comm%ygmin, self%comm%ygmax-self%comm%ygmin)
       enddo
    endif

  end subroutine drifters_modulo

  !============================================================================
#define _DIMS 2
#define drifters_set_field_XXX drifters_set_field_2d
#include "drifters_set_field.h"
#undef _DIMS
#undef drifters_set_field_XXX

  !============================================================================
#define _DIMS 3
#define drifters_set_field_XXX drifters_set_field_3d
#include "drifters_set_field.h"
#undef _DIMS
#undef drifters_set_field_XXX
  !============================================================================
  !> @brief Append new positions to NetCDF file.
  !!
  !> Use this method to append the new trajectory positions and the interpolated
  !! probe fields to a netCDF file.
  subroutine drifters_save(self, ermesg)
    type(drifters_type) :: self !< Opaque daata structure.
    character(len=*), intent(out) :: ermesg !< Error message (if any).

    integer nf, np

    ermesg = ''
    nf = size(self%input%field_names)
    np = self%core%np

    ! save to disk
    call drifters_io_write(self%io, self%time, np, self%core%nd, nf, &
         & self%core%ids, self%core%positions, &
         & fields=self%fields(:,1:np), ermesg=ermesg)

  end subroutine drifters_save
  !============================================================================

  !> @brief Distribute particles across PEs.
  !!
  !> Use this method after setting the domain boundaries
  !! (<TT>drifters_set_domain</TT>) to spread the particles across PE domains.
  subroutine drifters_distribute(self, ermesg)
    type(drifters_type) :: self !< Opaque handle.
    character(len=*), intent(out) :: ermesg !< Error message (if any).

    real x, y
    integer i, nptot, nd

    ermesg = ''
    nd = self%core%nd
    if(nd < 2) then
       ermesg = 'drifters_distribute: dimension must be >=2'
       return
    endif

    nptot = size(self%input%positions, 2)
    do i = 1, nptot
       x = self%input%positions(1,i)
       y = self%input%positions(2,i)
       if(x >= self%comm%xdmin .and. x <= self%comm%xdmax .and. &
        & y >= self%comm%ydmin .and. y <= self%comm%ydmax) then

          self%core%np = self%core%np + 1
          self%core%positions(1:nd, self%core%np) = self%input%positions(1:nd, i)
          self%core%ids(self%core%np)             = i

       endif
    enddo

  end subroutine drifters_distribute

  !============================================================================
  !> @brief Write restart file for drifters.
  !!
  !> Gather all the particle positions distributed
  !! across PE domains on root PE and save the data in netCDF file.
  subroutine drifters_write_restart(self, filename, &
       & x1, y1, geolon1, &
       & x2, y2, geolat2, &
       & root, mycomm, ermesg)
    ! gather all positions and ids and save the result in
    ! self%input data structure on PE "root", then write restart file

    type(drifters_type) :: self !< Opaque data structure.
    character(len=*), intent(in)  :: filename !< Restart file name.

    ! if these optional arguments are passed, the positions will
    ! mapped to lon/lat degrees and saved in the file.
    real, intent(in), optional    :: x1(:) !< Pseudo-longitude axis supporting longitudes.
    real, intent(in), optional    :: y1(:) !< Pseudo-latitude axis supporting longitudes.
    real, intent(in), optional    :: geolon1(:,:) !< Longitude array (x1, y1).
    real, intent(in), optional    :: x2(:) !< Pseudo-longitude axis supporting latitudes.
    real, intent(in), optional    :: y2(:) !< Pseudo-latitude axis supporting latitudes.
    real, intent(in), optional    :: geolat2(:,:) !< Latitudes array (x2, y2)


    integer, intent(in), optional :: root    !< root pe
    integer, intent(in), optional :: mycomm  !< MPI communicator
    character(len=*), intent(out) :: ermesg !< Error message (if any).

    integer :: np
    logical :: do_save_lonlat
    real, allocatable    ::  lons(:), lats(:)

    ermesg = ''

    np = self%core%np

    allocate(lons(np), lats(np))
    lons = -huge(1.)
    lats = -huge(1.)

    ! get lon/lat if asking for
    if(present(x1) .and. present(y1) .and. present(geolon1) .and. &
         & present(x2) .and. present(y2) .and. present(geolat2)) then
       do_save_lonlat = .TRUE.
    else
       do_save_lonlat = .FALSE.
    endif

    if(do_save_lonlat) then

       ! Interpolate positions onto geo longitudes/latitudes
       call drifters_positions2lonlat(self,   &
            & positions=self%core%positions(:,1:np), &
            & x1=x1, y1=y1, geolon1=geolon1,         &
            & x2=x2, y2=y2, geolat2=geolat2,         &
            & lons=lons, lats=lats, ermesg=ermesg)
       if(ermesg/='') return ! problems, bail off

    endif

    call drifters_comm_gather(self%comm, self%core, self%input, &
         & lons, lats, do_save_lonlat, &
         & filename, &
         & root, mycomm)

  end subroutine drifters_write_restart

  !============================================================================
#define _DIMS 2
#define drifters_compute_k_XXX drifters_computek2d
#include "drifters_compute_k.h"
#undef _DIMS
#undef drifters_compute_k_XXX

   !============================================================================
#define _DIMS 3
#define drifters_compute_k_XXX drifters_computek3d
#include "drifters_compute_k.h"
#undef _DIMS
#undef drifters_compute_k_XXX


 !============================================================================
  !> @brief Set velocity field axes.
  !! @details Velocity axis components may be located on different grids or cell faces. For instance, zonal (u)
  !!  and meridional (v) velcity components are staggered by half a cell size in Arakawa's C and D grids.
  !!  This call will set individual axes for each components do as to allow interpolation of the velocity
  !!  field on arbitrary positions.
  subroutine drifters_set_v_axes(self, component, x, y, z, ermesg)
    type(drifters_type) :: self !< Opaque data structure.
    character(len=*), intent(in)  :: component !< Velocity component: either 'u', 'v', or 'w'.
    real, intent(in)              :: x(:) !< X-axis.
    real, intent(in)              :: y(:) !< Y-axis.
    real, intent(in)              :: z(:) !< Z-axis.
    character(len=*), intent(out) :: ermesg !< Error message (if any).

    integer ier, nx, ny, nz

    ermesg = ''
    nx = size(x)
    ny = size(y)
    nz = size(z)
    select case (component(1:1))
       case ('u', 'U')
          if(nx > 0) then
             deallocate(self%xu, stat=ier)
             allocate(self%xu(nx))
             self%xu = x
             self%nx = max(self%nx, size(x))
          endif
          if(ny > 0) then
             deallocate(self%yu, stat=ier)
             allocate(self%yu(ny))
             self%yu = y
             self%ny = max(self%ny, size(y))
          endif
          if(nz > 0) then
             deallocate(self%zu, stat=ier)
             allocate(self%zu(nz))
             self%zu = z
          endif
      case ('v', 'V')
          if(nx > 0) then
             deallocate(self%xv, stat=ier)
             allocate(self%xv(nx))
             self%xv = x
             self%nx = max(self%nx, size(x))
          endif
          if(ny > 0) then
             deallocate(self%yv, stat=ier)
             allocate(self%yv(ny))
             self%yv = y
             self%ny = max(self%ny, size(y))
          endif
          if(nz > 0) then
             deallocate(self%zv, stat=ier)
             allocate(self%zv(nz))
             self%zv = z
          endif
      case ('w', 'W')
          if(nx > 0) then
             deallocate(self%xw, stat=ier)
             allocate(self%xw(nx))
             self%xw = x
             self%nx = max(self%nx, size(x))
          endif
          if(ny > 0) then
             deallocate(self%yw, stat=ier)
             allocate(self%yw(ny))
             self%yw = y
             self%ny = max(self%ny, size(y))
          endif
          if(nz > 0) then
             deallocate(self%zw, stat=ier)
             allocate(self%zw(nz))
             self%zw = z
          endif
      case default
         ermesg = 'drifters_set_v_axes: ERROR component must be "u", "v" or "w"'
    end select
  end subroutine drifters_set_v_axes

  !============================================================================
  !> @brief Set boundaries of "data" and "compute" domains
  !! @details Each particle will be tracked sol long is it is located in the data domain.
  subroutine drifters_set_domain_bounds(self, domain, backoff_x, backoff_y, ermesg)
    type(drifters_type) :: self !< Opaque data structure.
    _TYPE_DOMAIN2D      :: domain !< Instance of Domain2D (see mpp_domain)
    integer, intent(in) ::  backoff_x !< particles leaves domain when crossing ied-backoff_x
    integer, intent(in) ::  backoff_y !< particles leaves domain when crossing jed-backoff_y
    character(len=*), intent(out) :: ermesg !< Error message (if any).

    ermesg = ''

    if(.not.allocated(self%xu) .or. .not.allocated(self%yu)) then
       ermesg = 'drifters_set_domain_bounds: ERROR "u"-component axes not set'
       return
    endif
    call drifters_comm_set_domain(self%comm, domain, self%xu, self%yu, backoff_x, backoff_y)
    if(.not.allocated(self%xv) .or. .not.allocated(self%yv)) then
       ermesg = 'drifters_set_domain_bounds: ERROR "v"-component axes not set'
       return
    endif
    if(allocated(self%xw) .and. allocated(self%yw)) then
       call drifters_comm_set_domain(self%comm, domain, self%xv, self%yv, backoff_x, backoff_y)
    endif


  end subroutine drifters_set_domain_bounds

  !============================================================================
  !> @brief Interpolates positions onto longitude/latitude grid.
  !! @details In many cases, the integrated positions will not be longitudes  or latitudes. This call
  !!  can be ionvoked to recover the longitude/latitude positions from the "logical" positions.
  subroutine drifters_positions2lonlat(self, positions, &
       &                                        x1, y1, geolon1, &
       &                                        x2, y2, geolat2, &
       &                                        lons, lats, &
       &                                        ermesg)

    type(drifters_type) :: self !< Opaque data structure.
    ! Input positions
    real, intent(in)    :: positions(:,:) !< Logical positions.
    ! Input mesh
    real, intent(in)    :: x1(:) !< X-axis of "geolon1" field.
    real, intent(in)    :: y1(:) !< Y-axis of "geolon1" field.
    real, intent(in)    :: geolon1(:,:) !< Y-axis of "geolon1" field.
    real, intent(in)    :: x2(:) !< X-axis of "geolat2" field.
    real, intent(in)    :: y2(:) !< Y-axis of "geolat2" field.
    real, intent(in)    :: geolat2(:,:) !< Latitude field as an array of (x2, y2)
    ! Output lon/lat
    real, intent(out)   :: lons(:) !< Returned longitudes.
    real, intent(out)   :: lats(:) !< Returned latitudes.
    character(len=*), intent(out) :: ermesg !< Error message (if any).

    real    fvals(2**self%core%nd), ts(self%core%nd)
    integer np, ij(2), ip, ier, n1s(2), n2s(2), i, j, iertot
    character(len=10) :: n1_str, n2_str, np_str, iertot_str

    ermesg = ''
    lons = -huge(1.)
    lats = -huge(1.)

    ! check dimensions
    n1s = (/size(x1), size(y1)/)
    n2s = (/size(x2), size(y2)/)
    if(n1s(1) /= size(geolon1, 1) .or. n1s(2) /= size(geolon1, 2)) then
       ermesg = 'drifters_positions2geolonlat: ERROR incompatibles dims between (x1, y1, geolon1)'
       return
    endif
    if(n2s(1) /= size(geolat2, 1) .or. n2s(2) /= size(geolat2, 2)) then
       ermesg = 'drifters_positions2geolonlat: ERROR incompatibles dims between (x2, y2, geolat2)'
       return
    endif

    np = size(positions, 2)
    if(size(lons) < np .or. size(lats) < np) then
       write(np_str, '(i10)') np
       write(n1_str, '(i10)') size(lons)
       write(n2_str, '(i10)') size(lats)
       ermesg = 'drifters_positions2geolonlat: ERROR size of "lons" ('//trim(n1_str)// &
            & ') or "lats" ('//trim(n2_str)//') < '//trim(np_str)
       return
    endif

    ! Interpolate
    iertot = 0
    do ip = 1, np

       ! get longitude
       call cld_ntrp_locate_cell(x1, positions(1,ip), i, ier)
       iertot = iertot + ier
       call cld_ntrp_locate_cell(y1, positions(2,ip), j, ier)
       iertot = iertot + ier
       ij(1) = i; ij(2) = j;
       call cld_ntrp_get_cell_values(n1s, _FLATTEN(geolon1), ij, fvals, ier)
       iertot = iertot + ier
       ts(1) = (positions(1,ip) - x1(i))/(x1(i+1) - x1(i))
       ts(2) = (positions(2,ip) - y1(j))/(y1(j+1) - y1(j))
       call cld_ntrp_linear_cell_interp(fvals, ts, lons(ip), ier)
       iertot = iertot + ier

       ! get latitude
       call cld_ntrp_locate_cell(x2, positions(1,ip), i, ier)
       iertot = iertot + ier
       call cld_ntrp_locate_cell(y2, positions(2,ip), j, ier)
       iertot = iertot + ier
       ij(1) = i; ij(2) = j;
       call cld_ntrp_get_cell_values(n2s, _FLATTEN(geolat2), ij, fvals, ier)
       iertot = iertot + ier
       ts(1) = (positions(1,ip) - x2(i))/(x2(i+1) - x2(i))
       ts(2) = (positions(2,ip) - y2(j))/(y2(j+1) - y2(j))
       call cld_ntrp_linear_cell_interp(fvals, ts, lats(ip), ier)
       iertot = iertot + ier

    enddo

  if(iertot /= 0) then
     write(iertot_str, '(i10)') iertot
     ermesg = 'drifters_positions2geolonlat: ERROR '//trim(iertot_str)// &
          & ' interpolation errors (domain out of bounds?)'
  endif

  end subroutine drifters_positions2lonlat

  !============================================================================
  !> @brief Print Runge-Kutta check sums. Useful for debugging only.
  subroutine drifters_print_checksums(self, pe, ermesg)

    type(drifters_type) :: self !< Opaque handle.
    integer, intent(in), optional :: pe !< Processor element.
    character(len=*), intent(out) :: ermesg !< Error message (if any).

    integer, parameter :: i8 = selected_int_kind(13)
    integer(i8) :: mold, chksum_pos, chksum_k1, chksum_k2, chksum_k3, chksum_k4
    integer(i8) :: chksum_tot
    integer nd, np, me

    ermesg = ''

    if(.not. present(pe)) then
       me = _MPP_PE
    else
       me = pe
    endif

    if(me == _MPP_PE) then

       nd = self%core%nd
       np = self%core%np
       chksum_pos = transfer(sum(sum(self%core%positions(1:nd,1:np),1)), mold)
       chksum_k1  = transfer(sum(sum(self%rk4_k1(1:nd,1:np),1)), mold)
       chksum_k2  = transfer(sum(sum(self%rk4_k2(1:nd,1:np),1)), mold)
       chksum_k3  = transfer(sum(sum(self%rk4_k3(1:nd,1:np),1)), mold)
       chksum_k4  = transfer(sum(sum(self%rk4_k4(1:nd,1:np),1)), mold)
       chksum_tot = chksum_pos + chksum_k1 + chksum_k2 + chksum_k3 +chksum_k4

       print *,'==============drifters checksums=========================='
       print '(a,i25,a,i6,a,e15.7)','==positions: ', chksum_pos,  ' PE=', me, ' time = ', self%time
       print '(a,i25,a,i6,a,e15.7)','==k1       : ', chksum_k1,   ' PE=', me, ' time = ', self%time
       print '(a,i25,a,i6,a,e15.7)','==k2       : ', chksum_k2,   ' PE=', me, ' time = ', self%time
       print '(a,i25,a,i6,a,e15.7)','==k3       : ', chksum_k3,   ' PE=', me, ' time = ', self%time
       print '(a,i25,a,i6,a,e15.7)','==k4       : ', chksum_k4,   ' PE=', me, ' time = ', self%time
       print '(a,i25,a,i6,a,e15.7)','==total    : ', chksum_tot,  ' PE=', me, ' time = ', self%time

    endif

  end subroutine drifters_print_checksums

  subroutine drifters_reset_rk4(self, ermesg)
    type(drifters_type) :: self
    character(len=*), intent(out) :: ermesg

    integer ier, nd

    ermesg = ''

    if(size(self%rk4_k1, 2) < self%core%np) then
       deallocate(self%rk4_k1, stat=ier)
       allocate(self%rk4_k1(self%core%nd, self%core%npdim))
       self%rk4_k1 = 0
    endif
    if(size(self%rk4_k2, 2) < self%core%np) then
       deallocate(self%rk4_k2, stat=ier)
       allocate(self%rk4_k2(self%core%nd, self%core%npdim))
       self%rk4_k2 = 0
    endif
    if(size(self%rk4_k3, 2) < self%core%np) then
       deallocate(self%rk4_k3, stat=ier)
       allocate(self%rk4_k3(self%core%nd, self%core%npdim))
       self%rk4_k3 = 0
    endif
    if(size(self%rk4_k4, 2) < self%core%np) then
       deallocate(self%rk4_k4, stat=ier)
       allocate(self%rk4_k4(self%core%nd, self%core%npdim))
       self%rk4_k4 = 0
    endif

    if(size(self%remove) < self%core%np) then
       deallocate(self%remove, stat=ier)
       allocate(self%remove(self%core%npdim))
       self%remove = .FALSE.
    endif

    if(size(self%temp_pos, 2) < self%core%np) then
       deallocate(self%temp_pos, stat=ier)
       nd = size(self%input%velocity_names)
       allocate(self%temp_pos(nd, self%core%npdim))
       self%temp_pos = -huge(1.)
    endif

  end subroutine drifters_reset_rk4

end module drifters_mod
!> @}
! close documentation grouping
