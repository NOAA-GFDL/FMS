!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
! Check physical interpolation independently of source storage order.
program test_coordinate_order
  use fms_mod, only: fms_init, check_nml_error
  use mpp_mod, only: mpp_exit, input_nml_file, mpp_error, FATAL
  use platform_mod, only: r4_kind, r8_kind
  use horiz_interp_mod, only: horiz_interp_new, horiz_interp, horiz_interp_del, horiz_interp_type
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: k=TEST_FMS_KIND_, nx=8, ny=5, mx=16, my=7
  real(k), parameter :: period=2*acos(-1.0_k), tolerance=128*epsilon(1.0_k)
  character(8), parameter :: methods(*)=[character(8)::'bicubic','bilinear']
  character(9), parameter :: fields(*)=[character(9)::'constant','affine','periodic']
  character(10), parameter :: spacings(*)=[character(10)::'uniform','nonuniform']
  character(4), parameter :: interfaces(*)=[character(4)::'1d','2d']
  character(9), parameter :: orientations(*)=[character(9)::'ascending','longitude','latitude','both']

  type :: source_grid
    real(k) :: longitude(nx), latitude(ny), values(nx,ny)
  end type

  type :: destination_grid
    real(k) :: longitude(mx), latitude(my)
    real(k) :: longitude_2d(mx,my), latitude_2d(mx,my)
  end type

  character(32) :: scenario='cyclic', interface='1d', longitude_order='ascending'
  integer :: method, interface_index, failures=0, io, ierr
  namelist /test_coordinate_order_nml/ scenario, interface, longitude_order

  call fms_init
  read (input_nml_file, test_coordinate_order_nml, iostat=io)
  ierr = check_nml_error(io, 'test_coordinate_order_nml')

  if (interface/='1d' .and. interface/='2d') &
       call mpp_error(FATAL, 'test_coordinate_order_nml: invalid interface '//trim(interface)//' (expected 1d or 2d)')
  if (longitude_order/='ascending' .and. longitude_order/='longitude') &
       call mpp_error(FATAL, 'test_coordinate_order_nml: invalid longitude_order '//trim(longitude_order)// &
       ' (expected ascending or longitude)')

  select case (scenario)
  case ('cyclic','noncyclic')
    do interface_index=1,size(interfaces)
      do method=1,size(methods)
        call test_coordinate_reversal(methods(method), interfaces(interface_index), scenario=='cyclic')
      enddo
    enddo
  case ('compatibility')
    call test_bicubic_legacy_normalization(interface)
  case ('outside')
    call test_noncyclic_longitude_rejection(interface, longitude_order)
  case default
    call mpp_error(FATAL, 'test_coordinate_order_nml: invalid scenario '//trim(scenario)// &
         ' (expected cyclic, noncyclic, compatibility or outside)')
  end select

  call mpp_exit
  if (failures>0) stop 1

contains

  ! Run one method/interface pair on both spacings and all three fields.
  subroutine test_coordinate_reversal(method, interface, periodic, legacy)
    character(*), intent(in) :: method, interface
    logical, intent(in) :: periodic
    logical, optional, intent(in) :: legacy
    type(source_grid) :: source
    type(destination_grid) :: destination
    logical :: legacy_normalization
    integer :: spacing, field

    legacy_normalization=.false.
    if (present(legacy)) legacy_normalization=legacy
    do spacing=1,size(spacings)
      call make_source_grid(spacings(spacing), source)
      call make_destinations(source, interface, periodic, legacy_normalization, destination)
      do field=1,size(fields)
        call initialize_field(fields(field), source)
        call check_source_orderings(source, destination, method, interface, spacings(spacing), &
             fields(field), periodic, legacy_normalization)
      enddo
    enddo
  end subroutine

  ! Explicit false still wraps bicubic longitude, including 1-2*pi.
  subroutine test_bicubic_legacy_normalization(interface)
    character(*), intent(in) :: interface

    call test_coordinate_reversal('bicubic', interface, .true., legacy=.true.)
  end subroutine

  ! The shell driver verifies the specific fatal message from this call.
  subroutine test_noncyclic_longitude_rejection(interface, order)
    character(*), intent(in) :: interface, order
    type(source_grid) :: source, reordered
    type(destination_grid) :: destination
    real(k) :: result(mx,my)

    call make_source_grid('uniform', source)
    call initialize_field('constant', source)
    call make_destinations(source, interface, .false., .false., destination)
    destination%longitude(1)=source%longitude(1)-0.1_k
    destination%longitude_2d(1,:)=destination%longitude(1)
    call reorder_source(source, order, reordered)
    call interpolate(reordered, destination, 'bilinear', interface, .false., result)
  end subroutine

  subroutine make_source_grid(spacing, source)
    character(*), intent(in) :: spacing
    type(source_grid), intent(out) :: source
    real(k) :: fraction
    integer :: i

    do i=1,nx
      fraction=real(i-1,k)/nx
      if (spacing=='nonuniform') fraction=fraction**1.2_k
      source%longitude(i)=0.25_k+period*fraction
    enddo
    do i=1,ny
      fraction=real(i-1,k)/(ny-1)
      if (spacing=='nonuniform') fraction=fraction**1.3_k
      source%latitude(i)=-1.0_k+2*fraction
    enddo
  end subroutine

  subroutine initialize_field(field, source)
    character(*), intent(in) :: field
    type(source_grid), intent(inout) :: source
    real(k) :: longitude, latitude
    integer :: i,j

    do j=1,ny
      latitude=source%latitude(j)
      do i=1,nx
        longitude=source%longitude(i)
        select case (field)
        case ('constant')
          source%values(i,j)=2.0_k
        case ('affine')
          source%values(i,j)=0.5_k*longitude+0.75_k*latitude
        case ('periodic')
          source%values(i,j)=sin(longitude)*(1+latitude**2)+0.3_k*cos(2*longitude)+0.2_k*latitude
        end select
      enddo
    enddo
  end subroutine

  subroutine reorder_source(source, orientation, reordered)
    type(source_grid), intent(in) :: source
    character(*), intent(in) :: orientation
    type(source_grid), intent(out) :: reordered

    reordered=source
    if (orientation=='longitude' .or. orientation=='both') then
      reordered%longitude=source%longitude(nx:1:-1)
      reordered%values=source%values(nx:1:-1,:)
    endif
    if (orientation=='latitude' .or. orientation=='both') then
      reordered%latitude=source%latitude(ny:1:-1)
      reordered%values=reordered%values(:,ny:1:-1)
    endif
  end subroutine

  subroutine make_destinations(source, interface, periodic, legacy, destination)
    type(source_grid), intent(in) :: source
    character(*), intent(in) :: interface
    logical, intent(in) :: periodic, legacy
    type(destination_grid), intent(out) :: destination
    integer :: i,j

    ! First three points are source nodes. Columns 10–16 cover the seam and wrapped nodes.
    destination%longitude=[source%longitude(1),source%longitude(2),source%longitude(nx), &
         0.8_k,1.4_k,2.0_k,3.0_k,4.0_k,5.0_k,-0.125_k,0.125_k,period-0.125_k, &
         2*period-0.125_k,-period-0.125_k,source%longitude(1)+period,source%longitude(nx)-period]
    if (.not.periodic) destination%longitude(10:16)=[source%longitude(1)+0.1_k, &
         source%longitude(nx)-0.1_k,source%longitude(4:6),source%longitude(1),source%longitude(nx)]
    if (legacy) destination%longitude(4:5)=[1.0_k-period,1.0_k]

    ! Last two rows exercise clamping beyond each physical latitude endpoint.
    destination%latitude=[source%latitude(1),source%latitude(2),source%latitude(ny), &
         -0.6_k,0.6_k,-1.2_k,1.2_k]
    do j=1,my
      do i=1,mx
        destination%longitude_2d(i,j)=destination%longitude(i)
        destination%latitude_2d(i,j)=destination%latitude(j)
        if (interface=='2d' .and. i>3 .and. i<15 .and. j>3 .and. j<6) then
          destination%longitude_2d(i,j)=destination%longitude(i)+0.015_k*j
          destination%latitude_2d(i,j)=destination%latitude(j)+0.01_k*i
        endif
      enddo
    enddo
  end subroutine

  subroutine interpolate(source, destination, method, interface, src_modulo, result)
    type(source_grid), intent(in) :: source
    type(destination_grid), intent(in) :: destination
    character(*), intent(in) :: method, interface
    logical, intent(in) :: src_modulo
    real(k), intent(out) :: result(mx,my)
    type(horiz_interp_type) :: interpolation

    select case (interface)
    case ('1d')
      call horiz_interp_new(interpolation, source%longitude, source%latitude, &
           destination%longitude, destination%latitude, interp_method=trim(method), &
           src_modulo=src_modulo, grid_at_center=.true.)
    case ('2d')
      call horiz_interp_new(interpolation, source%longitude, source%latitude, &
           destination%longitude_2d, destination%latitude_2d, interp_method=trim(method), &
           src_modulo=src_modulo, grid_at_center=.true.)
    case default
      write(*,*) 'Unknown destination interface ', trim(interface)
      stop 1
    end select
    if (method=='bicubic') call check_periodic_endpoint_weights(interpolation, source, interface)
    call horiz_interp(interpolation, source%values, result)
    call horiz_interp_del(interpolation)
  end subroutine

  subroutine check_source_orderings(source, destination, method, interface, spacing, field, periodic, legacy)
    type(source_grid), intent(in) :: source
    type(destination_grid), intent(in) :: destination
    character(*), intent(in) :: method, interface, spacing, field
    logical, intent(in) :: periodic, legacy
    type(source_grid) :: reordered
    real(k) :: result(mx,my), ascending(mx,my), reference(mx,my), linear_reference(mx,my), analytic(mx,my)
    real(k) :: longitude, latitude, reference_limit
    character(160) :: context
    integer :: i,j,orientation

    do j=1,my
      do i=1,mx
        longitude=destination%longitude_2d(i,j)
        latitude=destination%latitude_2d(i,j)
        linear_reference(i,j)=bilinear_reference(source, longitude, latitude, periodic)
        reference(i,j)=linear_reference(i,j)
        if (method=='bicubic') reference(i,j)=bicubic_hermite_reference(source, longitude, latitude, periodic)
        analytic(i,j)=0.5_k*longitude+0.75_k*max(source%latitude(1),min(source%latitude(ny),latitude))
        if (method=='bicubic' .and. .not.periodic) &
             analytic(i,j)=affine_periodic_reference(source, longitude, latitude)
      enddo
    enddo
    ! Preserve the original reference scale for both algorithms.
    reference_limit=tolerance*max(1.0_k,maxval(abs(linear_reference)))

    do orientation=1,size(orientations)
      call reorder_source(source, orientations(orientation), reordered)
      call interpolate(reordered, destination, method, interface, periodic.and..not.legacy, result)
      write(context,'(a,i0,7(1x,a))') 'kind=',k,trim(method),interface,spacing,field,trim(scenario), &
           'reversal',trim(orientations(orientation))
      if (orientations(orientation)=='ascending') ascending=result
      call check_close('source reversal', result, ascending, tolerance*max(1.0_k,maxval(abs(ascending))), context)
      call check_close('independent reference', result, reference, reference_limit, context)
      call check_close('latitude clamping', result(:,6:7), reference(:,6:7), reference_limit, context)
      if (field=='constant') call check_close('constant field', result, 0*reference+2.0_k, reference_limit, context)
      if (field=='affine' .and. .not.periodic) &
           call check_close('affine samples with method endpoint slopes', result, analytic, reference_limit, context)
      if (periodic) call check_equivalent_longitudes(result, interface, legacy, reference_limit, context)
    enddo
  end subroutine

  ! On the perturbed 2D grid, rows 1–3 still have identical physical latitudes.
  subroutine check_equivalent_longitudes(result, interface, legacy, limit, context)
    real(k), intent(in) :: result(mx,my), limit
    character(*), intent(in) :: interface, context
    logical, intent(in) :: legacy
    integer :: last_row

    last_row=my
    if (interface=='2d') last_row=3
    call check_close('seam longitude plus one turn', result(10:10,:last_row), result(12:12,:last_row), limit, context)
    call check_close('seam longitude plus two turns', result(10:10,:last_row), result(13:13,:last_row), limit, context)
    call check_close('seam longitude minus one turn', result(10:10,:last_row), result(14:14,:last_row), limit, context)
    call check_close('wrapped first node', result(1:1,:last_row), result(15:15,:last_row), limit, context)
    call check_close('wrapped last node', result(3:3,:last_row), result(16:16,:last_row), limit, context)
    if (legacy) call check_close('explicit false at 1-2*pi', result(4:4,1:3), result(5:5,1:3), limit, context)
  end subroutine

  subroutine check_close(assertion, actual, expected, limit, context)
    character(*), intent(in) :: assertion, context
    real(k), intent(in) :: actual(:,:), expected(:,:), limit
    real(k) :: error
    integer :: point(2)

    if (.not.all(ieee_is_finite(actual))) then
      write(*,'(a)') trim(context)//' '//assertion//' returned nonfinite values'
      failures=failures+1
      return
    endif
    error=maxval(abs(actual-expected))
    if (error<=limit) return
    point=maxloc(abs(actual-expected))
    write(*,'(a,2(a,es14.6),a,2i4)') trim(context)//' '//assertion, &
         ' error=',error,' limit=',limit,' index=',point
    failures=failures+1
  end subroutine

  ! References use the strictly increasing physical source grid. Wrap only periodic
  ! longitudes and clamp latitude before finding the enclosing cell.
  subroutine reference_cell(source, longitude, latitude, periodic, west, east, south, north, dx, tx, ty)
    type(source_grid), intent(in) :: source
    real(k), intent(in) :: longitude, latitude
    logical, intent(in) :: periodic
    integer, intent(out) :: west,east,south,north
    real(k), intent(out) :: dx,tx,ty
    real(k) :: physical_lon, physical_lat, east_lon

    physical_lon=longitude
    if (periodic) physical_lon=source%longitude(1)+modulo(longitude-source%longitude(1),period)
    physical_lat=max(source%latitude(1),min(source%latitude(ny),latitude))

    if (physical_lon>source%longitude(nx)) then
      west=nx
      east=1
      east_lon=source%longitude(1)+period
    else
      west=lower_neighbor(source%longitude, physical_lon)
      east=west+1
      east_lon=source%longitude(east)
    endif
    south=lower_neighbor(source%latitude, physical_lat)
    north=south+1
    dx=east_lon-source%longitude(west)
    tx=(physical_lon-source%longitude(west))/dx
    ty=(physical_lat-source%latitude(south))/(source%latitude(north)-source%latitude(south))
  end subroutine

  integer function lower_neighbor(axis, coordinate) result(lower)
    real(k), intent(in) :: axis(:), coordinate

    do lower=1,size(axis)-2
      if (coordinate<=axis(lower+1)) return
    enddo
  end function

  real(k) function bilinear_reference(source, longitude, latitude, periodic) result(value)
    type(source_grid), intent(in) :: source
    real(k), intent(in) :: longitude, latitude
    logical, intent(in) :: periodic
    integer :: west,east,south,north
    real(k) :: dx,tx,ty, southern_value,northern_value

    call reference_cell(source, longitude, latitude, periodic, west,east,south,north,dx,tx,ty)
    southern_value=(1-tx)*source%values(west,south)+tx*source%values(east,south)
    northern_value=(1-tx)*source%values(west,north)+tx*source%values(east,north)
    value=(1-ty)*southern_value+ty*northern_value
  end function

  ! Two successive 1D Hermite evaluations test implementation consistency.
  ! Longitude slopes are periodic centered differences, including endpoints.
  ! Latitude slopes retain centered interior and one-sided endpoint differences,
  ! so agreement here does not measure exact analytical accuracy.
  real(k) function bicubic_hermite_reference(source, longitude, latitude, periodic) result(value)
    type(source_grid), intent(in) :: source
    real(k), intent(in) :: longitude, latitude
    logical, intent(in) :: periodic
    integer :: west,east,south,north,row
    real(k) :: dx,tx,ty, interpolated_rows(ny)

    call reference_cell(source, longitude, latitude, periodic, west,east,south,north,dx,tx,ty)
    do row=1,ny
      interpolated_rows(row)=hermite_value(source%values(west,row),source%values(east,row), &
           periodic_longitude_slope(source%longitude,source%values(:,row),west), &
           periodic_longitude_slope(source%longitude,source%values(:,row),east),dx,tx)
    enddo
    value=hermite_value(interpolated_rows(south),interpolated_rows(north), &
         slope(source%latitude,interpolated_rows,south),slope(source%latitude,interpolated_rows,north), &
         source%latitude(north)-source%latitude(south),ty)
  end function

  ! The reference grid is stored in increasing physical longitude order.
  real(k) function periodic_longitude_slope(axis, values, node) result(derivative)
    real(k), intent(in) :: axis(:), values(:)
    integer, intent(in) :: node
    integer :: count

    count=size(axis)
    if (node==1) then
      derivative=(values(2)-values(count))/(axis(2)-(axis(count)-period))
    else if (node==count) then
      derivative=(values(1)-values(count-1))/((axis(1)+period)-axis(count-1))
    else
      derivative=slope(axis,values,node)
    endif
  end function

  ! Affine longitude samples are discontinuous under periodic extension. For
  ! f=0.5*x+0.75*y, endpoint slopes differ from 0.5 by -0.5*period/span.
  ! Apply those two known Hermite corrections to the exact affine value.
  ! This checks every in-range sample, including the endpoint-adjacent cells.
  real(k) function affine_periodic_reference(source, longitude, latitude) result(value)
    type(source_grid), intent(in) :: source
    real(k), intent(in) :: longitude, latitude
    integer :: west, east
    real(k) :: width, fraction, slope_error

    west=lower_neighbor(source%longitude,longitude)
    east=west+1
    width=source%longitude(east)-source%longitude(west)
    fraction=(longitude-source%longitude(west))/width
    value=0.5_k*longitude+0.75_k*max(source%latitude(1),min(source%latitude(ny),latitude))
    if (west==1) then
      slope_error=-0.5_k*period/(source%longitude(2)-source%longitude(nx)+period)
      value=value+width*(fraction**3-2*fraction**2+fraction)*slope_error
    endif
    if (east==nx) then
      slope_error=-0.5_k*period/(source%longitude(1)+period-source%longitude(nx-1))
      value=value+width*(fraction**3-fraction**2)*slope_error
    endif
  end function

  ! Inspect both physical endpoints independently of storage direction. The
  ! interpolation reference separately checks the matching value numerators.
  subroutine check_periodic_endpoint_weights(interpolation, source, interface)
    type(horiz_interp_type), intent(in) :: interpolation
    type(source_grid), intent(in) :: source
    character(*), intent(in) :: interface
    integer :: west, east, step, endpoints(2), j
    real(k) :: span(2), latitude_span, actual(2,ny), expected(2,ny)
    character(160) :: context

    west=minloc(source%longitude,dim=1)
    east=maxloc(source%longitude,dim=1)
    step=1
    if (west==nx) step=-1
    endpoints=[west,east]
    span(1)=step*(source%longitude(west+step)-(source%longitude(east)-period))
    span(2)=step*((source%longitude(west)+period)-source%longitude(east-step))
    write(context,'(a,i0,2(1x,a),2(a,l1))') 'kind=',k,'bicubic',interface, &
         ' decreasing longitude=',step<0,' decreasing latitude=',source%latitude(1)>source%latitude(ny)
    do j=1,ny
      expected(:,j)=1/span
    enddo
    if (interpolation%horizInterpReals4_type%is_allocated) then
      actual=interpolation%horizInterpReals4_type%wti(endpoints,:,1)
    else
      actual=interpolation%horizInterpReals8_type%wti(endpoints,:,1)
    endif
    call check_close('periodic endpoint longitude weights',actual,expected, &
         tolerance*max(1.0_k,maxval(abs(expected))),context)
    do j=1,ny
      latitude_span=source%latitude(min(j+1,ny))-source%latitude(max(j-1,1))
      expected(:,j)=1/(span*latitude_span)
    enddo
    if (interpolation%horizInterpReals4_type%is_allocated) then
      actual=interpolation%horizInterpReals4_type%wti(endpoints,:,3)
    else
      actual=interpolation%horizInterpReals8_type%wti(endpoints,:,3)
    endif
    call check_close('periodic endpoint mixed weights',actual,expected, &
         tolerance*max(1.0_k,maxval(abs(expected))),context)
  end subroutine

  real(k) function slope(axis, values, node) result(derivative)
    real(k), intent(in) :: axis(:), values(:)
    integer, intent(in) :: node
    integer :: left,right

    left=max(node-1,1)
    right=min(node+1,size(axis))
    derivative=(values(right)-values(left))/(axis(right)-axis(left))
  end function

  real(k) function hermite_value(left, right, left_slope, right_slope, width, fraction) result(value)
    real(k), intent(in) :: left,right,left_slope,right_slope,width,fraction
    real(k) :: t, left_basis,right_basis,left_slope_basis,right_slope_basis

    t=fraction
    left_basis=2*t**3-3*t**2+1
    right_basis=-2*t**3+3*t**2
    left_slope_basis=t**3-2*t**2+t
    right_slope_basis=t**3-t**2
    value=left_basis*left+right_basis*right+width*(left_slope_basis*left_slope+right_slope_basis*right_slope)
  end function
end program
