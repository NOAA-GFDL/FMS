# Drifters

Alexander.Pletzer@noaa.gov

## Overview:

1. [What are drifters?](#1-what-are-drifters)
2. [The drifters' module](#2-the-drifters-module)
3. [What the drifters' module can and cannot do?](#3-what-the-drifters-module-can-and-cannot-do)
4. [How to invoke the drifters' module in your code](#4-how-to-invoke-the-drifters-module-in-your-code)
5. [Pre- and post-processing](#5-pre-and-post-processing)
6. [Loose ends](#6-loose-ends)

## 1. What are drifters?

Drifters are idealized point particles with positions that evolve in time according
to a prescribed velocity field, starting from some initial conditions. Drifters have
no mass, no energy, no size, and no friction and therefore have no impact on the
dynamics of the underlying system. The only feature that distinguishes a drifter
from another is its trajectory. This makes drifters ideal for tracking pollution
clouds and probing fields (e.g. temperature, salinity) along ocean currents.
Drifters can mimic real experiments such as the Argo floats
<http://www.metoffice.com/research/ocean/argo/ukfloats.html>.

## 2. The drifters' module

The drifters' module (drifters_mod) is a Fortran module built on top of low level
kernels. The code was designed to run in parallel within a coupled FMS model.

The general layout is:

##### Top layer:

* drifters.F90: exposes the API

##### Some low level kernels:

* cloud_interpolator.F90: (linear) interpolation of a structured field at random
positions
* drifters_comm.F90: updates positions across processor (PE) domains
* drifters_core.F90: handles the mechanics for adding/removing drifters
* drifters_input.F90: imports initial positions from a NetCDF file
* drifters_io.F90: saves data for postprocessing and restart

## 3. What the drifters' module can and cannot do?

The drifters' module will handle all the communication necessary to keep track of
positions within PE domains. This will involve adding/removing drifters when they
enter/leave a PE domain. As the number of drifters can vary greatly both between PE
domains and in the course of a simulation, the drifters' module will also manage
dynamically the memory for you.

The drifters' code was written to run efficiently by avoiding unnecessary data
copies. There is no software limit as to the maximum number of drifters and their
distribution. The drifters' module is ideally suited for long runs where the amount
of data required to advect particles is large, making it impractical to use a
post-processing approach. Such a situation occurs, in particular, when the velocity
fields must be saved at high sampling rate in order to compute accurately the
drifters' positions.

Although a Runge-Kutta method is provided to push drifters, users may want to
implement their own time-integration scheme in some cases. It is moreover possible to
combine time integration along some axes while prescribing the evolution in the
third axis, as would be the case for robotic floats that move up and down in
predetermined way, for instance.

There are a number of basic assumptions which could make the drifters' module
ill-suited for some tasks. First and foremost, it is assumed that the motion of
drifters is not erratic but follows deterministic trajectories. Furthermore,
drifters should not cross both compute and data domain boundaries within less
than a time step. This limitation is imposed by the Runge-Kutta integration
scheme, which must be able to complete, within a time step, a trajectory
calculation that starts inside the compute domain. Therefore, the drifters,
as they are presently modeled, are unlikely to work for very fast objects.
This constraint also puts a upper limit to the domain decomposition, although
it can often be remedied by increasing the number of ghost nodes.

Another fundamental assumption is that the (e.g. velocity) fields are structured,
on a per PE domain basis. There is no support for locally nested or unstructured
meshes. Meshes need not be smooth and continuous across PE domains, however.

If the higher level module drifters_mod cannot do what you need, the chances are
the kernels will be able to. For instance, the upper level assumes, at present,
two-dimensional (x,y) positions, and a two-dimension domain decomposition. The
cloud_interpolator kernel on the other hand makes no such assumption; it can be
used to interpolate structured data in arbitrary number of dimensions and with no
dependence on FMS code.


## 4. How to invoke the drifters' module in your code

In order to implement drifters in your code, you will need to get access to the
velocity fields at all times. You will also need to settle on a definition of a
position, which could be the longitude, latitude, elevation/depth triplet. The
drifters' code does not care what your choice of position is. Typically, you will
want, however, to choose positions that are aligned to the velocity axes, and to
your coordinate system. It is perfectly fine to use a warped coordinate system
(e.g. tripolar), which avoids the problem of geometric singularity (pole).

Naturally, it is important for the velocity field's units to be consistent with
both positions and time. Often, some metric coefficients are required for this
transformation. The drifters' module does not have any support for specific meshes
per se. Each component of the velocity can sit on different nodes, as in the case
of Arakawa C and D meshes.

You will also need to know what the boundaries of your "compute" and "data" domains are,
when running in parallel. Your compute domain, expressed in terms of `xcmin`, `xcmax`,
`ycmin` and `ycmax` should densely cover the entire domain; that is
```Fortran
xcmin[pe_east]  = xcmax[pe]
ycmin[pe_north] = ycmax[pe]
etc.
```

The data domain should be larger than the compute domain. The larger the velocities
or the larger the time steps, the more ghost nodes are required. You should also
pay attention to the fact that *all* components of your velocity need to fall inside
the data domain. If you're using a staggered grid (Arakawa's B, C grids) you should
reduce the max data domain by `dx`, resp. `dy`, in every direction.

The implementation of drifters involves typically the following steps, which will
be assumed to run on different PEs:

```Fortran
  use drifters_mod
  ...

  ! declare drifters object
  type(drifters_type) :: drfts ! drifters' object

  ! constructor
  call drifters_new(drfts, &
       & input_file ='INPUT/drifters_inp.nc'  , &
       & output_file='drifters_out.nc', &
       & ermesg=ermesg)

  ! set the initial time and dt
  drfts%time = t0
  drfts%dt   = dt


  ! set the PE domain boundaries. Xmin_comp/ymin_comp, xmax_comp/ymax_comp
  ! refer to the "compute" domain, which should cover densily the domain:
  ! xcmax[pe] = xcmin[pe_east]
  ! ycmax[pe] = ycmin[pe_north]
  ! Xmin_data/ymin_data, xmax_data/ymax_data refer to the "data" domain, which
  ! should be larger than the compute domain and therefore overlap:
  ! xdmax[pe] > xdmin[pe_east]
  ! ydmax[pe] > ydmin[pe_north]
  ! Particles in the overlap regions are tracked by several PEs.

  call drifters_set_domain(drfts, &
       & xmin_comp=xcmin, xmax_comp=xcmax, &
       & ymin_comp=ycmin, ymax_comp=ycmax, &
       & xmin_data=xdmin, xmax_data=xdmax, &
       & ymin_data=ydmin, ymax_data=ydmax, &
       & xmin_glob=xmin , xmax_glob=xmax,  & ! this will set periodicity in x
       & ermesg=ermesg)

  ! set neighboring PEs [domain2d is of type(domain2d)]

  call drifters_set_pe_neighbors(drfts, domain=domain2d, ermesg=ermesg)

  ! set the velocities axes. Each velocity can have different axes.

  call drifters_set_v_axes(drfts, component='u', &
       & x=x_u, y=y_u, z=z_u, ermesg=ermesg)

  call drifters_set_v_axes(drfts, component='v', &
       & x=x_v, y=y_v, z=z_v, ermesg=ermesg)

  ...

  ! Distribute the drifters across PEs
  call drifters_distribute(drfts, ermesg)

  ! Push the drifters. u_comp, v_comp etc are provided by the host code

  call drifters_push(drfts, u=u_comp, v=v_comp, w=w_comp, ermesg=ermesg)

  do while (Time < Time_end)

  ! Update the velocities
  ....

  ! Push the drifters
  call drifters_push(drfts, u=u_comp, v=v_comp, w=w_comp, ermesg=ermesg)

  ! check if RK4 integration is complete
  if(drfts%rk4_step==1) then

       ! interpolate fields

       call drifters_set_field(drfts, index_field=1, x=x1, y=y1, z=z1, &
       &    data=temp, ermesg=ermesg)

       .... ! set other fields here

       ! save data

       call drifters_save(drfts, ermesg=ermesg)

  endif

  enddo

  ...
  ! write restart file, optionally with lon/lat data coordinates

  call drifters_write_restart(drfts, filename='RESTART/drifters_inp.nc', &
       & x1=x_u, y1=y_v, geolon1=grid%geolonq, &
       & x2=x_u, y2=y_v, geolat2=grid%geolatq, &
       & ermesg=ermesg)  


  ! destroy

  call drifters_del(drfts, ermesg=ermesg)
```

## 5. Pre- and post-processing

A typical NetCDF input file contains:

```
netcdf drifters_inp {
dimensions:
  nd = 2 ; // number of dimensions (2 or 3)
  np = 5 ; // number of particles
variables:
  double positions(np, nd) ;
  positions:names = "i j" ;
  positions:units = "- -" ;
  int ids(np) ;

// global attributes:
  :velocity_names = "u v" ;
  :field_names = "lon lat temp salt" ;
  :field_units = "deg_N deg_E Celsius ppu" ;
  :time_units = "seconds" ;
  :title = "example of input data for drifters" ;
data:

 positions =
  240, 100,
  240, 110,
  240, 120,
  240, 130,
  240, 140 ;

 ids = 1, 2, 3, 4, 5 ; // must range from 1 to np, in any order
}
```

In this case, the 2-d initial conditions consist of 5 particles at positions
(240,100), (240, 110), .... Here, the positions are defined as i, j global indices.
The fields to interpolate are the longitudes, latitudes, temperature and salinity.

The output files have names `<output>.nc.<PE>` where PE ranges 0...number of domains-1.
Each PE writes to an independent NetCDF file. These files contain a history of
positions and field values:

```
netcdf drifters_out.nc {
dimensions:
  it_id = UNLIMITED ; // (10 currently)
  nf = 4 ;
  nd = 2 ;
variables:
  int index_time(it_id) ;
  double time(it_id) ;
  int ids(it_id) ;
  double positions(it_id, nd) ;
  positions:name_1 = "i" ;
  positions:name_2 = "j" ;
  positions:unit_1 = "-" ;
  positions:unit_2 = "-" ;
  double fields(it_id, nf) ;
  fields:name_1 = "lon" ;
  fields:name_2 = "lat" ;
  fields:name_3 = "temp" ;
  fields:name_4 = "salt" ;
  fields:unit_1 = "deg_N" ;
  fields:unit_2 = "deg_E" ;
  fields:unit_3 = "Celsius" ;
  fields:unit_4 = "ppu" ;

// global attributes:
  :time_units = "seconds" ;
data:

 index_time = 0, 0, 1, 1, 2, 2, 3, 3, 4, 4 ;

 time = 17280, 17280, 34560, 34560, 51840, 51840, 69120, 69120, 86400, 86400 ;

 ids = 1, 2, 1, 2, 1, 2, 1, 2, 1, 2 ;

 positions =
  239.992427824695, 99.9982201799336,
  239.992780388996, 110.000731652201,
...
 fields =
  -52.5075721753049, 11.5056520209176, 26.5841099585423, 35.4729937970157,
...
}
```

A tool, drifters_combine, has been written to combine all the NetCDF files into a
single file, which contains all the drifters data and is suitable for post-processing,
including visualization using Ncvtk 1.7 or later.

## 6. Loose ends

* What to do with fields that are outside the valid_range (e.g. on continents in the
case of floats). Use piece-wise linear interpolation or fill with missing values?
* How to make the drifters Mosaic compliant?
* Allow drifters to turn on/off at specific time? A predetermined life cycle?
* Extract initial positions from databases.
* How best to visualize drifters?

