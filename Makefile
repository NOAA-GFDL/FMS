#CC = cc
#CPPDEFS = -DINTERNAL_FILE_NML -g -Duse_libMPI -Duse_netCDF -D_FILE_VERSION="foo"
#CPPFLAGS = -D__IFC $(shell nc-config --cflags)
#CFLAGS := -msse2 -sox -traceback -O0 -g -ftrapuv -qopenmp
OTHERFLAGS =
OTHER_CFLAGS =
#FC = ftn
#FPPFLAGS = -fpp -Wp,-w $(shell nf-config --fflags)
#FFLAGS := -msse2 -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -i4 -r8 -nowarn -sox -traceback -g -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fpe0 -ftrapuv -qopenmp
OTHER_FFLAGS =

OBJ = mpp_data.o mpp_domains.o monin_obukhov_kernel.o atmos_ocean_fluxes.o \
      fms_io_utils.o affinity.o time_interp.o memutils.o astronomy.o fft.o monin_obukhov.o drifters_comm.o \
      tridiagonal.o station_data.o diag_integral.o drifters_input.o time_manager.o horiz_interp_type.o \
      mpp_pset.o horiz_interp_conserve.o amip_interp.o quicksort.o block_control.o \
      xbt_drop_rate_adjust.o horiz_interp_spherical.o fms2_io.o fm_util.o \
      netcdf_io.o gradient_c2l.o mpp.o topography.o fms.o column_diagnostics.o axis_utils.o \
      random_numbers.o mpp_parameter.o gradient.o diag_util.o gaussian_topog.o \
      drifters_core.o threadloc.o stock_constants.o interpolator.o mpp_utilities.o diag_axis.o \
      field_manager.o drifters.o fms_netcdf_domain_io.o cloud_interpolator.o legacy.o \
      fms_netcdf_unstructured_domain_io.o sat_vapor_pres.o mpp_memutils.o diag_grid.o \
      diag_output.o coupler_types.o diag_manifest.o get_cal_time.o \
      xgrid.o diag_data.o oda_types.o oda_core_ecda.o diag_table.o read_mosaic.o \
      sat_vapor_pres_k.o interp.o time_interp_external.o \
      ensemble_manager.o MersenneTwister.o horiz_interp.o data_override.o drifters_io.o mosaic_util.o \
      platform.o horiz_interp_bilinear.o memuse.o tracer_manager.o nsclock.o fft99.o \
      write_ocean_data.o create_xgrid.o grid.o constants.o mpp_efp.o oda_core.o fms_io.o \
      mpp_io.o mosaic.o horiz_interp_bicubic.o diag_manager.o

.DEFAULT:
	-echo $@ does not exist.

all: libfms.a

libfms.a: $(OBJ)
	$(AR) $(ARFLAGS) libfms.a $(OBJ)

MersenneTwister.o: random_numbers/MersenneTwister.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c random_numbers/MersenneTwister.F90

affinity.o: mpp/affinity.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c mpp/affinity.c

amip_interp.o: amip_interp/amip_interp.F90 include/file_version.h time_interp.o time_manager.o get_cal_time.o mpp_io.o horiz_interp.o fms.o fms_io.o constants.o platform.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude amip_interp/amip_interp.F90

astronomy.o: astronomy/astronomy.F90 include/file_version.h fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude astronomy/astronomy.F90

atmos_ocean_fluxes.o: coupler/atmos_ocean_fluxes.F90 include/file_version.h mpp.o fms.o coupler_types.o field_manager.o fm_util.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude coupler/atmos_ocean_fluxes.F90

axis_utils.o: axis_utils/axis_utils.F90 include/file_version.h mpp_io.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude axis_utils/axis_utils.F90

blackboxio.o: fms2_io/blackboxio.F90 fms_io_utils.o netcdf_io.o fms_netcdf_domain_io.o fms_netcdf_unstructured_domain_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c fms2_io/blackboxio.F90

block_control.o: block_control/block_control.F90 include/fms_platform.h mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude block_control/block_control.F90

cloud_interpolator.o: drifters/cloud_interpolator.F90 include/fms_platform.h include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude	drifters/cloud_interpolator.F90

column_diagnostics.o: column_diagnostics/column_diagnostics.F90 include/file_version.h mpp_io.o fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude column_diagnostics/column_diagnostics.F90

constants.o: constants/constants.F90 include/file_version.h platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude constants/constants.F90

coupler_types.o: coupler/coupler_types.F90 include/file_version.h fms.o fms_io.o time_manager.o diag_manager.o data_override.o mpp_domains.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude coupler/coupler_types.F90

create_xgrid.o: mosaic/create_xgrid.c mosaic/mosaic_util.h mosaic/create_xgrid.h mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -Imosaic mosaic/create_xgrid.c

data_override.o: data_override/data_override.F90 include/fms_platform.h include/file_version.h platform.o constants.o mpp_io.o mpp.o horiz_interp.o time_interp_external.o fms_io.o fms.o axis_utils.o mpp_domains.o time_manager.o diag_manager.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude data_override/data_override.F90

diag_axis.o: diag_manager/diag_axis.F90 include/fms_platform.h include/file_version.h mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude diag_manager/diag_axis.F90

diag_data.o: diag_manager/diag_data.F90 include/fms_platform.h include/file_version.h time_manager.o mpp_domains.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude diag_manager/diag_data.F90

diag_grid.o: diag_manager/diag_grid.F90 include/fms_platform.h include/file_version.h constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude diag_manager/diag_grid.F90

diag_integral.o: diag_integral/diag_integral.F90 include/fms_platform.h time_manager.o mpp.o fms.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude diag_integral/diag_integral.F90

diag_manager.o: diag_manager/diag_manager.F90 include/fms_platform.h include/file_version.h time_manager.o mpp_io.o mpp.o fms.o fms_io.o diag_axis.o diag_util.o diag_data.o diag_table.o diag_output.o diag_grid.o diag_manifest.o constants.o mpp_domains.o mpp_parameter.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude diag_manager/diag_manager.F90

diag_manifest.o: diag_manager/diag_manifest.F90 diag_data.o mpp.o fms.o fms_io.o time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c diag_manager/diag_manifest.F90

diag_output.o: diag_manager/diag_output.F90 include/fms_platform.h include/file_version.h mpp_io.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude diag_manager/diag_output.F90

diag_table.o: diag_manager/diag_table.F90 mpp_io.o mpp.o fms.o time_manager.o constants.o diag_data.o diag_util.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c diag_manager/diag_table.F90

diag_util.o: diag_manager/diag_util.F90 include/fms_platform.h include/file_version.h diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o fms_io.o mpp_domains.o time_manager.o mpp_io.o mpp.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude diag_manager/diag_util.F90

drifters.o: drifters/drifters.F90 drifters/fms_switches.h include/fms_platform.h include/file_version.h drifters/drifters_push.h drifters/drifters_set_field.h drifters/drifters_compute_k.h mpp.o mpp_domains.o drifters_core.o drifters_input.o drifters_io.o drifters_comm.o cloud_interpolator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Idrifters -Iinclude drifters/drifters.F90

drifters_comm.o: drifters/drifters_comm.F90 drifters/fms_switches.h include/fms_platform.h mpp.o mpp_domains.o drifters_core.o drifters_input.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Idrifters -Iinclude drifters/drifters_comm.F90

drifters_core.o: drifters/drifters_core.F90 include/fms_platform.h include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude drifters/drifters_core.F90

drifters_input.o: drifters/drifters_input.F90 include/fms_platform.h include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude drifters/drifters_input.F90

drifters_io.o: drifters/drifters_io.F90 include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude drifters/drifters_io.F90

ensemble_manager.o: coupler/ensemble_manager.F90 fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c coupler/ensemble_manager.F90

fft.o: fft/fft.F90 include/file_version.h platform.o fms.o fft99.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude fft/fft.F90

fft99.o: fft/fft99.F90 constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c fft/fft99.F90

field_manager.o: field_manager/field_manager.F90 include/file_version.h field_manager/parse.inc mpp.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude -Ifield_manager field_manager/field_manager.F90

fm_util.o: field_manager/fm_util.F90 include/file_version.h field_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude field_manager/fm_util.F90

fms.o: fms/fms.F90 include/file_version.h mpp.o mpp_domains.o mpp_io.o fms_io.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude fms/fms.F90

fms2_io.o: fms2_io/fms2_io.F90 fms_io_utils.o netcdf_io.o fms_netcdf_domain_io.o fms_netcdf_unstructured_domain_io.o blackboxio.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c fms2_io/fms2_io.F90

fms_io.o: fms/fms_io.F90 include/fms_platform.h include/file_version.h fms/read_data_2d.inc fms/read_data_3d.inc fms/read_data_4d.inc fms/write_data.inc fms/fms_io_unstructured_register_restart_axis.inc fms/fms_io_unstructured_setup_one_field.inc fms/fms_io_unstructured_register_restart_field.inc fms/fms_io_unstructured_save_restart.inc fms/fms_io_unstructured_read.inc fms/fms_io_unstructured_get_file_name.inc fms/fms_io_unstructured_get_file_unit.inc fms/fms_io_unstructured_file_unit.inc fms/fms_io_unstructured_get_field_size.inc fms/fms_io_unstructured_field_exist.inc mpp_io.o mpp_domains.o mpp.o platform.o mpp_parameter.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude -Ifms fms/fms_io.F90

fms_io_utils.o: fms2_io/fms_io_utils.F90 fms2_io/include/array_utils.inc fms2_io/include/array_utils_char.inc fms2_io/include/get_data_type_string.inc fms2_io/include/get_checksum.inc mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Ifms2_io/include fms2_io/fms_io_utils.F90

fms_netcdf_domain_io.o: fms2_io/fms_netcdf_domain_io.F90 fms2_io/include/register_domain_restart_variable.inc fms2_io/include/domain_read.inc fms2_io/include/domain_write.inc mpp.o mpp_domains.o fms_io_utils.o netcdf_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Ifms2_io/include fms2_io/fms_netcdf_domain_io.F90

fms_netcdf_unstructured_domain_io.o: fms2_io/fms_netcdf_unstructured_domain_io.F90 fms2_io/include/register_unstructured_domain_restart_variable.inc fms2_io/include/unstructured_domain_read.inc fms2_io/include/unstructured_domain_write.inc mpp_domains.o fms_io_utils.o netcdf_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Ifms2_io/include fms2_io/fms_netcdf_unstructured_domain_io.F90

gaussian_topog.o: topography/gaussian_topog.F90 include/file_version.h fms.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude topography/gaussian_topog.F90

get_cal_time.o: time_manager/get_cal_time.F90 include/file_version.h fms.o time_manager.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude time_manager/get_cal_time.F90

gradient.o: mosaic/gradient.F90 include/file_version.h mpp.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude mosaic/gradient.F90

gradient_c2l.o: mosaic/gradient_c2l.c mosaic/constant.h mosaic/mosaic_util.h mosaic/gradient_c2l.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -Imosaic mosaic/gradient_c2l.c

grid.o: mosaic/grid.F90 include/file_version.h mpp.o constants.o fms.o fms_io.o mosaic.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude mosaic/grid.F90

horiz_interp.o: horiz_interp/horiz_interp.F90 include/file_version.h fms.o mpp.o constants.o horiz_interp_type.o horiz_interp_conserve.o horiz_interp_bilinear.o horiz_interp_bicubic.o horiz_interp_spherical.o mpp_io.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude horiz_interp/horiz_interp.F90

horiz_interp_bicubic.o: horiz_interp/horiz_interp_bicubic.F90 include/file_version.h mpp.o fms.o horiz_interp_type.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude horiz_interp/horiz_interp_bicubic.F90

horiz_interp_bilinear.o: horiz_interp/horiz_interp_bilinear.F90 include/file_version.h mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude horiz_interp/horiz_interp_bilinear.F90

horiz_interp_conserve.o: horiz_interp/horiz_interp_conserve.F90 include/file_version.h mpp.o fms.o fms_io.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude horiz_interp/horiz_interp_conserve.F90

horiz_interp_spherical.o: horiz_interp/horiz_interp_spherical.F90 include/file_version.h mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude horiz_interp/horiz_interp_spherical.F90

horiz_interp_type.o: horiz_interp/horiz_interp_type.F90 mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c horiz_interp/horiz_interp_type.F90

interp.o: mosaic/interp.c mosaic/mosaic_util.h mosaic/interp.h mosaic/create_xgrid.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -Imosaic mosaic/interp.c

interpolator.o: interpolator/interpolator.F90 include/fms_platform.h include/file_version.h mpp.o mpp_io.o mpp_domains.o diag_manager.o fms.o horiz_interp.o time_manager.o time_interp.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude interpolator/interpolator.F90

legacy.o: fms2_io/legacy.F90 fms2_io.o fms_io_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c fms2_io/legacy.F90

memuse.o: memutils/memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c memutils/memuse.c

memutils.o: memutils/memutils.F90 mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c memutils/memutils.F90

monin_obukhov.o: monin_obukhov/monin_obukhov.F90 constants.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c monin_obukhov/monin_obukhov.F90

monin_obukhov_kernel.o: monin_obukhov/monin_obukhov_kernel.F90 include/fms_platform.h monin_obukhov/monin_obukhov_interfaces.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude -Imonin_obukhov monin_obukhov/monin_obukhov_kernel.F90

mosaic.o: mosaic/mosaic.F90 include/file_version.h fms.o mpp.o mpp_io.o fms_io.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude mosaic/mosaic.F90

mosaic_util.o: mosaic/mosaic_util.c mosaic/mosaic_util.h mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -Imosaic mosaic/mosaic_util.c

mpp.o: mpp/mpp.F90 include/fms_platform.h include/file_version.h mpp/include/system_clock.h mpp/include/mpp_util.inc mpp/include/mpp_util_sma.inc mpp/include/mpp_util_mpi.inc mpp/include/mpp_util_nocomm.inc mpp/include/mpp_error_a_a.h mpp/include/mpp_error_a_s.h mpp/include/mpp_error_s_a.h mpp/include/mpp_error_s_s.h mpp/include/mpp_comm.inc mpp/include/mpp_comm_sma.inc mpp/include/mpp_transmit_sma.h mpp/include/mpp_transmit.inc mpp/include/mpp_reduce_sma.h mpp/include/mpp_sum_sma.h mpp/include/mpp_sum.inc mpp/include/mpp_alltoall_sma.h mpp/include/mpp_type_sma.h mpp/include/mpp_comm_mpi.inc mpp/include/mpp_transmit_mpi.h mpp/include/mpp_reduce_mpi.h mpp/include/mpp_sum_mpi.h mpp/include/mpp_sum_mpi_ad.h mpp/include/mpp_sum_ad.inc mpp/include/mpp_alltoall_mpi.h mpp/include/mpp_type_mpi.h mpp/include/mpp_comm_nocomm.inc mpp/include/mpp_transmit_nocomm.h mpp/include/mpp_reduce_nocomm.h mpp/include/mpp_sum_nocomm.h mpp/include/mpp_alltoall_nocomm.h mpp/include/mpp_type_nocomm.h mpp/include/mpp_chksum_int.h mpp/include/mpp_chksum_scalar.h mpp/include/mpp_chksum.h mpp/include/mpp_gather.h mpp/include/mpp_scatter.h mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude -Impp/include mpp/mpp.F90

mpp_data.o: mpp/mpp_data.F90 include/fms_platform.h include/file_version.h mpp/include/mpp_data_sma.inc mpp/include/mpp_data_mpi.inc mpp/include/mpp_data_nocomm.inc mpp_parameter.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude -Impp/include mpp/mpp_data.F90

mpp_domains.o: mpp/mpp_domains.F90 include/fms_platform.h include/file_version.h mpp/include/mpp_define_nest_domains.inc mpp/include/mpp_domains_util.inc mpp/include/mpp_domains_comm.inc mpp/include/mpp_domains_define.inc mpp/include/mpp_domains_misc.inc mpp/include/mpp_update_domains2D.h mpp/include/mpp_update_domains2D_nonblock.h mpp/include/mpp_do_update_nonblock.h mpp/include/mpp_do_updateV_nonblock.h mpp/include/mpp_do_update.h mpp/include/mpp_do_updateV.h mpp/include/mpp_do_check.h mpp/include/mpp_do_checkV.h mpp/include/mpp_update_nest_domains.h mpp/include/mpp_do_update_nest.h mpp/include/mpp_update_domains2D_ad.h mpp/include/mpp_do_update_ad.h mpp/include/mpp_do_updateV_ad.h mpp/include/mpp_do_redistribute.h mpp/include/mpp_get_boundary.h mpp/include/mpp_get_boundary_ad.h mpp/include/mpp_do_get_boundary.h mpp/include/mpp_do_get_boundary_ad.h mpp/include/mpp_group_update.h mpp/include/group_update_pack.inc mpp/include/group_update_unpack.inc mpp/include/mpp_domains_reduce.inc mpp/include/mpp_global_reduce.h mpp/include/mpp_global_sum.h mpp/include/mpp_global_sum_tl.h mpp/include/mpp_global_sum_ad.h mpp/include/mpp_global_field.h mpp/include/mpp_global_field_ad.h mpp/include/mpp_do_global_field.h mpp/include/mpp_do_global_field_ad.h mpp/include/mpp_unstruct_domain.inc mpp/include/mpp_unstruct_pass_data.h mpp/include/mpp_global_field_ug.h mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_pset.o mpp_efp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude -Impp/include mpp/mpp_domains.F90

mpp_efp.o: mpp/mpp_efp.F90 include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude mpp/mpp_efp.F90

mpp_io.o: mpp/mpp_io.F90 include/fms_platform.h include/file_version.h mpp/include/mpp_io_util.inc mpp/include/mpp_io_misc.inc mpp/include/mpp_io_connect.inc mpp/include/mpp_io_read.inc mpp/include/mpp_read_2Ddecomp.h mpp/include/mpp_read_compressed.h mpp/include/mpp_read_distributed_ascii.inc mpp/include/mpp_read_distributed_ascii.h mpp/include/mpp_io_write.inc mpp/include/mpp_write_2Ddecomp.h mpp/include/mpp_write_compressed.h mpp/include/mpp_write_unlimited_axis.h mpp/include/mpp_write.h mpp/include/mpp_io_unstructured_write.inc mpp/include/mpp_io_unstructured_read.inc mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude -Impp/include mpp/mpp_io.F90

mpp_memutils.o: mpp/mpp_memutils.F90 mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c mpp/mpp_memutils.F90

mpp_parameter.o: mpp/mpp_parameter.F90 include/fms_platform.h include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude mpp/mpp_parameter.F90

mpp_pset.o: mpp/mpp_pset.F90 include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude mpp/mpp_pset.F90

mpp_utilities.o: mpp/mpp_utilities.F90 include/file_version.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude mpp/mpp_utilities.F90

netcdf_io.o: fms2_io/netcdf_io.F90 fms2_io/include/netcdf_add_restart_variable.inc fms2_io/include/netcdf_read_data.inc fms2_io/include/netcdf_write_data.inc fms2_io/include/register_global_attribute.inc fms2_io/include/register_variable_attribute.inc fms2_io/include/get_global_attribute.inc fms2_io/include/get_variable_attribute.inc fms2_io/include/compressed_write.inc fms2_io/include/compressed_read.inc mpp.o fms_io_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Ifms2_io/include fms2_io/netcdf_io.F90

nsclock.o: mpp/nsclock.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c mpp/nsclock.c

oda_core.o: oda_tools/oda_core.F90 fms.o mpp.o mpp_domains.o time_manager.o get_cal_time.o axis_utils.o constants.o oda_types.o write_ocean_data.o mpp_io.o field_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c oda_tools/oda_core.F90

oda_core_ecda.o: oda_tools/oda_core_ecda.F90 fms.o mpp.o mpp_io.o mpp_domains.o mpp_memutils.o time_manager.o get_cal_time.o axis_utils.o horiz_interp_type.o horiz_interp_bilinear.o constants.o oda_types.o xbt_drop_rate_adjust.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c oda_tools/oda_core_ecda.F90

oda_types.o: oda_tools/oda_types.F90 time_manager.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c oda_tools/oda_types.F90

platform.o: platform/platform.F90 include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude platform/platform.F90

quicksort.o: drifters/quicksort.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c drifters/quicksort.F90

random_numbers.o: random_numbers/random_numbers.F90 MersenneTwister.o time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c random_numbers/random_numbers.F90

read_mosaic.o: mosaic/read_mosaic.c mosaic/read_mosaic.h mosaic/constant.h mosaic/mosaic_util.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -Imosaic mosaic/read_mosaic.c

sat_vapor_pres.o: sat_vapor_pres/sat_vapor_pres.F90 include/file_version.h constants.o fms.o mpp_io.o mpp.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude sat_vapor_pres/sat_vapor_pres.F90

sat_vapor_pres_k.o: sat_vapor_pres/sat_vapor_pres_k.F90 include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude sat_vapor_pres/sat_vapor_pres_k.F90

station_data.o: station_data/station_data.F90 include/file_version.h axis_utils.o mpp_io.o fms.o mpp.o mpp_domains.o diag_axis.o diag_output.o diag_manager.o diag_util.o time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude station_data/station_data.F90

stock_constants.o: exchange/stock_constants.F90 include/file_version.h mpp.o fms.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude exchange/stock_constants.F90

threadloc.o: mpp/threadloc.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c mpp/threadloc.c

time_interp.o: time_interp/time_interp.F90 include/file_version.h time_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude time_interp/time_interp.F90

time_interp_external.o: time_interp/time_interp_external.F90 include/fms_platform.h include/file_version.h fms.o mpp.o mpp_io.o time_manager.o get_cal_time.o mpp_domains.o time_interp.o axis_utils.o platform.o horiz_interp.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude time_interp/time_interp_external.F90

time_manager.o: time_manager/time_manager.F90 include/fms_platform.h include/file_version.h constants.o fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude time_manager/time_manager.F90

topography.o: topography/topography.F90 include/file_version.h gaussian_topog.o horiz_interp.o fms.o fms_io.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude topography/topography.F90

tracer_manager.o: tracer_manager/tracer_manager.F90 include/file_version.h mpp.o mpp_io.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude tracer_manager/tracer_manager.F90

tridiagonal.o: tridiagonal/tridiagonal.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c tridiagonal/tridiagonal.F90

write_ocean_data.o: oda_tools/write_ocean_data.F90 mpp_io.o mpp.o oda_types.o time_manager.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c oda_tools/write_ocean_data.F90

xbt_drop_rate_adjust.o: oda_tools/xbt_drop_rate_adjust.f90 oda_types.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c oda_tools/xbt_drop_rate_adjust.f90

xgrid.o: exchange/xgrid.F90 include/fms_platform.h include/file_version.h fms.o fms_io.o mpp.o mpp_domains.o mpp_io.o constants.o mosaic.o stock_constants.o gradient.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -Iinclude exchange/xgrid.F90

clean: neat
	-rm -f .libfms.a.cppdefs $(OBJ) libfms.a *.mod *__genmod.f90

neat:
	-rm -f $(TMPFILES)
