#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with the FV3 dynamical core.
#* If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

SHELL = /bin/sh

include ../../NEMS/src/conf/configure.nems

LIBRARY  = libfms.a
FMS_INSTALL=$(realpath $(PWD)/..)/FMS_INSTALL

FFLAGS   += -I../include -I../mpp/include -I../fms

SRCS_f   =

SRCS_f90 = \
		   ../oda_tools/xbt_drop_rate_adjust.f90

SRCS_F   =

SRCS_F90 = \
		   ../amip_interp/amip_interp.F90                                     \
		   ../astronomy/astronomy.F90                                         \
		   ../axis_utils/axis_utils.F90                                       \
		   ../block_control/block_control.F90                                 \
		   ../column_diagnostics/column_diagnostics.F90                       \
		   ../coupler/atmos_ocean_fluxes.F90                                  \
		   ../coupler/coupler_types.F90                                       \
		   ../coupler/ensemble_manager.F90                                    \
		   ../data_override/data_override.F90                                 \
		   ../diag_manager/diag_axis.F90                                      \
		   ../diag_manager/diag_data.F90                                      \
		   ../diag_manager/diag_grid.F90                                      \
		   ../diag_manager/diag_manifest.F90                                  \
		   ../diag_manager/diag_manager.F90                                   \
		   ../diag_manager/diag_output.F90                                    \
		   ../diag_manager/diag_table.F90                                     \
		   ../diag_manager/diag_util.F90                                      \
		   ../drifters/cloud_interpolator.F90                                 \
		   ../drifters/drifters.F90                                           \
		   ../drifters/drifters_comm.F90                                      \
		   ../drifters/drifters_core.F90                                      \
		   ../drifters/drifters_input.F90                                     \
		   ../drifters/drifters_io.F90                                        \
		   ../drifters/quicksort.F90                                          \
		   ../exchange/stock_constants.F90                                    \
		   ../exchange/test_xgrid.F90                                         \
		   ../exchange/xgrid.F90                                              \
		   ../fft/fft.F90                                                     \
		   ../fft/fft99.F90                                                   \
		   ../field_manager/field_manager.F90                                 \
		   ../field_manager/fm_util.F90                                       \
		   ../fms/fms.F90                                                     \
		   ../fms/fms_io.F90                                                  \
		   ../fms/test_fms_io.F90                                             \
		   ../constants/constants.F90                                         \
		   ../horiz_interp/horiz_interp.F90                                   \
		   ../horiz_interp/horiz_interp_bicubic.F90                           \
		   ../horiz_interp/horiz_interp_bilinear.F90                          \
		   ../horiz_interp/horiz_interp_conserve.F90                          \
		   ../horiz_interp/horiz_interp_spherical.F90                         \
		   ../horiz_interp/horiz_interp_type.F90                              \
		   ../horiz_interp/test_horiz_interp.F90                              \
		   ../interpolator/interpolator.F90                                   \
		   ../memutils/memutils.F90                                           \
		   ../mosaic/gradient.F90                                             \
		   ../mosaic/grid.F90                                                 \
		   ../mosaic/mosaic.F90                                               \
		   ../mpp/mpp.F90                                                     \
		   ../mpp/mpp_data.F90                                                \
		   ../mpp/mpp_domains.F90                                             \
		   ../mpp/mpp_efp.F90                                                 \
		   ../mpp/mpp_io.F90                                                  \
		   ../mpp/mpp_memutils.F90                                            \
		   ../mpp/mpp_parameter.F90                                           \
		   ../mpp/mpp_pset.F90                                                \
		   ../mpp/mpp_utilities.F90                                           \
		   ../mpp/test_mpp.F90                                                \
		   ../mpp/test_mpp_domains.F90                                        \
		   ../mpp/test_mpp_io.F90                                             \
		   ../mpp/test_mpp_pset.F90                                           \
		   ../oda_tools/oda_core.F90                                          \
		   ../oda_tools/oda_core_ecda.F90                                     \
		   ../oda_tools/oda_types.F90                                         \
		   ../oda_tools/write_ocean_data.F90                                  \
		   ../platform/platform.F90                                           \
		   ../random_numbers/MersenneTwister.F90                              \
		   ../random_numbers/random_numbers.F90                               \
		   ../sat_vapor_pres/sat_vapor_pres.F90                               \
		   ../sat_vapor_pres/sat_vapor_pres_k.F90                             \
		   ../station_data/station_data.F90                                   \
		   ../time_interp/time_interp.F90                                     \
		   ../time_interp/time_interp_external.F90                            \
		   ../time_manager/get_cal_time.F90                                   \
		   ../time_manager/time_manager.F90                                   \
		   ../topography/gaussian_topog.F90                                   \
		   ../topography/topography.F90                                       \
		   ../tracer_manager/tracer_manager.F90                               \
		   ../tridiagonal/tridiagonal.F90

SRCS_c   = \
		   ../memutils/memuse.c                                               \
		   ../mosaic/create_xgrid.c                                           \
		   ../mosaic/gradient_c2l.c                                           \
		   ../mosaic/interp.c                                                 \
		   ../mosaic/mosaic_util.c                                            \
		   ../mosaic/read_mosaic.c                                            \
		   ../mpp/affinity.c                                                  \
		   ../mpp/nsclock.c                                                   \
		   ../mpp/threadloc.c

DEPEND_FILES = $(SRCS_f) $(SRCS_f90) $(SRCS_F) $(SRCS_F90)

OBJS_f   = $(SRCS_f:.f=.o)
OBJS_f90 = $(SRCS_f90:.f90=.o)
OBJS_F   = $(SRCS_F:.F=.o)
OBJS_F90 = $(SRCS_F90:.F90=.o)
OBJS_c   = $(SRCS_c:.c=.o)

OBJS = $(OBJS_f) $(OBJS_f90) $(OBJS_F) $(OBJS_F90) $(OBJS_c)

all: esmf_make_fragment
	cp -fp *.a *.mod ../include/*.h "$(FMS_INSTALL)"/.

esmf_make_fragment: depend $(LIBRARY)
	mkdir -p "$(FMS_INSTALL)"
	@echo "# ESMF self-describing build dependency makefile fragment"              > fms.mk
	@echo "# src location /lfs3/projects/hfv3gfs/Samuel.Trahan/fms-component/FV3" >> fms.mk
	@echo ""                                                                      >> fms.mk
	@echo "ESMF_DEP_FRONT     ="                                                  >> fms.mk
	@echo "ESMF_DEP_INCPATH   = $(FMS_INSTALL)"                                   >> fms.mk
	@echo "ESMF_DEP_CMPL_OBJS ="                                                  >> fms.mk
	@echo "ESMF_DEP_LINK_OBJS = $(FMS_INSTALL)/libfms.a"                          >> fms.mk
	@echo "ESMF_DEP_SHRD_PATH ="                                                  >> fms.mk
	@echo "ESMF_DEP_SHRD_LIBS ="                                                  >> fms.mk
	mv fms.mk "$(FMS_INSTALL)/."

$(LIBRARY): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

.PHONY: clean esmf_make_fragment all
clean:
	@echo "Cleaning fms ... "
	@echo
	cd .. ; \
	ls -1 */*.a */*.o */*.mod */depend \
	  | grep -vE 'INSTALL|\.git' | xargs rm -f || true ; \
	rm -rf FMS_INSTALL || true

MKDEPENDS = ./mkDepends.pl
include ./make.rules

# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
    -include depend
endif

