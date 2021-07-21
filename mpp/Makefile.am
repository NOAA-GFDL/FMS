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
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# This is an automake file for the mpp directory of the MPP
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = \
  -I$(top_srcdir)/include \
  -I$(top_srcdir)/mpp/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

noinst_LTLIBRARIES = libmpp.la

libmpp_la_SOURCES = \
  mpp_parameter.F90 \
  mpp.F90 \
  mpp_data.F90 \
  mpp_utilities.F90 \
  mpp_domains.F90 \
  mpp_io.F90 \
  mpp_efp.F90 \
  mpp_memutils.F90 \
  mpp_memuse.c \
  mpp_pset.F90 \
  include/group_update_pack.inc \
  include/group_update_unpack.inc \
  include/mpp_alltoall_mpi.h \
  include/mpp_alltoall_nocomm.h \
  include/mpp_chksum.h \
  include/mpp_chksum_int.h \
  include/mpp_chksum_scalar.h \
  include/mpp_comm.inc \
  include/mpp_comm_mpi.inc \
  include/mpp_comm_nocomm.inc \
  include/mpp_data_mpi.inc \
  include/mpp_data_nocomm.inc \
  include/mpp_define_nest_domains.inc \
  include/mpp_do_check.h \
  include/mpp_do_checkV.h \
  include/mpp_do_get_boundary.h \
  include/mpp_do_get_boundary_ad.h \
  include/mpp_do_global_field.h \
  include/mpp_do_global_field_ad.h \
  include/mpp_do_redistribute.h \
  include/mpp_do_update.h \
  include/mpp_do_updateV.h \
  include/mpp_do_updateV_ad.h \
  include/mpp_do_updateV_nonblock.h \
  include/mpp_do_update_ad.h \
  include/mpp_do_update_nest.h \
  include/mpp_do_update_nonblock.h \
  include/mpp_domains_comm.inc \
  include/mpp_domains_define.inc \
  include/mpp_domains_misc.inc \
  include/mpp_domains_reduce.inc \
  include/mpp_domains_util.inc \
  include/mpp_error_a_a.h \
  include/mpp_error_a_s.h \
  include/mpp_error_s_a.h \
  include/mpp_error_s_s.h \
  include/mpp_gather.h \
  include/mpp_get_boundary.h \
  include/mpp_get_boundary_ad.h \
  include/mpp_global_field.h \
  include/mpp_global_field_ad.h \
  include/mpp_global_field_ug.h \
  include/mpp_global_reduce.h \
  include/mpp_global_sum.h \
  include/mpp_global_sum_ad.h \
  include/mpp_global_sum_tl.h \
  include/mpp_group_update.h \
  include/mpp_io_connect.inc \
  include/mpp_io_misc.inc \
  include/mpp_io_read.inc \
  include/mpp_io_unstructured_read.inc \
  include/mpp_io_unstructured_write.inc \
  include/mpp_io_util.inc \
  include/mpp_io_write.inc \
  include/mpp_read_2Ddecomp.h \
  include/mpp_read_compressed.h \
  include/mpp_read_distributed_ascii.h \
  include/mpp_read_distributed_ascii.inc \
  include/mpp_reduce_mpi.h \
  include/mpp_reduce_nocomm.h \
  include/mpp_scatter.h \
  include/mpp_sum.inc \
  include/mpp_sum_ad.inc \
  include/mpp_sum_mpi.h \
  include/mpp_sum_mpi_ad.h \
  include/mpp_sum_nocomm.h \
  include/mpp_sum_nocomm_ad.h \
  include/mpp_transmit.inc \
  include/mpp_transmit_mpi.h \
  include/mpp_transmit_nocomm.h \
  include/mpp_type_mpi.h \
  include/mpp_type_nocomm.h \
  include/mpp_unstruct_domain.inc \
  include/mpp_unstruct_pass_data.h \
  include/mpp_update_domains2D.h \
  include/mpp_update_domains2D_ad.h \
  include/mpp_update_domains2D_nonblock.h \
  include/mpp_update_nest_domains.h \
  include/mpp_util.inc \
  include/mpp_util_mpi.inc \
  include/mpp_util_nocomm.inc \
  include/mpp_write.h \
  include/mpp_write_2Ddecomp.h \
  include/mpp_write_compressed.h \
  include/mpp_write_unlimited_axis.h \
  include/system_clock.h

mpp_mod.$(FC_MODEXT): \
  mpp_parameter_mod.$(FC_MODEXT) \
  mpp_data_mod.$(FC_MODEXT) \
  include/mpp_util.inc \
  include/mpp_error_a_a.h \
  include/mpp_error_a_s.h \
  include/mpp_error_s_a.h \
  include/mpp_error_s_s.h \
  include/mpp_util_mpi.inc \
  include/mpp_util_nocomm.inc \
  include/mpp_comm.inc \
  include/mpp_chksum.h \
  include/mpp_chksum_int.h \
  include/mpp_chksum_scalar.h \
  include/mpp_comm_mpi.inc \
  include/mpp_alltoall_mpi.h \
  include/mpp_reduce_mpi.h \
  include/mpp_sum_mpi.h \
  include/mpp_sum.inc \
  include/mpp_sum_mpi_ad.h \
  include/mpp_sum_ad.inc \
  include/mpp_transmit_mpi.h \
  include/mpp_transmit.inc \
  include/mpp_type_mpi.h \
  include/mpp_comm_nocomm.inc \
  include/mpp_alltoall_nocomm.h \
  include/mpp_reduce_nocomm.h \
  include/mpp_sum_nocomm.h \
  include/mpp_sum_nocomm_ad.h \
  include/mpp_transmit_nocomm.h \
  include/mpp_type_nocomm.h \
  include/mpp_gather.h \
  include/mpp_scatter.h \
  include/system_clock.h
mpp_data_mod.$(FC_MODEXT): \
  mpp_parameter_mod.$(FC_MODEXT) \
  include/mpp_data_mpi.inc \
  include/mpp_data_nocomm.inc
mpp_utilities_mod.$(FC_MODEXT): \
  mpp_mod.$(FC_MODEXT) \
  mpp_efp_mod.$(FC_MODEXT)
mpp_domains_mod.$(FC_MODEXT): \
  mpp_data_mod.$(FC_MODEXT) \
  mpp_parameter_mod.$(FC_MODEXT) \
  mpp_mod.$(FC_MODEXT) \
  mpp_memutils_mod.$(FC_MODEXT) \
  mpp_efp_mod.$(FC_MODEXT) \
  include/mpp_define_nest_domains.inc \
  include/mpp_domains_util.inc \
  include/mpp_domains_comm.inc \
  include/mpp_domains_define.inc \
  include/mpp_domains_misc.inc \
  include/mpp_do_check.h \
  include/mpp_do_checkV.h \
  include/mpp_do_get_boundary.h \
  include/mpp_do_get_boundary_ad.h \
  include/mpp_do_redistribute.h \
  include/mpp_do_update.h \
  include/mpp_do_updateV.h \
  include/mpp_do_updateV_ad.h \
  include/mpp_do_updateV_nonblock.h \
  include/mpp_do_update_ad.h \
  include/mpp_do_update_nest.h \
  include/mpp_do_update_nonblock.h \
  include/mpp_get_boundary.h \
  include/mpp_get_boundary_ad.h \
  include/mpp_group_update.h \
  include/group_update_pack.inc \
  include/group_update_unpack.inc \
  include/mpp_update_domains2D.h \
  include/mpp_update_domains2D_ad.h \
  include/mpp_update_domains2D_nonblock.h \
  include/mpp_update_nest_domains.h \
  include/mpp_domains_reduce.inc \
  include/mpp_do_global_field.h \
  include/mpp_do_global_field_ad.h \
  include/mpp_global_field.h \
  include/mpp_global_field_ad.h \
  include/mpp_global_reduce.h \
  include/mpp_global_sum.h \
  include/mpp_global_sum_ad.h \
  include/mpp_global_sum_tl.h \
  include/mpp_unstruct_domain.inc \
  include/mpp_global_field_ug.h \
  include/mpp_unstruct_pass_data.h
mpp_io_mod.$(FC_MODEXT): \
  mpp_parameter_mod.$(FC_MODEXT) \
  mpp_mod.$(FC_MODEXT) \
  mpp_domains_mod.$(FC_MODEXT) \
  include/mpp_io_util.inc \
  include/mpp_io_misc.inc \
  include/mpp_io_connect.inc \
  include/mpp_io_read.inc \
  include/mpp_read_2Ddecomp.h \
  include/mpp_read_compressed.h \
  include/mpp_read_distributed_ascii.inc \
  include/mpp_read_distributed_ascii.h \
  include/mpp_io_write.inc \
  include/mpp_write.h \
  include/mpp_write_2Ddecomp.h \
  include/mpp_write_compressed.h \
  include/mpp_write_unlimited_axis.h \
  include/mpp_io_unstructured_write.inc \
  include/mpp_io_unstructured_read.inc
mpp_efp_mod.$(FC_MODEXT): mpp_parameter_mod.$(FC_MODEXT) mpp_mod.$(FC_MODEXT)
mpp_memutils_mod.$(FC_MODEXT): mpp_mod.$(FC_MODEXT)
mpp_pset_mod.$(FC_MODEXT): mpp_mod.$(FC_MODEXT)

# Mod files are built and then installed as headers.
BUILT_SOURCES = \
  mpp_parameter_mod.$(FC_MODEXT) \
  mpp_data_mod.$(FC_MODEXT) \
  mpp_mod.$(FC_MODEXT) \
  mpp_utilities_mod.$(FC_MODEXT) \
  mpp_memutils_mod.$(FC_MODEXT) \
  mpp_efp_mod.$(FC_MODEXT) \
  mpp_domains_mod.$(FC_MODEXT) \
  mpp_io_mod.$(FC_MODEXT) \
  mpp_pset_mod.$(FC_MODEXT)
nodist_include_HEADERS = $(BUILT_SOURCES)

include $(top_srcdir)/mkmods.mk
