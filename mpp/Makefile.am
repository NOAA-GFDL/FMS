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

# Include .fh and .mod files.
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
  include/mpp_alltoall_mpi.fh \
  include/mpp_alltoall_nocomm.fh \
  include/mpp_chksum.fh \
  include/mpp_chksum_int.fh \
  include/mpp_chksum_scalar.fh \
  include/mpp_comm.inc \
  include/mpp_comm_mpi.inc \
  include/mpp_comm_nocomm.inc \
  include/mpp_data_mpi.inc \
  include/mpp_data_nocomm.inc \
  include/mpp_define_nest_domains.inc \
  include/mpp_do_check.fh \
  include/mpp_do_checkV.fh \
  include/mpp_do_get_boundary.fh \
  include/mpp_do_get_boundary_ad.fh \
  include/mpp_do_global_field.fh \
  include/mpp_do_global_field_ad.fh \
  include/mpp_do_redistribute.fh \
  include/mpp_do_update.fh \
  include/mpp_do_updateV.fh \
  include/mpp_do_updateV_ad.fh \
  include/mpp_do_updateV_nonblock.fh \
  include/mpp_do_update_ad.fh \
  include/mpp_do_update_nest.fh \
  include/mpp_do_update_nonblock.fh \
  include/mpp_domains_comm.inc \
  include/mpp_domains_define.inc \
  include/mpp_domains_misc.inc \
  include/mpp_domains_reduce.inc \
  include/mpp_domains_util.inc \
  include/mpp_error_a_a.fh \
  include/mpp_error_a_s.fh \
  include/mpp_error_s_a.fh \
  include/mpp_error_s_s.fh \
  include/mpp_gather.fh \
  include/mpp_get_boundary.fh \
  include/mpp_get_boundary_ad.fh \
  include/mpp_global_field.fh \
  include/mpp_global_field_ad.fh \
  include/mpp_global_field_ug.fh \
  include/mpp_global_reduce.fh \
  include/mpp_global_sum.fh \
  include/mpp_global_sum_ad.fh \
  include/mpp_global_sum_tl.fh \
  include/mpp_group_update.fh \
  include/mpp_io_connect.inc \
  include/mpp_io_misc.inc \
  include/mpp_io_read.inc \
  include/mpp_io_unstructured_read.inc \
  include/mpp_io_unstructured_write.inc \
  include/mpp_io_util.inc \
  include/mpp_io_write.inc \
  include/mpp_read_2Ddecomp.fh \
  include/mpp_read_compressed.fh \
  include/mpp_read_distributed_ascii.fh \
  include/mpp_read_distributed_ascii.inc \
  include/mpp_reduce_mpi.fh \
  include/mpp_reduce_nocomm.fh \
  include/mpp_scatter.fh \
  include/mpp_sum.inc \
  include/mpp_sum_ad.inc \
  include/mpp_sum_mpi.fh \
  include/mpp_sum_mpi_ad.fh \
  include/mpp_sum_nocomm.fh \
  include/mpp_sum_nocomm_ad.fh \
  include/mpp_transmit.inc \
  include/mpp_transmit_mpi.fh \
  include/mpp_transmit_nocomm.fh \
  include/mpp_type_mpi.fh \
  include/mpp_type_nocomm.fh \
  include/mpp_unstruct_domain.inc \
  include/mpp_unstruct_pass_data.fh \
  include/mpp_update_domains2D.fh \
  include/mpp_update_domains2D_ad.fh \
  include/mpp_update_domains2D_nonblock.fh \
  include/mpp_update_nest_domains.fh \
  include/mpp_util.inc \
  include/mpp_util_mpi.inc \
  include/mpp_util_nocomm.inc \
  include/mpp_write.fh \
  include/mpp_write_2Ddecomp.fh \
  include/mpp_write_compressed.fh \
  include/mpp_write_unlimited_axis.fh \
  include/system_clock.fh

mpp_mod.$(FC_MODEXT): \
  mpp_parameter_mod.$(FC_MODEXT) \
  mpp_data_mod.$(FC_MODEXT) \
  include/mpp_util.inc \
  include/mpp_error_a_a.fh \
  include/mpp_error_a_s.fh \
  include/mpp_error_s_a.fh \
  include/mpp_error_s_s.fh \
  include/mpp_util_mpi.inc \
  include/mpp_util_nocomm.inc \
  include/mpp_comm.inc \
  include/mpp_chksum.fh \
  include/mpp_chksum_int.fh \
  include/mpp_chksum_scalar.fh \
  include/mpp_comm_mpi.inc \
  include/mpp_alltoall_mpi.fh \
  include/mpp_reduce_mpi.fh \
  include/mpp_sum_mpi.fh \
  include/mpp_sum.inc \
  include/mpp_sum_mpi_ad.fh \
  include/mpp_sum_ad.inc \
  include/mpp_transmit_mpi.fh \
  include/mpp_transmit.inc \
  include/mpp_type_mpi.fh \
  include/mpp_comm_nocomm.inc \
  include/mpp_alltoall_nocomm.fh \
  include/mpp_reduce_nocomm.fh \
  include/mpp_sum_nocomm.fh \
  include/mpp_sum_nocomm_ad.fh \
  include/mpp_transmit_nocomm.fh \
  include/mpp_type_nocomm.fh \
  include/mpp_gather.fh \
  include/mpp_scatter.fh \
  include/system_clock.fh
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
  include/mpp_do_check.fh \
  include/mpp_do_checkV.fh \
  include/mpp_do_get_boundary.fh \
  include/mpp_do_get_boundary_ad.fh \
  include/mpp_do_redistribute.fh \
  include/mpp_do_update.fh \
  include/mpp_do_updateV.fh \
  include/mpp_do_updateV_ad.fh \
  include/mpp_do_updateV_nonblock.fh \
  include/mpp_do_update_ad.fh \
  include/mpp_do_update_nest.fh \
  include/mpp_do_update_nonblock.fh \
  include/mpp_get_boundary.fh \
  include/mpp_get_boundary_ad.fh \
  include/mpp_group_update.fh \
  include/group_update_pack.inc \
  include/group_update_unpack.inc \
  include/mpp_update_domains2D.fh \
  include/mpp_update_domains2D_ad.fh \
  include/mpp_update_domains2D_nonblock.fh \
  include/mpp_update_nest_domains.fh \
  include/mpp_domains_reduce.inc \
  include/mpp_do_global_field.fh \
  include/mpp_do_global_field_ad.fh \
  include/mpp_global_field.fh \
  include/mpp_global_field_ad.fh \
  include/mpp_global_reduce.fh \
  include/mpp_global_sum.fh \
  include/mpp_global_sum_ad.fh \
  include/mpp_global_sum_tl.fh \
  include/mpp_unstruct_domain.inc \
  include/mpp_global_field_ug.fh \
  include/mpp_unstruct_pass_data.fh
mpp_io_mod.$(FC_MODEXT): \
  mpp_parameter_mod.$(FC_MODEXT) \
  mpp_mod.$(FC_MODEXT) \
  mpp_domains_mod.$(FC_MODEXT) \
  include/mpp_io_util.inc \
  include/mpp_io_misc.inc \
  include/mpp_io_connect.inc \
  include/mpp_io_read.inc \
  include/mpp_read_2Ddecomp.fh \
  include/mpp_read_compressed.fh \
  include/mpp_read_distributed_ascii.inc \
  include/mpp_read_distributed_ascii.fh \
  include/mpp_io_write.inc \
  include/mpp_write.fh \
  include/mpp_write_2Ddecomp.fh \
  include/mpp_write_compressed.fh \
  include/mpp_write_unlimited_axis.fh \
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
