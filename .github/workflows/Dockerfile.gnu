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
# FMS CI image recipefile for GNU
# Runs on centos stream (builder has same base from redhat registry)
#
# arguments to specify versions to build can be given to docker or changed here (--build-arg name=val)
FROM spack/rockylinux9:latest as builder

ARG gcc_version=12.3.0
ARG netcdfc_version=4.9.0
ARG netcdff_version=4.6.0
ARG libyaml_version=0.2.5
ARG mpich_version=4.0.2

COPY spack.env /opt/deps/spack.env

# perl's download kept timing out
RUN sed -i 's/connect_timeout: 10/connect_timeout: 600/' /opt/spack/etc/spack/defaults/config.yaml && \
    spack install gcc@${gcc_version}                          && \
    source /opt/spack/share/spack/setup-env.sh                && \
    spack load gcc@${gcc_version}                             && \
    spack compiler find                                       && \
    sed "s/COMPILER/gcc@$gcc_version/" /opt/deps/spack.env > spack.yaml && \
    sed -i "s/NETCDF_C_VERSION/$netcdfc_version/" spack.yaml  && \
    sed -i "s/NETCDF_F_VERSION/$netcdff_version/" spack.yaml  && \
    sed -i "s/LIBYAML_VERSION/$libyaml_version/" spack.yaml   && \
    sed -i "s/MPI_LIB/mpich@$mpich_version/" spack.yaml       && \
    spack env activate -d .                                   && \
    spack -e . concretize -f > /opt/deps/deps.log             && \
    spack install --fail-fast

# copy built software to base from first image
FROM rockylinux:9

COPY --from=builder /opt/view/ /opt/view/
COPY --from=builder /opt/deps/ /opt/deps/

# input files used with --enable-input-tests
# need to be on the dev boxes if building
COPY ./fms_test_input /home/unit_tests_input

RUN dnf install -y autoconf make automake m4 libtool pkg-config zip

ENV FC="mpifort"
ENV CC="mpicc"
ENV MPICH_FC="/opt/view/bin/gfortran"
ENV MPICH_CC="/opt/view/bin/gcc"
ENV FCFLAGS="-I/opt/view/include"
ENV CFLAGS="-I/opt/view/include"
ENV LDFLAGS="-L/opt/view/lib"
ENV LD_LIBRARY_PATH="/opt/view/lib:/opt/view/lib64:/usr/local/lib:/usr/local/lib64"
ENV PATH="/opt/view/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin"
