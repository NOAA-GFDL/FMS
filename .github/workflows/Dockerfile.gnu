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
FROM spack/centos-stream:v0.19.1 as builder

ARG gcc_version=12.2.0
ARG netcdfc_version=4.9.0
ARG netcdff_version=4.6.0
ARG libyaml_version=0.2.5
ARG mpich_version=4.0.2

COPY spack.env /opt/deps/spack.env

RUN spack install gcc@${gcc_version}                          && \
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

RUN find -L /opt/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-arcive\|x-sharedlib' | \
    awk -F: '{print $1}' |  xargs strip -s

# copy built software to base from first image
FROM quay.io/centos/centos:stream

COPY --from=builder /opt/deps/ /opt/deps/
COPY --from=builder /opt/spack/opt/spack/linux-centos8-haswell/gcc-8.5.0/  /opt/spack/opt/spack/linux-centos8-haswell/gcc-8.5.0/

# input files used with --enable-input-tests
# need to be on the dev boxes if building
COPY /home/fms_test_input /home/fms_test_input

# extends glob patterns for exceptions
SHELL ["/bin/bash", "-O", "extglob", "-c"]

RUN ln -s /opt/deps/linux-centos8-haswell/gcc-12.2.0/*/lib/!(pkgconfig|cmake) /usr/local/lib && \
    ln -s /opt/deps/linux-centos8-haswell/gcc-12.2.0/*/bin/* /usr/local/bin && \
    ln -s /opt/deps/linux-centos8-haswell/gcc-12.2.0/*/include/* /usr/local/include && \
    dnf install -y autoconf automake make cmake binutils glibc-devel

ENV FC="mpifort"
ENV CC="mpicc"
ENV FCFLAGS="-I/usr/local/include"
ENV CFLAGS="-I/usr/local/include"
