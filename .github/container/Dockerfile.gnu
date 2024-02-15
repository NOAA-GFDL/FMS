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
# Build stage with Spack pre-installed and ready to be used
FROM spack/rockylinux9:latest as builder


# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir /opt/spack-environment \
&&  (echo spack: \
&&   echo '  specs:' \
&&   echo '  - gcc@13.2.0' \
&&   echo '  - mpich@4.1.2' \
&&   echo '  - netcdf-c@4.9.2 ^hdf5@1.14.2' \
&&   echo '  - netcdf-fortran@4.6.1' \
&&   echo '  - libyaml@0.2.5' \
&&   echo '  concretizer:' \
&&   echo '    unify: true' \
&&   echo '    reuse: true' \
&&   echo '  packages:' \
&&   echo '    all:' \
&&   echo '      compiler: [gcc@13.2.0]' \
&&   echo '  # template file for any extra steps' \
&&   echo '  config:' \
&&   echo '    template_dirs:' \
&&   echo '    - ./' \
&&   echo '    install_tree: /opt/software' \
&&   echo '  view: /opt/views/view') > /opt/spack-environment/spack.yaml

# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && spack env activate . && spack install --fail-fast && spack gc -y

# Strip all the binaries
RUN find -L /opt/views/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-archive\|x-sharedlib' | \
    awk -F: '{print $1}' | xargs strip

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . > activate.sh



# Bare OS image to run the installed executables
FROM docker.io/rockylinux:9

COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software

# paths.view is a symlink, so copy the parent to avoid dereferencing and duplicating it
COPY --from=builder /opt/views /opt/views

RUN { \
      echo '#!/bin/sh' \
      && echo '.' /opt/spack-environment/activate.sh \
      && echo 'exec "$@"'; \
    } > /entrypoint.sh \
&& chmod a+x /entrypoint.sh \
&& ln -s /opt/views/view /opt/view


RUN dnf update -y && dnf install -y epel-release && dnf update -y \
 && dnf install -y autoconf make automake m4 libtool pkg-config zip diffutils git libgomp cmake \
 && rm -rf /var/cache/dnf && dnf clean all

# Copy input files for certain unit tests
COPY ./fms_test_input /home/unit_tests_input
# Fix libgcc link
RUN ranlib -U /opt/software/linux-rocky9-skylake/gcc-*/gcc-*/lib/gcc/x86_64-pc-linux-gnu/*/libgcc.a
# Any env vars needed, most of these are also set by the entrypoint but it is not run for github actions
ENV MPICH_FC=/opt/views/view/bin/gfortran
ENV MPICH_CC=/opt/views/view/bin/gcc
ENV FC=/opt/views/view/bin/mpifort
ENV CC=/opt/views/view/bin/mpicc
ENV FCFLAGS="-I/opt/views/view/include"
ENV CFLAGS="-I/opt/views/view/include"
ENV LDFLAGS="-L/opt/views/view/lib"
ENV LD_LIBRARY_PATH="/opt/views/view/lib64:/opt/views/view/lib"
ENV PATH="/root/.local/bin:/root/bin:/opt/views/view/bin:/opt/spack/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
LABEL "maintainer"="Ryan Mulhall <Ryan.Mulhall@noaa.gov>"
LABEL "copyright"="2024 GFDL"
LABEL "license"="LGPL v3+"
LABEL "gov.noaa.gfdl.version"="1.0.0"
LABEL "vendor"="Geophysical Fluid Dynamics Laboratory"
LABEL "gov.noaa.gfdl.release-date"="2024-02-06"
ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash" ]

