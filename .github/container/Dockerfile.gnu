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
 && dnf install -y autoconf make automake m4 libtool pkg-config zip diffutils git \
 && rm -rf /var/cache/dnf && dnf clean all

RUN ranlib -U /opt/software/linux-rocky9-skylake/gcc-*/gcc-*/lib/gcc/x86_64-pc-linux-gnu/*/libgcc.a
ENV MPICH_FC=gfortran
ENV MPICH_CC=gcc
LABEL "maintainer"="Ryan Mulhall <Ryan.Mulhall@noaa.gov>"
LABEL "copyright"="2024 GFDL"
LABEL "license"="LGPL v3+"
LABEL "gov.noaa.gfdl.version"="1.0.0"
LABEL "vendor"="Geophysical Fluid Dynamics Laboratory"
LABEL "gov.noaa.gfdl.release-date"="2024-01-24"
ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash" ]

