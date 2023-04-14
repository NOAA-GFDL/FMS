# base from spack org
# needs a shell cmd to run interactively
# TODO could be a little slimmer (remove gcc base+deps if possible)
FROM spack/centos-stream:v0.19.1

# arguments to specify versions to build (--build-arg name=val)
ARG gcc_version=12.2.0
ARG netcdfc_version=4.9.0
ARG netcdff_version=4.6.0
ARG libyaml_version=0.2.5
ARG mpich_version=4.0.2

# copy over spack environment template yaml
COPY spack.env .

# install comp version needed via spack
# set up spack template yaml file with desired versions
# install all dependencies 
RUN spack install gcc@${gcc_version}                          && \
    source /opt/spack/share/spack/setup-env.sh                && \
    spack load gcc@${gcc_version}                             && \
    spack compiler find                                       && \
    sed "s/COMPILER/gcc@$gcc_version/" spack.env > spack.yaml && \
    sed -i "s/NETCDF_C_VERSION/$netcdfc_version/" spack.yaml  && \
    sed -i "s/NETCDF_F_VERSION/$netcdff_version/" spack.yaml  && \
    sed -i "s/LIBYAML_VERSION/$libyaml_version/" spack.yaml   && \
    sed -i "s/MPI_LIB/mpich@$mpich_version/" spack.yaml       && \
    spack env activate -d .                                   && \
    spack -e . concretize -f > deps.log                       && \
    spack env activate . && spack install

# flag/lib vars need nc/nf-config to get path with spack hash
# can't run commands in ENV so sets them in a entrypoint script
# shell expansions are escaped in order to run in the entrypoint (except for gcc_version)
# TODO is there a spack command to get install directory? would be better than hardcoding libyaml paths
RUN echo "#!/usr/bin/bash"                                                                   > /root/env.sh && \
    echo "source /opt/spack/share/spack/setup-env.sh"                                        >>/root/env.sh && \
    echo "spack load gcc@${gcc_version} mpich netcdf-c netcdf-fortran libyaml"               >>/root/env.sh && \
    echo "export LIBYAML_PATH=\$(find /opt/spack/opt/spack/linux-centos8-haswell/gcc-${gcc_version} -maxdepth 1 -type d -name libyaml\*)" >> /root/env.sh && \
    echo "export FC=mpifort CC=mpicc"                                                        >>/root/env.sh && \
    echo "export FCFLAGS=\"\$(nf-config --fflags) -g\""                                      >>/root/env.sh && \
    echo "export CFLAGS=\"\$(nc-config --cflags) -I\${LIBYAML_PATH}/include -g\""                                       >>/root/env.sh && \
    echo "export LDFLAGS=\"\$(nc-config --libs) \$(nf-config --flibs) -L\${LIBYAML_PATH}/lib \""                     >>/root/env.sh && \
    echo "export LD_LIBRARY_PATH=\"\$(nc-config --prefix)/lib:\$(nf-config --prefix)/lib:\${LD_LIBRARY_PATH}\"" >>/root/env.sh && \
    echo "export LD_LIBRARY_PATH=\"\${LIBYAML_PATH}/lib:\${LD_LIBRARY_PATH}\"" >>/root/env.sh && \
    echo "ulimit -s unlimited"                                                               >>/root/env.sh && \
    echo "exec \"\$@\"" >>/root/env.sh && \
    chmod +x /root/env.sh

ENTRYPOINT ["/root/env.sh"]
