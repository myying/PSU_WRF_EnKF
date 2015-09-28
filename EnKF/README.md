About the derivation and work flow of the algorithm, see comments in enkf.f

How to compile:

I. using PGI compiler
---------------------

Make sure you have pgi compiler, mpich and netcdf correctly configured in your environment.

We recommend you having the mpich platform and netcdf compiled by pgi as well.

The environment variables PATH should have the search paths of PGI, NETCDF, and MPICH binaries.
and LD_LIBRARY_PATH should have the search paths of PGI, NETCDF, and MPICH shared libraries.

The following is a sample .cshrc file configured for using pgi, mvapich2 and netcdf (on mc1.met.psu.edu):
    setenv PGI /usr/global/pgi
    setenv PATH $PGI/linux86-64/10.6/bin:$PATH

    setenv MPICH /usr/mpi/pgi/mvapich2-1.5.1
    setenv PATH $MPICH/bin:$PATH

    setenv NETCDF /usr/global/netcdf
    setenv PATH $NETCDF/bin:$PATH
    setenv LD_LIBRARY_PATH $NETCDF/lib:$LD_LIBRARY_PATH

To compile:
    cd src
    cp Makefile.pgi Makefile
    make clean
    make


II. using Intel compiler
------------------------

Similar with PGI compiler, the environment should be correctly configured for intel, mpich and netcdf.

The following is a sample .bashrc file configured for using intel, mvapich2 and netcdf (on stampede.tacc.utexas.edu):
    export NETCDF="/opt/apps/intel13/netcdf/3.6.3"
    export PATH="$NETCDF/bin:$PATH"
    export LD_LIBRARY_PATH="$NETCDF/lib:$LD_LIBRARY_PATH"

    export INTEL="/opt/apps/intel/13/composer_xe_2013.1.117"
    export PATH="$INTEL/bin/intel64:$PATH"
    export LD_LIBRARY_PATH="$INTEL/tbb/lib/intel64:$INTEL/mkl/lib/intel64:$INTEL/ipp/lib/intel64:$INTEL/mpirt/lib/intel64:/opt/intel/mic/myo/lib:/opt/intel/mic/coi/host-linux-release/lib:$INTEL/compiler/lib/intel64:$LD_LIBRARY_PATH"

    export MPICH="/opt/apps/intel13/mvapich2/1.9"
    export PATH="$MPICH/bin:$PATH"
    export LD_LIBRARY_PATH="$MPICH/lib:$MPICH/lib/shared:$LD_LIBRARY_PATH"

To compile:
    cd src
    cp Makefile.intel Makefile
    make clean
    make


Executables:
------------
    If build is successful, you can find the following executables:

    enkf.mpi  - the main program of EnKF

    Also some utility programs:
    enkf_gfs_hybrid_inflation.exe
    replace_perturbationmean_by_initial.exe
    replace_geo_by_initial.exe
    replace_xfmean_with_gfs.exe 
    soairborne_to_3dvar.exe
    so_to_3dvar.exe

    If you only need enkf.mpi, simply use "make clean; make enkf.mpi"
