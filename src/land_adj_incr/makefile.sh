#!/bin/ksh -l
set -x
#-----------------------------------------------------
#-use standard module.
#-----------------------------------------------------

export FCMP=mpiifort #ftn    #ifort
export INCS="-I${netcdf_fortran_ROOT}/include"          # netcdf modules does not 
export FFLAGS="$INCS -O3 -fp-model precise -r8 -convert big_endian -traceback -g -diag-disable=10448"
export LIBSM="-L${netcdf_fortran_ROOT}/lib -lnetcdff"


make -f Makefile_soilinc clean
make -f Makefile_soilinc
make -f Makefile_soilinc install
make -f Makefile_soilinc clean
