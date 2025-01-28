#!/bin/bash 
#
OPT=/opt/homebrew/opt

#musthave to deactivate conda
conda deactivate

/Users/petrumilev/Documents/phd/fxc_implementation/yambo/configure \
  CC=gcc-14 \
  CPP="gcc-14 -E -P" \
  MPICC=mpicc \
  FC=gfortran-14 \
  F77=gfortran-14 \
  FPP="gfortran-14 -E -P" \
  MPIFC=mpifort \
  --enable-keep-src \
  --enable-mpi \
  --enable-open-mp \
  --enable-time-profile \
  --enable-memory-profile \
  --enable-msgs-comps \
  --with-blas-libs="-L$OPT/openblas/lib -lopenblas" \
  --with-lapack-libs="-L$OPT/openblas/lib -lopenblas" \
  --with-libxc-path=$OPT/libxc \
  --with-fft-path=$OPT/fftw \
  --with-hdf5-path=$OPT/hdf5-mpi \
  --enable-par-linalg \
  --with-scalapack-libs="-L$OPT/scalapack/lib -lscalapack" \
  --with-blacs-libs="-L$OPT/scalapack/lib -lscalapack" \
  --enable-slepc-linalg \
  --with-extlibs-path=$HOME/opt/ext-libs

