# mkdir build && cd build && sh ../configure_cmake.sh
export CC=gcc
export CXX=g++
export MPI_C=mpicc
export MPI_CXX=mpicxx
export FFTW3_ROOT=$HOME/local/gnu
export HDF5_ROOT=$HOME/local/gnu
export GSL_ROOT=$HOME/local/gnu

cmake ..