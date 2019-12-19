# monofonIC

High order LPT/QPT tool for single resolution simulations

## Build Instructions
Clone code including submodules (currently only CLASS is used as a submodule):

    git clone --recurse-submodules https://ohahn@bitbucket.org/ohahn/monofonic.git


Create build directory, configure, and build:

    mkdir monofonic/build; cd monofonic/build
	
    ccmake ..
	
    make

this should create an executable in the build directory. 
There is an example parameter file 'example.conf' in the main directory

If you run into problems with CMake not being able to find your local FFTW3 or HDF5 installation, it is best to give the path directly as

    FFTW3_ROOT=<path> HDF5_ROOT=<path> ccmake ..

make sure to delete previous files generated by CMake before reconfiguring like this.