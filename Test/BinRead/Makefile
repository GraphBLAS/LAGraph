#-------------------------------------------------------------------------------
# Test/BinRead/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
# http://graphblas.org  See LAGraph/LICENSE for license.

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph/Test/BinRead/bc_test,

# relies on cmake to do the actual build.  Use the CMAKE_OPTIONS argument to
# this Makefile to pass options to cmake.

# Install LAGraph and GraphBLAS before trying to compile.

JOBS ?= 1

# build and run test
default: compile
	./build/mtx2bin west0067.mtx west0067.grb
	./build/binread west0067.grb
	./build/binread west0067.grb west0067_int32.grb int32
	./build/binread west0067_int32.grb

# build test
compile:
	( cd build ; cmake $(CMAKE_OPTIONS) .. ; $(MAKE) --jobs=$(JOBS) )

# just run cmake; do not compile
cmake:
	( cd build ; cmake $(CMAKE_OPTIONS) .. ; )

clean: distclean

purge: distclean

# remove all files not in the distribution
distclean:
	rm -rf build/* west0067.grb west0067_int32.grb

