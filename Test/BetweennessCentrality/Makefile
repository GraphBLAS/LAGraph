#-------------------------------------------------------------------------------
# Test/BetweennessCentrality/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
# http://graphblas.org  See LAGraph/LICENSE for license.

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph/Test/BetweennessCentrality/bc_test,

# relies on cmake to do the actual build.  Use the CMAKE_OPTIONS argument to
# this Makefile to pass options to cmake.

# Install LAGraph and GraphBLAS before trying to compile.

JOBS ?= 1

# build and run test
default: compile
	./build/bc_test < sample.mtx 2>stderr.txt
	./build/bc_test < sample2.mtx 2>>stderr.txt
	./build/bc_test < west0067.mtx 2>>stderr.txt
	./build/gap_bc_test < LFAT5.mtx
	./build/gap_bc_test < west0067.mtx
	./build/gap_bc_test_allthreads < west0067.mtx
#	./build/bc_test < jagmesh7.mtx 2>>stderr.txt
#	./build/gap_bc_test ~/indochina-2004.grb ~/indochina-sources.mtx

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
	rm -rf build/* A.mtx v_batch v_brandes stderr.txt bc_test

