#-------------------------------------------------------------------------------
# Test2/TriangleCount/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph/Test2/TriangleCount/test_tc

# Install LAGraph and GraphBLAS before trying to compile this test.
JOBS ?= 1

# build and run test
default: compile
	./build/test_tc < ../../Matrices/cover.mtx
	./build/test_tc < ../../Matrices/west0067.mtx
	./build/test_tc < ../../Matrices/jagmesh7.mtx

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
	rm -rf build/* test_tc

