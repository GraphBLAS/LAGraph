#-------------------------------------------------------------------------------
# Test2/ConnectedComponents/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph/Test2/ConnectedComponents/test_cc

# Install LAGraph and GraphBLAS before trying to compile this test.
JOBS ?= 1

# build and run test
default: compile
	./build/test_cc ../../Matrices/cover.mtx
	./build/test_cc ../../Matrices/west0067.mtx
	./build/test_cc ../../Matrices/jagmesh7.mtx
	./build/test_cc ../../Matrices/ldbc-wcc-example.mtx
	./build/test_cc ../../Matrices/tree-example.mtx

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
	rm -rf build/* test_cc

