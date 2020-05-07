#-------------------------------------------------------------------------------
# Test/MatrixExtractKeepDimensions/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
# http://graphblas.org  See LAGraph/LICENSE for license.

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph/Test/MatrixExtractKeepDimensions/matrix_extract_keep_dimensions_test

# relies on cmake to do the actual build.  Use the CMAKE_OPTIONS argument to
# this Makefile to pass options to cmake.

# Install LAGraph and GraphBLAS before trying to compile LAGraph.

JOBS ?= 1

# build and run test
default: compile
	./build/matrix_extract_keep_dimensions_test ../../Data/ldbc-directed-example-unweighted.mtx
	./build/matrix_extract_keep_dimensions_test ../../Data/ldbc-undirected-example-unweighted.mtx

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
	rm -rf build/* matrix_extract_keep_dimensions_test
