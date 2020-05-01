%-------------------------------------------------------------------------------
# Test/MatrixExtractKeepDimensions/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
# http://graphblas.org  See LAGraph/LICENSE for license.

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph/Test/MatrixExtractKeepDimensions/matrix_extract_keep_dimensions

# relies on cmake to do the actual build.  Use the CMAKE_OPTIONS argument to
# this Makefile to pass options to cmake.

# Install LAGraph and GraphBLAS before trying to compile LAGraph.

JOBS ?= 1

# build and run test
default: compile
	./build/matrix_extract_keep_dimensions ../../Data/ldbc-directed-example-unweighted.mtx 0
	./build/matrix_extract_keep_dimensions ../../Data/ldbc-undirected-example-unweighted.mtx 1

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
	rm -rf build/* A.mtx matrix_extract_keep_dimensions stderr.txt
