#-------------------------------------------------------------------------------
# test/PageRank/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
#
# See additional acknowledgments in the LICENSE file,
# or contact permission@sei.cmu.edu for the full terms.

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph/test/PageRank/test_gappagerank

# Install LAGraph and GraphBLAS before trying to compile this test.
JOBS ?= 1

# build and run test
default: compile
	./build/test_gappagerank < ../data/cover.mtx
	./build/test_gappagerank < ../data/west0067.mtx
	./build/test_gappagerank < ../data/jagmesh7.mtx

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
	rm -rf build/* test_gappagerank

