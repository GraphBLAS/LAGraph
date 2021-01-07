#-------------------------------------------------------------------------------
# Test2/BFS/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph/Test2/BFS/test_bfs

# Install LAGraph and GraphBLAS before trying to compile LAGraph.
JOBS ?= 1

# build and run test
default: compile
	./build/test_bfs < ../../Matrices/cover.mtx
	./build/test_bfs < ../../Matrices/west0067.mtx
	./build/test_bfs < ../../Matrices/jagmesh7.mtx

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
	rm -rf build/* test_bfs

