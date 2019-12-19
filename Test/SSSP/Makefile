#-------------------------------------------------------------------------------
# Test/SSSP/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
# http://graphblas.org  See LAGraph/LICENSE for license.

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph/Test/SSSP/sssp_test,

# relies on cmake to do the actual build.  Use the CMAKE_OPTIONS argument to
# this Makefile to pass options to cmake.

# Install LAGraph and GraphBLAS before trying to compile.

JOBS ?= 1

# build and run test
default: compile
	#./build/sssp_test < jagmesh7.mtx 3 > jagmesh7.out #2>stderr.txt
	#./build/sssp_check jagmesh7.out jagmesh7.check
	./build/sssp_test < cover.mtx 3 #> cover.out #2>stderr.txt
	./build/sssp_test /raid/hyper/gapbs/GAP/GAP-kron/GAP-kron.grb /raid/hyper/gapbs/GAP/GAP-kron/GAP-kron_sources.mtx 2 > kron.out
	./build/sssp_test /raid/hyper/gapbs/GAP/GAP-road/GAP-road.grb /raid/hyper/gapbs/GAP/GAP-road/GAP-road_sources.mtx 50000 > road.out
	./build/sssp_test /raid/hyper/gapbs/GAP/GAP-twitter/GAP-twitter.grb /raid/hyper/gapbs/GAP/GAP-twitter/GAP-twitter_sources.mtx 3 > twitter.out
	./build/sssp_test /raid/hyper/gapbs/GAP/GAP-urand/GAP-urand.grb /raid/hyper/gapbs/GAP/GAP-urand/GAP-urand_sources.mtx 3 > urand.out
	./build/sssp_test /raid/hyper/gapbs/GAP/GAP-web/GAP-web.grb /raid/hyper/gapbs/GAP/GAP-web/GAP-web_sources.mtx 3 > web.out

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
	rm -rf build/* stderr.txt sssp_test sssp_exe *.out

