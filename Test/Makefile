#-------------------------------------------------------------------------------
# LAGraph/Test/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
# http://graphblas.org  See LAGraph/LICENSE for license.

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph, relies on cmake to do the actual build.  Use
# the CMAKE_OPTIONS argument to this Makefile to pass options to cmake.

# Install GraphBLAS before trying to compile LAGraph.

JOBS ?= 1

# build all the tests
default:
	( cd MatrixMarket ; $(MAKE) --jobs=$(JOBS) )
	( cd BFS          ; $(MAKE) --jobs=$(JOBS) )

clean: distclean

purge: distclean

# remove all files not in the distribution
distclean:
	( cd MatrixMarket ; $(MAKE) distclean )
	( cd BFS          ; $(MAKE) distclean )

