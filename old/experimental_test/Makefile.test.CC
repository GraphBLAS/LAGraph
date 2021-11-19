#-------------------------------------------------------------------------------
# Test/CC/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
# http://graphblas.org  See LAGraph/LICENSE for license.

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph/Test/CC/cctest,

# relies on cmake to do the actual build.  Use the CMAKE_OPTIONS argument to
# this Makefile to pass options to cmake.

# Install LAGraph and GraphBLAS before trying to compile LAGraph.

JOBS ?= 1

# build and run test
default:
	../build/experimental/CC_test ../data/west0067.mtx
	../build/experimental/CC_test ../data/ldbc-wcc-example.mtx
	../build/experimental/CC_test ../data/tree-example.mtx
