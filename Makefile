#-------------------------------------------------------------------------------
# LAGraph/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
#
# See additional acknowledgments in the LICENSE file,
# or contact permission@sei.cmu.edu for the full terms.

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph, relies on cmake to do the actual build.  Use
# the CMAKE_OPTIONS argument to this Makefile to pass options to cmake.

# to compile with gcc, use:

#   make CC=gcc-11 CXX=g++-11

JOBS ?= 8

default: library

# just build the dynamic library, not the demos
library:
	( cd build ; cmake $(CMAKE_OPTIONS) .. ; $(MAKE) --jobs=$(JOBS) )

# run the tests
test: library
	( cd build ; make test )

# just run cmake; do not compile
docmake:
	( cd build ; cmake $(CMAKE_OPTIONS) .. ; )

# installs LAGraph to the install location defined by cmake, usually
# /usr/local/lib and /usr/local/include
install:
	( cd build ; $(MAKE) install )

# remove any installed libraries and #include files
uninstall:
	- xargs rm < build/install_manifest.txt

clean: distclean

purge: distclean

# remove all files not in the distribution
distclean:
	rm -rf build/* 

