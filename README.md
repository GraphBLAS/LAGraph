[![Build Status](https://travis-ci.org/GraphBLAS/LAGraph.svg?branch=master)](https://travis-ci.org/GraphBLAS/LAGraph)

# LAGraph

This is a library plus a test harness for collecting algorithms that use the
GraphBLAS.  It contains the following files and folders:

    CMakeLists.txt: a CMake script for compiling.  Do not run cmake in this
        top-level directory.  Do "make" here, which does the build in the
        ./build directory:

	( cd build ; cmake .. ; make )

    Doc: documentation

    Include: contains the LAGraph.h file

    LICENSE: BSD 2-clause license

    Makefile: a simple Makefile that relies on CMake to build LAGraph.

    README.md: this file

    Source: stable source code for the LAGraph library

        * Algorithms: graph algorithms such as BFS, connected components,
            centrality, etc, will go here

        * Utilities: read/write a graph from a file, etc, will go here...

    Experimental: draft code under development: do not benchmark without
        asking the LAGraph authors first

        * Algorithms: graph algorithms such as BFS, connected components,
            centrality, etc

        * Utilities: read/write a graph from a file, etc

    Test: main programs that test LAGraph.  To run the tests, first compile
        GraphBLAS and LAGraph, and then do "make tests" in this directory.

    build: initially empty

To link against GraphBLAS, first install whatever GraphBLAS library you wish to
use.  LAGraph will use -lgraphblas and will include the GraphBLAS.h file
from its installed location.  Alternatively, the CMakeLists.txt script can use
a relative directory:

    ../GraphBLAS: any GraphBLAS implementation.

So that LAGraph and GraphBLAS reside in the same parent folder.  The include
file for GraphBLAS will be assumed to appear in ../GraphBLAS/Include, and the
compiled GraphBLAS library is assumed to appear in ../GraphBLAS/build.  If you
use a GraphBLAS library that uses a different structure, then edit the
CMakeLists.txt file to point to right location.

Authors: (... list them here)

