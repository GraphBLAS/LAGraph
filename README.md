[![Build Status](https://github.com/GraphBLAS/LAGraph/workflows/LAGraph%20CI/badge.svg)](https://github.com/GraphBLAS/LAGraph/actions)

# LAGraph

LAGraph is a draft library plus a test harness for collecting algorithms that
use the GraphBLAS.

The library can be built and tested with any GraphBLAS library.
Currently, if SuiteSparse:GraphBLAS is used, v6.0.0 or later is required.

Since it's a draft, it contains are many draft/experimental codes with known
sub-par performance.  Performance of the best methods is highly sensitive on
which version of SuiteSparse:GraphBLAS is being used, as well.  No one other
than the authors of this code are aware of which methods are the best, and how
to achieve that performance.

Thus, do not benchmark LAGraph on your own without asking the authors first.

To build do the following from the top level directory:
```
$ mkdir build
$ cd build
$ cmake ..
$ make [-j#]
$ make test
```

To compile with a non-default compiler:
```
$ CC=gcc-11 cmake ..
```

To compile with test coverage:
```
cd build
cmake -DCOVERAGE=1 .. ; make -j8 ; make test_coverage
```

The test coverage of the latest [CI build](https://github.com/GraphBLAS/LAGraph/actions) is deployed to <https://graphblas.github.io/LAGraph/>.

Then open this file in your browser:
```
build/test_coverage/index.html
```

To run the GAP benchmarks, see the instructions in these files:
```
./BENCHMARK_INSTRUCTIONS__README1st.txt
./src/demo/README.md
```

# The rest of this README is outdated and needs to be updated.

LAGraph contains the following files and folders:

    CMakeLists.txt: a CMake script for compiling.  Do not run cmake in this
        top-level directory.  Do "make" here, which does the build in the
        ./build directory:

	( cd build ; cmake .. ; make )

    Doc: documentation

    Include: contains the LAGraph.h file

    LICENSE: BSD 2-clause license

    Makefile: a simple Makefile that relies on CMake to build LAGraph.

    README.md: this file

    Source: stable source code for the LAGraph library:  this is currently
        empty.

        * Algorithms: graph algorithms such as BFS, connected components,
            centrality, etc, will go here

        * Utilities: read/write a graph from a file, etc, will go here...

    Experimental*: draft code under development: do not benchmark without
        asking the LAGraph authors first

        * Algorithms: draft graph algorithms such as BFS, connected components,
            centrality, etc

        * Utilities: draft utilities go here

    Test*: main programs that test LAGraph.  To run the tests, first compile
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

