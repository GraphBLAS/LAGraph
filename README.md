[![Build Status](https://github.com/GraphBLAS/LAGraph/workflows/LAGraph%20CI/badge.svg)](https://github.com/GraphBLAS/LAGraph/actions)

# LAGraph

LAGraph is a draft library plus a test harness for collecting algorithms that
use GraphBLAS.

Currently, SuiteSparse:GraphBLAS v6.0.2 or later is required.  However, use the
latest stable release of SuiteSparse:GraphBLAS for best results.
See <https://github.com/DrTimothyAldenDavis/GraphBLAS>

Since LAGraph is a draft, it contains are many draft/experimental codes with
known sub-par performance.  Performance of the best methods is highly sensitive
on which version of SuiteSparse:GraphBLAS is being used, as well.  No one other
than the authors of this code are aware of which methods are the best, and how
to achieve that performance.

Thus, do not benchmark LAGraph on your own without asking the authors first.

To build do the following from the top level directory:
```
cd build
cmake ..
make
make test
```

To install, first build LAGraph and then do (in the build directory):
```
sudo make install
```

To compile with a non-default compiler:
```
CC=gcc-11 cmake ..
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

To run the GAP benchmarks, see the instructions in this file:
```
./src/benchmark/README.md
```

# LAGraph contents

LAGraph contains the following files and folders:

    CMakeLists.txt: a CMake script for compiling:

	( cd build ; cmake .. ; make )

    cmake:  helper scripts for CMake to find GraphBLAS and to provide
        test coverage

    data: small test matrices for the continuous integration tests

    deps: 3rd party dependencies

    doc: documentation

    include: contains the LAGraph.h and LAGraphX.h files

    LICENSE: BSD 2-clause license

    README.md: this file

    src: stable source code for the LAGraph library (LAGraph.h)

        * algorithms: graph algorithms such as BFS, connected components,
            centrality, etc

        * utilities: read/write a graph from a file, etc

    experimental*: draft code under development: (LAGraphX.h)
        do not benchmark without asking the LAGraph authors first

        * algorithms: draft graph algorithms such as Maximal Independent Set

        * utilities: draft utilities

    build: initially empty

# LAGraph and GraphBLAS

To link against GraphBLAS, first install whatever GraphBLAS library you wish to
use.  LAGraph will use -lgraphblas and will include the GraphBLAS.h file
from its installed location.  Alternatively, the CMakeLists.txt script can use
a relative directory:

    ../GraphBLAS: any GraphBLAS implementation.

So that LAGraph and GraphBLAS reside in the same parent folder.  The include
file for GraphBLAS will be assumed to appear in ../GraphBLAS/Include, and the
compiled GraphBLAS library is assumed to appear in ../GraphBLAS/build.  The
CMake should find GraphBLAS, but if you use a GraphBLAS library that uses a
different structure, then edit the CMakeLists.txt file to point to right
location.

# Authors

    Tim Davis, Texas A&M University
    Scott McMillan, SEI, Carnegie Mellon University
    Gabor Szarnyas
    Jinhao Chen, Texas A&M University
    Michel Pelletier, Graphegon
    Scott Kolodziej, Texas A&M University
    Yongzhe Zhang, SOKENDAI, Japan
    Marton Elekes
    Balint Hegyi
    Tim Mattson, Intel
    Mohsen Aznaveh, Texas A&M University
    James Kitchen, Anaconda
    Aydin Buluc, Lawrence Berkeley National Lab
    Janos B. Antal
    Roi Lipman, Redis
    Erik Welch, Anaconda
    Carl Yang
    Tze Meng Low,
    Florentin Dorre

