Installation
============

LAGraph is available at `<https://github.com/GraphBLAS/LAGraph>`_.
Be sure to check out the reorg branch (for now).
LAGraph requires SuiteSparse:GraphBLAS, available at `<https://github.com/DrTimothyAldenDavis/GraphBLAS>`_.

To compile and install LAGraph, you must first compile and install a recent
version of SuiteSparse:GraphBLAS.  Place LAGraph and GraphBLAS in the same
folder, side-by-side.  Compile and (optionally) install SuiteSparse:GraphBLAS
(see the documentation in SuiteSparse:GraphBLAS for details).
At least on Linux or Mac, if GraphBLAS is not installed system-wide,
LAGraph can find it if GraphBLAS appears in the same folder as LAGraph,
so you do not need system privileges to use GraphBLAS.

Next, in Linux or Mac, run these commands::

    cd LAGraph/build
    cmake ..
    make
    make test

If you have system admin privileges, you can then install it::

    sudo make install

On Windows ... (FIXME: add this)

