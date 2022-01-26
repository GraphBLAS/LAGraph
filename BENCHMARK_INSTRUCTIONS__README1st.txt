FIXME: remove this file when reaching v1.0

################################################################################
################################################################################
##                                                                            ##
##  To benchmark LAGraph: please contact Tim Davis at davis@tamu.edu first.   ##
##                                                                            ##
################################################################################
################################################################################

LAGraph is a draft package, and its performance is not yet stable.  It includes
many draft algorithms that are sometimes posted on github in debug mode, or
with known suboptimal performance.  We ask that you not benchmark LAGraph on
your own without contacting the authors to make sure you have the right
version, and the right version of SuiteSparse:GraphBLAS to go with it.

If you run in vanilla mode, by compiling LAGraph with

    cmake -DLG_VANILLA=1 ..

Then performance can be quite low since in this case LAGraph does not use
any SuiteSparse:GraphBLAS GxB* extensions.  We are still developing the
pure GrB* implementations of these algorithms.

However, assuming things are stable, follow the instructions in
LAGraph/src/demo/README.md to run the GAP benchmarks.

