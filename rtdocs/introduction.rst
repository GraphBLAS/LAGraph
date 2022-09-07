Introduction
============

FIXME: discuss the GraphBLAS idea:  matrix operations with semirings.
GraphBLAS, however, doesn't have any graph algorithms.

discuss the LAGraph idea:  graph algorithms using GraphBLAS,
with several additional data structures.  In particular the LAGraph_Graph.
Discuss its properties, and cached properties, and the differences
between Basic (``LAGraph_*``) and Advanced (``LAGr_*``).

LAGraph and GraphBLAS work together; you can have a mix of
LAGraph_Graph objects and GraphBLAS objects together in the same
program.  LAGraph methods often return GrB_Matrix or GrB_Vector
results.  An LAGraph_Graph can be modified by the user application
directly (it is not an opaque object), by including GrB_Matrix or
GrB_Vector objects into it.

