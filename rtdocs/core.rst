Core Objects
============

Graph Object
------------

The fundamental object in LAGraph is the `LAGraph_Graph`.

.. doxygentypedef:: LAGraph_Graph

.. doxygenstruct:: LAGraph_Graph_struct
    :members:
    :undoc-members:

Program Functions
-----------------

.. doxygenfunction:: LAGraph_Version

.. doxygenfunction:: LAGraph_Init

.. doxygenfunction:: LAGr_Init

.. doxygenfunction:: LAGraph_Finalize

.. doxygenfunction:: LAGraph_GetNumThreads

.. doxygenfunction:: LAGraph_SetNumThreads

Graph Functions
---------------

.. doxygenfunction:: LAGraph_New

.. doxygenfunction:: LAGraph_Delete

.. doxygenfunction:: LAGraph_DeleteCached

.. doxygenfunction:: LAGraph_Cached_AT

.. doxygenfunction:: LAGraph_Cached_IsSymmetricStructure

.. doxygenfunction:: LAGraph_Cached_OutDegree

.. doxygenfunction:: LAGraph_Cached_InDegree

.. doxygenfunction:: LAGraph_Cached_NSelfEdges

.. doxygenfunction:: LAGraph_Cached_EMin

.. doxygenfunction:: LAGraph_Cached_EMax

.. doxygenfunction:: LAGraph_DeleteSelfEdges

.. doxygenfunction:: LAGraph_CheckGraph

Input/Output Functions
----------------------

.. doxygenfunction:: LAGraph_MMRead

.. doxygenfunction:: LAGraph_MMWrite

Error handling
--------------

FIXME: Discuss the msg string, and return values.

.. doxygendefine:: LAGRAPH_RETURN_VALUES

.. doxygendefine:: LAGRAPH_TRY

.. doxygendefine:: LAGRAPH_INVALID_GRAPH

.. doxygendefine:: GRB_TRY

Enums
-----

.. doxygenenum:: LAGraph_Kind

.. doxygenenum:: LAGraph_Boolean

.. doxygenenum:: LAGraph_State

Pre-defined semirings
---------------------

LAGraph adds the following pre-defined semirings.  They are created
by `LAGr_Init` or `LAGraph_Init`, and freed by `LAGraph_Finalize`.

.. doxygenvariable:: LAGraph_plus_first_int8


