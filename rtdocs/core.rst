FIXME TITLE
============

FIXME: intro.  startup

Program Functions
-----------------

.. doxygenfunction:: LAGraph_Init

.. doxygenfunction:: LAGr_Init

.. doxygenfunction:: LAGraph_Finalize

.. doxygenfunction:: LAGraph_Version

.. doxygenfunction:: LAGraph_GetNumThreads

.. doxygenfunction:: LAGraph_SetNumThreads

Error handling
--------------

FIXME: overview here.

.. doxygendefine:: LAGRAPH_RETURN_VALUES


    .. doxygendefine:: LAGRAPH_RETURN_VALUES

    .. doxygendefine:: LAGRAPH_INVALID_GRAPH

    .. doxygendefine:: LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED

    .. doxygendefine:: LAGRAPH_IO_ERROR

    .. doxygendefine:: LAGRAPH_NOT_CACHED

    .. doxygendefine:: LAGRAPH_NO_SELF_EDGES_ALLOWED

    .. doxygendefine:: LAGRAPH_CONVERGENCE_FAILURE

    .. doxygendefine:: LAGRAPH_CACHE_NOT_NEEDED

.. doxygendefine::  LAGRAPH_MSG_LEN

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


