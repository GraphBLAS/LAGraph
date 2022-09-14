LAGraph context and error handling
==================================

The sections below describe a set of functions that manage the LAGraph context
within a user application, and discuss how errors are handled.


LAGraph Context Functions
-------------------------

.. doxygenfunction:: LAGraph_Init

.. doxygenfunction:: LAGr_Init

.. doxygenfunction:: LAGraph_Finalize

.. doxygenfunction:: LAGraph_Version

.. doxygenfunction:: LAGraph_GetNumThreads

.. doxygenfunction:: LAGraph_SetNumThreads

Error handling
--------------

.. doxygendefine:: LAGRAPH_RETURN_VALUES

.. doxygendefine:: LAGRAPH_MSG_LEN

.. doxygendefine:: LAGRAPH_TRY

.. doxygendefine:: GRB_TRY

