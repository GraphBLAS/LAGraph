#[=======================================================================[.rst:
FindGraphBLAS
--------

Find the native GRAPHBLAS includes and library.

IMPORTED Targets
^^^^^^^^^^^^^^^^

This module defines :prop_tgt:`IMPORTED` target ``GRAPHBLAS::GRAPHBLAS``, if
GRAPHBLAS has been found.

Result Variables
^^^^^^^^^^^^^^^^

This module defines the following variables:

::

  GRAPHBLAS_INCLUDE_DIRS   - where to find GraphBLAS.h, etc.
  GRAPHBLAS_LIBRARIES      - List of libraries when using GraphBLAS.
  GRAPHBLAS_FOUND          - True if GraphBLAS found.

::

Hints
^^^^^

A user may set ``GRAPHBLAS_ROOT`` to a GraphBLAS installation root to tell this
module where to look.
#]=======================================================================]


# NB: this is built around assumptions about one particular GraphBLAS
# installation. As other installations become available changes to
# this will likely be required.

# "Include" for a well known installation with unconventional naming.
find_path(
  GRAPHBLAS_INCLUDE_DIR
  NAMES GraphBLAS.h
  PATHS GRAPHBLAS_ROOT ENV GRAPHBLAS_ROOT
  PATH_SUFFIXES include Include
  )

# "build" for a well known installation with unconventional naming.
find_library(
  GRAPHBLAS_LIBRARY
  NAMES graphblas
  PATHS GRAPHBLAS_ROOT ENV GRAPHBLAS_ROOT
  PATH_SUFFIXES lib build
  )

# get version of .so using REALPATH
get_filename_component(GRAPHBLAS_LIBRARY ${GRAPHBLAS_LIBRARY} REALPATH)
string(
  REGEX MATCH "[0-9]+.[0-9]+.[0-9]+"
  GRAPHBLAS_VERSION
  ${GRAPHBLAS_LIBRARY}
  )
set(GRAPHBLAS_LIBRARIES ${GRAPHBLAS_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  GraphBLAS
  REQUIRED_VARS GRAPHBLAS_LIBRARIES GRAPHBLAS_INCLUDE_DIR
  VERSION_VAR GRAPHBLAS_VERSION
  )

mark_as_advanced(
  GRAPHBLAS_INCLUDE_DIR
  GRAPHBLAS_LIBRARY
  GRAPHBLAS_LIBRARIES
  )
