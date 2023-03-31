#[=======================================================================[.rst:
FindLAGraphX
--------

LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
SPDX-License-Identifier: BSD-2-Clause
See additional acknowledgments in the LICENSE file,
or contact permission@sei.cmu.edu for the full terms.

Find the native LAGRAPHX includes and library.

IMPORTED Targets
^^^^^^^^^^^^^^^^

This module defines :prop_tgt:`IMPORTED` target ``LAGRAPHX::LAGRAPHX``, if
LAGRAPHX has been found.

Result Variables
^^^^^^^^^^^^^^^^

This module defines the following variables:

::

  LAGRAPHX_INCLUDE_DIR    - where to find LAGraph.h, etc.
  LAGRAPHX_LIBRARY        - LAGraph library
  LAGRAPHX_LIBRARIES      - List of libraries when using LAGraph.
  LAGRAPHX_FOUND          - True if LAGraph found.

::

Hints
^^^^^

A user may set ``LAGRAPHX_ROOT`` to a LAGraph installation root to tell this
module where to look.

Otherwise, the first place searched is in ../LAGraph, relative to the current
source directory.  That is, if the user application and LAGraph reside in the
same parent folder, side-by-side, and if it contains LAGraph/include/LAGraph.h
file and LAGraph/build/lib/liblagraph.so (or dylib, etc), then that version is
used.  This takes precedence over the system-wide installation of LAGraph,
which might be an older version.  This method gives the user the ability to
compile LAGraph with their own copy of LAGraph, ignoring the system-wide
version.

#]=======================================================================]

# "include" for LAGraphX
find_path(
  LAGRAPHX_INCLUDE_DIR
  NAMES LAGraphX.h
  HINTS ${LAGRAPHX_ROOT}
  HINTS ENV LAGRAPHX_ROOT
  HINTS ${CMAKE_SOURCE_DIR}/../LAGraph
  PATHS LAGRAPHX_ROOT ENV LAGRAPHX_ROOT
  PATH_SUFFIXES include Include
  )

# "build/lib" for LAGraph
message ( STATUS "Look in " ${CMAKE_SOURCE_DIR}/../LAGraph/build )
find_library(
  LAGRAPHX_LIBRARY
  NAMES lagraphx
  HINTS ${LAGRAPHX_ROOT}
  HINTS ENV LAGRAPHX_ROOT
  HINTS ${CMAKE_SOURCE_DIR}/../LAGraph
  PATHS LAGRAPHX_ROOT ENV LAGRAPHX_ROOT
  PATH_SUFFIXES lib build
  )

# get version of .so using REALPATH
get_filename_component(LAGRAPHX_LIBRARY ${LAGRAPHX_LIBRARY} REALPATH)
string(
  REGEX MATCH "[0-9]+.[0-9]+.[0-9]+"
  LAGRAPHX_VERSION
  ${LAGRAPHX_LIBRARY}
  )
set(LAGRAPHX_LIBRARIES ${LAGRAPHX_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  LAGraphX
  REQUIRED_VARS LAGRAPHX_LIBRARIES LAGRAPHX_INCLUDE_DIR
  VERSION_VAR LAGRAPHX_VERSION
  )

mark_as_advanced(
  LAGRAPHX_INCLUDE_DIR
  LAGRAPHX_LIBRARY
  LAGRAPHX_LIBRARIES
  )

if ( LAGRAPHX_FOUND )
    message ( STATUS "LAGraphX include dir: " ${LAGRAPHX_INCLUDE_DIR} )
    message ( STATUS "LAGraphX library:     " ${LAGRAPHX_LIBRARY} )
    message ( STATUS "LAGraphX version:     " ${LAGRAPHX_VERSION} )
else ( )
    message ( STATUS "LAGraphX not found" )
endif ( )

