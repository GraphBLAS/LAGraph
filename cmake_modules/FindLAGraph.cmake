#[=======================================================================[.rst:
FindLAGraph
--------

LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
SPDX-License-Identifier: BSD-2-Clause
See additional acknowledgments in the LICENSE file,
or contact permission@sei.cmu.edu for the full terms.

Find the native LAGRAPH includes and library.

IMPORTED Targets
^^^^^^^^^^^^^^^^

This module defines :prop_tgt:`IMPORTED` target ``LAGRAPH::LAGRAPH``, if
LAGRAPH has been found.

Result Variables
^^^^^^^^^^^^^^^^

This module defines the following variables:

::

  LAGRAPH_INCLUDE_DIR    - where to find LAGraph.h, etc.
  LAGRAPH_LIBRARY        - LAGraph library
  LAGRAPH_LIBRARIES      - List of libraries when using LAGraph.
  LAGRAPH_FOUND          - True if LAGraph found.

::

Hints
^^^^^

A user may set ``LAGRAPH_ROOT`` to a LAGraph installation root to tell this
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

# "include" for LAGraph
find_path(
  LAGRAPH_INCLUDE_DIR
  NAMES LAGraph.h
  HINTS ${CMAKE_SOURCE_DIR}/../LAGraph
  PATHS LAGRAPH_ROOT ENV LAGRAPH_ROOT
  PATH_SUFFIXES include Include
  )

# "build/lib" for LAGraph
message ( STATUS "Look in " ${CMAKE_SOURCE_DIR}/../LAGraph/build )
find_library(
  LAGRAPH_LIBRARY
  NAMES lagraph
  HINTS ${CMAKE_SOURCE_DIR}/../LAGraph
  PATHS LAGRAPH_ROOT ENV LAGRAPH_ROOT
  PATH_SUFFIXES lib build
  )

# get version of .so using REALPATH
get_filename_component(LAGRAPH_LIBRARY ${LAGRAPH_LIBRARY} REALPATH)
string(
  REGEX MATCH "[0-9]+.[0-9]+.[0-9]+"
  LAGRAPH_VERSION
  ${LAGRAPH_LIBRARY}
  )
set(LAGRAPH_LIBRARIES ${LAGRAPH_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  LAGraph
  REQUIRED_VARS LAGRAPH_LIBRARIES LAGRAPH_INCLUDE_DIR
  VERSION_VAR LAGRAPH_VERSION
  )

mark_as_advanced(
  LAGRAPH_INCLUDE_DIR
  LAGRAPH_LIBRARY
  LAGRAPH_LIBRARIES
  )

if ( LAGRAPH_FOUND )
    message ( STATUS "LAGraph include dir: " ${LAGRAPH_INCLUDE_DIR} )
    message ( STATUS "LAGraph library:     " ${LAGRAPH_LIBRARY} )
    message ( STATUS "LAGraph version:     " ${LAGRAPH_VERSION} )
else ( )
    message ( STATUS "LAGraph not found" )
endif ( )

