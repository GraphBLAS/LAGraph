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
  LAGRAPH_LIBRARY        - dynamic LAGraph library
  LAGRAPH_STATIC         - static LAGraph library
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
  HINTS ${LAGRAPH_ROOT}
  HINTS ENV LAGRAPH_ROOT
  HINTS ${CMAKE_SOURCE_DIR}/../LAGraph
  PATH_SUFFIXES include Include
  )

# dynamic LAGraph library
find_library(
  LAGRAPH_LIBRARY
  NAMES lagraph
  HINTS ${LAGRAPH_ROOT}
  HINTS ENV LAGRAPH_ROOT
  HINTS ${CMAKE_SOURCE_DIR}/../LAGraph
  PATH_SUFFIXES lib build
  )

if ( MSVC )
    set ( STATIC_SUFFIX .lib )
else ( )
    set ( STATIC_SUFFIX .a )
endif ( )

# static LAGraph library
set ( save ${CMAKE_FIND_LIBRARY_SUFFIXES} )
set ( CMAKE_FIND_LIBRARY_SUFFIXES ${STATIC_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES} )
find_library ( LAGRAPH_STATIC
  NAMES lagraph
  HINTS ${LAGRAPH_ROOT}
  HINTS ENV LAGRAPH_ROOT
  HINTS ${CMAKE_SOURCE_DIR}/../LAGraph
  PATH_SUFFIXES lib build
  )
set ( CMAKE_FIND_LIBRARY_SUFFIXES ${save} )

# get version of the library from the dynamic library name
get_filename_component(LAGRAPH_LIBRARY ${LAGRAPH_LIBRARY} REALPATH)
string(
  REGEX MATCH "[0-9]+.[0-9]+.[0-9]+"
  LAGRAPH_VERSION
  ${LAGRAPH_LIBRARY}
  )

if ( NOT LAGRAPH_VERSION )
    # if the version does not appear in the filename, read the include file.
    # The LAGraph.h file includes the following lines:
    #   #define LAGRAPH_VERSION_MAJOR  1
    #   #define LAGRAPH_VERSION_MINOR  0
    #   #define LAGRAPH_VERSION_UPDATE 1
    file ( STRINGS ${LAGRAPH_INCLUDE_DIR}/LAGraph.h LAGRAPH_VER_MAJOR_STRING
        REGEX "define LAGRAPH_VERSION_MAJOR " )
    file ( STRINGS ${LAGRAPH_INCLUDE_DIR}/LAGraph.h LAGRAPH_VER_MINOR_STRING
        REGEX "define LAGRAPH_VERSION_MINOR " )
    file ( STRINGS ${LAGRAPH_INCLUDE_DIR}/LAGraph.h LAGRAPH_VER_UPDATE_STRING
        REGEX "define LAGRAPH_VERSION_UPDATE " )
    message ( STATUS "major from LAGraph.h:  ${LAGRAPH_VER_MAJOR_STRING}" )
    message ( STATUS "minor from LAGraph.h:  ${LAGRAPH_VER_MINOR_STRING}" )
    message ( STATUS "update from LAGraph.h: ${LAGRAPH_VER_UPDATE_STRING}" )
    string ( REGEX MATCH "[0-9]+" LAGRAPH_VER_MAJOR ${LAGRAPH_VER_MAJOR_STRING} )
    string ( REGEX MATCH "[0-9]+" LAGRAPH_VER_MINOR ${LAGRAPH_VER_MINOR_STRING} )
    string ( REGEX MATCH "[0-9]+" LAGRAPH_VER_UPDATE ${LAGRAPH_VER_UPDATE_STRING} )
    set ( LAGRAPH_VERSION "${LAGRAPH_VER_MAJOR}.${LAGRAPH_VER_MINOR}.${LAGRAPH_VER_UPDATE}")
endif ( )

set ( LAGRAPH_LIBRARIES ${LAGRAPH_LIBRARY} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  LAGraph
  REQUIRED_VARS LAGRAPH_LIBRARIES LAGRAPH_INCLUDE_DIR
  VERSION_VAR LAGRAPH_VERSION
  )

mark_as_advanced(
  LAGRAPH_INCLUDE_DIR
  LAGRAPH_LIBRARY
  LAGRAPH_STATIC
  LAGRAPH_LIBRARIES
  )

if ( LAGRAPH_FOUND )
    message ( STATUS "LAGraph version: " ${LAGRAPH_VERSION} )
    message ( STATUS "LAGraph include: " ${LAGRAPH_INCLUDE_DIR} )
    message ( STATUS "LAGraph library: " ${LAGRAPH_LIBRARY} )
    message ( STATUS "LAGraph static:: " ${LAGRAPH_STATIC} )
else ( )
    message ( STATUS "LAGraph not found" )
endif ( )

