#-------------------------------------------------------------------------------
# LAGraph/CMakeLists.txt:  cmake script for LAGraph
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
#
# See additional acknowledgments in the LICENSE file,
# or contact permission@sei.cmu.edu for the full terms.

#-------------------------------------------------------------------------------

# CMakeLists.txt: instructions for cmake to build LAGraph.  An ANSI C11
# compiler is required.  First, install any GraphBLAS library.  Alternatively,
# use ../GraphBLAS (see comments below).
#
# To compile the LAGraph library, all of its tests and demos:
#
#   mkdir buld
#   cd build
#   cmake ..
#   make
#
# If that fails for any reason, make sure your compiler supports ANSI C11.  You
# could try changing your compiler for the cmake command, for example:
#
#   CC=icc cmake ..
#   CC=xlc cmake ..
#   CC=gcc cmake ..

#-------------------------------------------------------------------------------
# get the version
#-------------------------------------------------------------------------------

# cmake 3.13 is preferred.
cmake_minimum_required ( VERSION 3.1 )

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

if ( CMAKE_VERSION VERSION_GREATER "3.0" )
    cmake_policy ( SET CMP0042 NEW )
    cmake_policy ( SET CMP0048 NEW )
endif ( )

# version of LAGraph
set ( LAGraph_DATE "April 1, 2021" )
set ( LAGraph_VERSION_MAJOR 0 )
set ( LAGraph_VERSION_MINOR 9 )
set ( LAGraph_VERSION_SUB   9 )

project ( lagraph
          VERSION "${LAGraph_VERSION_MAJOR}.${LAGraph_VERSION_MINOR}.${LAGraph_VERSION_SUB}" )

#-------------------------------------------------------------------------------
# determine build type
#-------------------------------------------------------------------------------

# For development only, not for end-users:
# set ( CMAKE_BUILD_TYPE Debug )

if ( NOT CMAKE_BUILD_TYPE )
    set ( CMAKE_BUILD_TYPE Release )
endif ( )

#-------------------------------------------------------------------------------
# Find the GraphBLAS installation
#-------------------------------------------------------------------------------

# If GraphBLAS is not in a standard installation location, either
#     export GRAPHBLAS_ROOT <path>
# or
#     GRAPHBLAS_ROOT=<path> cmake ..
# or uncomment the next line:
set ( ENV{GRAPHBLAS_ROOT} ${CMAKE_SOURCE_DIR}/../GraphBLAS )

find_package (GraphBLAS REQUIRED MODULE)
include_directories (${GRAPHBLAS_INCLUDE_DIR})


#-------------------------------------------------------------------------------
# Is this really necessary?
#-------------------------------------------------------------------------------

include ( GNUInstallDirs )

#-------------------------------------------------------------------------------
# determine what user threading model to use
#-------------------------------------------------------------------------------

include ( FindOpenMP  )
include ( FindThreads )

#-------------------------------------------------------------------------------
# report status
#-------------------------------------------------------------------------------

message ( STATUS "CMAKE build type:          " ${CMAKE_BUILD_TYPE} )

if ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    message ( STATUS "CMAKE C Flags debug:       " ${CMAKE_C_FLAGS_DEBUG} )
else ( )
    message ( STATUS "CMAKE C Flags release:     " ${CMAKE_C_FLAGS_RELEASE} )
endif ( )

message ( STATUS "CMAKE compiler ID:         " ${CMAKE_C_COMPILER_ID} )
message ( STATUS "CMAKE thread library:      " ${CMAKE_THREAD_LIBS_INIT} )
message ( STATUS "CMAKE have pthreads:       " ${CMAKE_USE_PTHREADS_INIT}  )
message ( STATUS "CMAKE have Win32 pthreads: " ${CMAKE_USE_WIN32_THREADS_INIT} )
message ( STATUS "CMAKE have OpenMP:         " ${OPENMP_FOUND} )

#-------------------------------------------------------------------------------
# include directories for lagraph library
#-------------------------------------------------------------------------------

#set ( CMAKE_INCLUDE_CURRENT_DIR ON )

include_directories ( ${CMAKE_SOURCE_DIR}/include
                      ${CMAKE_SOURCE_DIR}/src/utility
                      ${CMAKE_SOURCE_DIR}/test/include )
#                      ${CMAKE_SOURCE_DIR}/experimental/include
#                      ${CMAKE_SOURCE_DIR}/experimental/utility
#                      ${CMAKE_SOURCE_DIR}/experimental/utility/Template )

#-------------------------------------------------------------------------------
# compiler options
#-------------------------------------------------------------------------------

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

# check which compiler is being used.  If you need to make
# compiler-specific modifications, here is the place to do it.
if ( "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    # cmake 2.8 workaround: gcc needs to be told to do ANSI C11.
    # cmake 3.0 doesn't have this problem.
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -std=c11 -lm -Wno-pragmas " )
    # check all warnings:
#   set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -Wextra -Wpedantic -Werror " )
    # set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -g " )
    #if ( CMAKE_C_COMPILER_VERSION VERSION_LESS 4.9 )
    #    message ( FATAL_ERROR "gcc version must be at least 4.9" )
    #endif ( )
elseif ( "${CMAKE_C_COMPILER_ID}" STREQUAL "Intel" )
    # options for icc: also needs -std=c11
    # note that -g can be used, for VTune.  Comment out the following line
    # to compile without -g.
    # set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -g " )
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -std=c11 " )
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -qopt-malloc-options=3" )
    # check all warnings:
#   set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -w3 " )
    #if ( CMAKE_C_COMPILER_VERSION VERSION_LESS 18.0 )
    #    message ( FATAL_ERROR "icc version must be at least 18.0" )
    #endif ( )
elseif ( "${CMAKE_C_COMPILER_ID}" STREQUAL "Clang" )
    # options for clang
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -O3 -std=c11 " )
    #if ( CMAKE_C_COMPILER_VERSION VERSION_LESS 3.3 )
    #    message ( FATAL_ERROR "clang version must be at least 3.3" )
    #endif ( )
elseif ( "${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC" )
    # options for MicroSoft Visual Studio
elseif ( "${CMAKE_C_COMPILER_ID}" STREQUAL "PGI" )
    # options for PGI pgcc compiler
    set ( CMAKE_C_FLAGS    "${CMAKE_C_FLAGS} -acc -Minfo=accel -Mcuda -Mnoopenmp -noswitcherror -c11 -lm -fPIC " )
    set ( CMAKE_CXX_FLAGS  "${CMAKE_C_FLAGS} -acc -Minfo=accel -Mcuda -Mnoopenmp -D__GCC_ATOMIC_TEST_AND_SET_TRUEVAL=1 -noswitcherror --c++11 -lm -fPIC " )
endif ( )

if ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}" )
else ( )
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}" )
endif ( )

#-------------------------------------------------------------------------------
# select the threading library 
#-------------------------------------------------------------------------------

if ( USER_OPENMP )
    # user insists on OpenMP synchronization inside LAGraph
    message ( STATUS "cmake -DUSER_OPENMP=1: insisting on using OpenMP" )
    set ( USE_OPENMP true )
elseif ( USER_POSIX )
    # user insists on POSIX synchronization inside LAGraph
    message ( STATUS "cmake -DUSER_POSIX=1: insisting on using POSIX" )
    set ( USE_POSIX true )
elseif ( USER_NONE )
    message ( STATUS "cmake -DUSER_NONE=1: insisting on using no threading" )
    set ( USE_NONE true )
else ( )
    # default: automatic selection
    if ( OPENMP_FOUND )
        set ( USE_OPENMP true )
        message ( STATUS "Automatic selection of synchronization method: openmp" )
    elseif ( CMAKE_USE_PTHREADS_INIT )
        set ( USE_POSIX true )
        message ( STATUS "Automatic selection of synchronization method: pthreads" )
    endif ( )
endif ( )

# add the threading library
# if ( USE_OPENMP )
#     # use OpenMP
#     message ( STATUS "Using OpenMP to synchronize user threads" )
#     target_link_libraries ( lagraph        m ${OpenMP_C_LIBRARIES} )
#     target_link_libraries ( lagraph_static m ${OpenMP_C_LIBRARIES} )
#     set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DUSER_OPENMP_THREADS " )
# elseif ( USE_POSIX )
#     # use POSIX
#     message ( STATUS "Using POSIX pthreads to synchronize user threads" )
#     target_link_libraries ( lagraph        m )
#     target_link_libraries ( lagraph_static m )
#     set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -pthread -DUSER_POSIX_THREADS " )
# else ( )
#     # use no threading at all
#     message ( WARNING "No support for user threads; LAGraph will not be thread-safe" )
#     target_link_libraries ( lagraph        m )
#     target_link_libraries ( lagraph_static m )
#     set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -DUSER_NO_THREADS " )
# endif ( )

# if ( CMAKE_THREAD_LIBS_INIT )
#     target_link_libraries ( lagraph        ${CMAKE_THREAD_LIBS_INIT} )
#     target_link_libraries ( lagraph_static ${CMAKE_THREAD_LIBS_INIT} )
# endif ( )

# target_link_libraries ( lagraph        ${GRAPHBLAS_LIBRARIES} )
# target_link_libraries ( lagraph_static ${GRAPHBLAS_LIBRARIES} )

# if ( CMAKE_USE_PTHREADS_INIT )
#     set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -DHAVE_PTHREADS " )
# endif ( )

# if ( CMAKE_USE_WIN32_THREADS_INIT )
#     set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -DHAVE_WINDOWS_THREADS " )
# endif ( )

#-------------------------------------------------------------------------------
# print final C flags
#-------------------------------------------------------------------------------

message ( STATUS "CMAKE C flags: " ${CMAKE_C_FLAGS} )

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# CTEST cannot find the binaries when in this location
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

enable_testing()

add_subdirectory(src)
#add_subdirectory(experimental)

