#-------------------------------------------------------------------------------
# LAGraph/CMakeLists.txt:  cmake script for LAGraph
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
# See additional acknowledgments in the LICENSE file,
# or contact permission@sei.cmu.edu for the full terms.

#-------------------------------------------------------------------------------

# CMakeLists.txt: instructions for cmake to build LAGraph.  An ANSI C11
# compiler is required.  First, install any GraphBLAS library.  Alternatively,
# use ../GraphBLAS (see comments below).
#
# To compile the LAGraph library and its tests and benchmarks, and run the
# tests:
#
#   cd build
#   cmake ..
#   make
#   make test
#
# If that fails for any reason, make sure your compiler supports ANSI C11.  You
# could try changing your compiler for the cmake command, for example:
#
#   CC=icc cmake ..
#   CC=xlc cmake ..
#   CC=gcc cmake ..

#-------------------------------------------------------------------------------
# link directories for lagraph library
#-------------------------------------------------------------------------------

# ensure the tests look for liblagraph.so and liblagraphx.so in the build
# folder, before looking elsewhere
# link_directories ( BEFORE ${CMAKE_BINARY_DIR} )
# link_directories ( ${CMAKE_BINARY_DIR} )
# message ( STATUS "binary dir:     : " ${CMAKE_BINARY_DIR} )
# message ( STATUS "link directories: " ${CMAKE_LINK_DIRECTORIES} )
# message ( STATUS "library path:   : " ${CMAKE_LIBRARY_PATH} )

#-------------------------------------------------------------------------------
# get the version
#-------------------------------------------------------------------------------

cmake_minimum_required ( VERSION 3.13 )

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

if ( CMAKE_VERSION VERSION_GREATER "3.0" )
    cmake_policy ( SET CMP0042 NEW )
    cmake_policy ( SET CMP0048 NEW )
endif ( )
set ( CMAKE_MACOSX_RPATH TRUE )

# version of LAGraph
# FIXME: use config to create include/LAGraph.h from these definitions
set ( LAGraph_DATE "Mar 11, 2022" )
set ( LAGraph_VERSION_MAJOR 0 )
set ( LAGraph_VERSION_MINOR 9 )
set ( LAGraph_VERSION_SUB   13 )

project ( lagraph
          VERSION "${LAGraph_VERSION_MAJOR}.${LAGraph_VERSION_MINOR}.${LAGraph_VERSION_SUB}" )

#-------------------------------------------------------------------------------
# code coverage and build type
#-------------------------------------------------------------------------------

# To compile with test coverage:
#   cd build
#   cmake -DCOVERAGE=1 ..
#   make -j8
#   make test_coverage

if ( "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
	if ( COVERAGE )
                message ( STATUS "============== Code coverage enabled ===============" )
		set ( CMAKE_BUILD_TYPE Debug )
                # On the Mac, you need gcov-11 from homebrew (part of gcc-11):
                # and uncomment this line:
                # set ( GCOV_PATH /usr/local/bin/gcov-11)
		include ( CodeCoverage )

		append_coverage_compiler_flags ( )

		# turn off optimization for non-skewed coverage reports
		set ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -DCOVERAGE" )
		set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0" )

		setup_target_for_coverage_lcov (
			NAME test_coverage
			EXECUTABLE ctest
			DEPENDENCIES ${PROJECT_NAME}
			BASE_DIRECTORY "."
			NO_DEMANGLE TRUE
                        EXCLUDE "*/benchmark/*" "deps/json*/*"
                            "src/test/include/acutest.h"
		)
	endif ( )
endif ( )

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
message ( STATUS "GraphBLAS root:              " ${GRAPHBLAS_ROOT} )
find_package (GraphBLAS REQUIRED MODULE)

#-------------------------------------------------------------------------------
# Is this really necessary?
#-------------------------------------------------------------------------------

# Tim: yes, if we want a "sudo make install" to work.
# With also need to use the "install ( ... )" command.
# See the last part of the SuiteSparse:GraphBLAS CMakeLists.txt

include ( GNUInstallDirs )

#-------------------------------------------------------------------------------
# determine what user threading model to use
#-------------------------------------------------------------------------------

if ( COVERAGE )
    message ( STATUS "OpenMP disabled for test coverage" )
else ( )
    include ( FindOpenMP  )
    include ( FindThreads )
    if ( OPENMP_FOUND )
        set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} " )
    endif ( )
endif ( )

#-------------------------------------------------------------------------------
# report status
#-------------------------------------------------------------------------------

message ( STATUS "CMAKE build type:          " ${CMAKE_BUILD_TYPE} )
message ( STATUS "CMAKE source directory:    " ${CMAKE_SOURCE_DIR} )
message ( STATUS "CMAKE build directory:     " ${CMAKE_BINARY_DIR} )

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

# The GraphBLAS include directory is placed last.  If it is /usr/local/include
# and there is already a /usr/local/include/LAGraph.h present there, then that
# LAGraph.h may be an older version of LAGraph.

include_directories ( ${CMAKE_SOURCE_DIR}/include
                      ${CMAKE_SOURCE_DIR}/src/utility
                      ${CMAKE_SOURCE_DIR}/test/include
                      ${GRAPHBLAS_INCLUDE_DIR} )
include_directories("/usr/local/include")

set ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DLGDIR=${CMAKE_SOURCE_DIR}")

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
# force use of vanilla code
#-------------------------------------------------------------------------------

if ( LAGRAPH_VANILLA )
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -DLAGRAPH_VANILLA=1 " )
    message ( STATUS "Vanilla build: not relying on any GxB extensions" )
else ( )
    message ( STATUS "GxB build: relying on SuiteSparse GxB extensions" )
endif ( )

#-------------------------------------------------------------------------------
# print final C flags
#-------------------------------------------------------------------------------

message ( STATUS "CMAKE C flags: " ${CMAKE_C_FLAGS} )

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# CTEST cannot find the binaries when in this location
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

enable_testing()

add_subdirectory(src)
add_subdirectory(experimental)

#-------------------------------------------------------------------------------
# to install use "sudo make install"
#-------------------------------------------------------------------------------

# only install the dynamic  library
install ( TARGETS lagraph lagraphx
         LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
         PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
         ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR} )
