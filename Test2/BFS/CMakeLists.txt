#-------------------------------------------------------------------------------
# LAGraph/Test2/BFS/CMakeLists.txt:  cmake script for LAGraph test
#-------------------------------------------------------------------------------

# cmake 3.0 is preferred.
cmake_minimum_required ( VERSION 2.8.12 )
if ( CMAKE_VERSION VERSION_GREATER "3.0" )
    cmake_policy ( SET CMP0042 NEW )
    cmake_policy ( SET CMP0048 NEW )
endif ( )

project ( test_bfs )

include ( GNUInstallDirs )

# For development only, not for end-users:
# set ( CMAKE_BUILD_TYPE Debug )
if ( NOT CMAKE_BUILD_TYPE )
    set ( CMAKE_BUILD_TYPE Release )
endif ( )

# link_directories ( /usr/local/lib )
link_directories ( ../../../GraphBLAS/build )
link_directories ( ../../build )

# include_directories ( /usr/local/include )
include_directories ( ../../../GraphBLAS/Include )
include_directories ( ../../Include )

# find OpenMP
include ( FindOpenMP  )

# report status
message ( STATUS "CMAKE build type:  " ${CMAKE_BUILD_TYPE} )
message ( STATUS "CMAKE compiler ID: " ${CMAKE_C_COMPILER_ID} )
message ( STATUS "CMAKE have OpenMP: " ${OPENMP_FOUND} )

if ( OPENMP_FOUND )
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
endif ( )

# include directories
set ( CMAKE_INCLUDE_CURRENT_DIR ON )

#-------------------------------------------------------------------------------
# compiler options
#-------------------------------------------------------------------------------

# check which compiler is being used.  If you need to make
# compiler-specific modifications, here is the place to do it.
if ( "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    # cmake 2.8 workaround: gcc needs to be told to do ANSI C11.
    # cmake 3.0 doesn't have this problem.
    # set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -g " )
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -std=c11 -lm -Wno-pragmas " )
    # check all warnings:
#   set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -Wextra -Wpedantic -Werror " )
    if ( CMAKE_C_COMPILER_VERSION VERSION_LESS 4.9 )
        message ( FATAL_ERROR "gcc version must be at least 4.9" )
    endif ( )
elseif ( "${CMAKE_C_COMPILER_ID}" STREQUAL "Intel" )
    # options for icc: also needs -std=c11
    # set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -g " )
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -std=c11 " )
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -qopt-malloc-options=3" )
    # check all warnings:
#   set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -w3 " )
    if ( CMAKE_C_COMPILER_VERSION VERSION_LESS 18.0 )
        message ( FATAL_ERROR "icc version must be at least 18.0" )
    endif ( )
elseif ( "${CMAKE_C_COMPILER_ID}" STREQUAL "Clang" )
    # options for clang
    if ( CMAKE_C_COMPILER_VERSION VERSION_LESS 3.3 )
        message ( FATAL_ERROR "clang version must be at least 3.3" )
    endif ( )
elseif ( "${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC" )
    # options for MicroSoft Visual Studio
endif ( )

if ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}" )
else ( )
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}" )
endif ( )

#-------------------------------------------------------------------------------
# select the threading library 
#-------------------------------------------------------------------------------

# print final C flags
message ( STATUS "CMAKE C flags: " ${CMAKE_C_FLAGS} )

#-------------------------------------------------------------------------------
# Demo programs
#-------------------------------------------------------------------------------

file ( GLOB SOURCES "*.c" )

add_executable ( test_bfs test_bfs.c )

# Libraries required for Demo programs
target_link_libraries ( test_bfs graphblas lagraph m )

