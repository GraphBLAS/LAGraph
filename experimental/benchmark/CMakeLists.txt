#-------------------------------------------------------------------------------
# LAGraph/experimental/benchmark/CMakeLists.txt:  cmake script for LAGraph/experiment/benchmark
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
# See additional acknowledgments in the LICENSE file,
# or contact permission@sei.cmu.edu for the full terms.

#-------------------------------------------------------------------------------

include_directories ( ${CMAKE_SOURCE_DIR}/src/test/include
    ${CMAKE_SOURCE_DIR}/src/benchmark
    ${CMAKE_SOURCE_DIR}/experimental/test/include )

file( GLOB DEMO_SOURCES LIST_DIRECTORIES false *_demo.c )
foreach( demosourcefile ${DEMO_SOURCES} )
    get_filename_component(justname ${demosourcefile} NAME)
    string( REPLACE ".c" "" demoname ${justname} )
#   message("Adding: ${demoname}")
    add_executable( ${demoname} ${demosourcefile} )
    target_link_libraries( ${demoname} lagraphx lagraph lagraphtest lagraphxtest ${GRAPHBLAS_LIBRARIES} )
    target_link_directories( ${demoname} BEFORE PUBLIC ${CMAKE_SOURCE_DIR}/build )
endforeach( demosourcefile ${DEMO_SOURCES} )
