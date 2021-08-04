#-------------------------------------------------------------------------------
# LAGraph/experimental/benchmark/CMakeLists.txt:  cmake script for LAGraph/experiment/benchmark
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
#
# See additional acknowledgments in the LICENSE file,
# or contact permission@sei.cmu.edu for the full terms.

#-------------------------------------------------------------------------------

file( GLOB DEMO_SOURCES LIST_DIRECTORIES false *_demo.c )
foreach( demosourcefile ${DEMO_SOURCES} )
    get_filename_component(justname ${demosourcefile} NAME)
    string( REPLACE ".c" "" demoname ${justname} )
    message("Adding: ${demoname}")
    add_executable( ${demoname} ${demosourcefile} )
    target_link_libraries( ${demoname} lagraphx lagraph ${GRAPHBLAS_LIBRARIES} )
endforeach( demosourcefile ${DEMO_SOURCES} )
