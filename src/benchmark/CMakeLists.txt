#-------------------------------------------------------------------------------
# LAGraph/src/demo/CMakeLists.txt:  cmake script for LAGraph/src/demo
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
# See additional acknowledgments in the LICENSE file,
# or contact permission@sei.cmu.edu for the full terms.

#-------------------------------------------------------------------------------

include_directories ( ${CMAKE_SOURCE_DIR}/src/test/include
    ${CMAKE_SOURCE_DIR}/src/algorithm )

set(DEMO_SOURCES
    bc_demo.c
    bfs_demo.c
    cc_demo.c
    gappagerank_demo.c
    mtx2bin_demo.c
    ss2_demo.c
    sssp_demo.c
    tc_demo.c)

if(ENABLE_NVIDIA)
    list(APPEND DEMO_SOURCES tc_gpu_demo.c)
endif()

foreach( demosourcefile ${DEMO_SOURCES} )
    get_filename_component(justname ${demosourcefile} NAME)
    string( REPLACE ".c" "" demoname ${justname} )
#   message("Adding: ${demoname}")
    add_executable( ${demoname} ${demosourcefile} )
    target_link_libraries( ${demoname} lagraph lagraphx lagraphtest ${GRAPHBLAS_LIBRARIES} )
    target_link_directories( ${demoname} BEFORE PUBLIC ${CMAKE_SOURCE_DIR}/build )
endforeach( demosourcefile ${DEMO_SOURCES} )
