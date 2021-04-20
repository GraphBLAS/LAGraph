#-------------------------------------------------------------------------------
# LAGraph/src/demo/CMakeLists.txt:  cmake script for LAGraph/src/demo
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
#
# See additional acknowledgments in the LICENSE file,
# or contact permission@sei.cmu.edu for the full terms.

#-------------------------------------------------------------------------------

file( GLOB DEMO_SOURCES LIST_DIRECTORIES false *.c )
foreach( testsourcefile ${DEMO_SOURCES} )
    get_filename_component(justname ${testsourcefile} NAME)
    string( REPLACE ".c" "" testname ${justname} )
    message("Adding: ${testname}")
    add_executable( ${testname} ${testsourcefile})
    target_link_libraries( ${testname} lagraph ${GRAPHBLAS_LIBRARIES} )
endforeach( testsourcefile ${DEMO_SOURCES} )
