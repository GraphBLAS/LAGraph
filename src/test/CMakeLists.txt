#-------------------------------------------------------------------------------
# LAGraph/src/test/CMakeLists.txt:  cmake script for LAGraph/src/test
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
#
# See additional acknowledgments in the LICENSE file,
# or contact permission@sei.cmu.edu for the full terms.

#-------------------------------------------------------------------------------
# This will only build tests from files with the name "test_*.c"
#-------------------------------------------------------------------------------

include_directories ( ${CMAKE_SOURCE_DIR}/src/test/include )

file( GLOB TEST_SOURCES LIST_DIRECTORIES false test_*.c )
foreach( testsourcefile ${TEST_SOURCES} )
    get_filename_component(justname ${testsourcefile} NAME)
    string( REPLACE ".c" "" testname ${justname} )
    message("Adding: ${testname}")
    add_executable( ${testname} ${testsourcefile})
    target_link_libraries( ${testname} lagraph ${GRAPHBLAS_LIBRARIES} )
    string( REPLACE "test_" "ctest_" ctestname ${testname})
    add_test( ${ctestname} ${testname})
endforeach( testsourcefile ${TEST_SOURCES} )
