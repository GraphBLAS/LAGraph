#-------------------------------------------------------------------------------
# LAGraph/experimental/CMakeLists.txt:  cmake script for LAGraph/experimental
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
#
# See additional acknowledgments in the LICENSE file,
# or contact permission@sei.cmu.edu for the full terms.

#-------------------------------------------------------------------------------

include_directories ( ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src/utility utility )

file ( GLOB LAGRAPHX_LIB_SOURCES "utility/*.c" "algorithm/*.c" )

#-------------------------------------------------------------------------------
# dynamic lagraphx library properties
#-------------------------------------------------------------------------------

add_library ( lagraphx SHARED ${LAGRAPHX_LIB_SOURCES} )
SET_TARGET_PROPERTIES ( lagraphx PROPERTIES
    VERSION ${LAGraph_VERSION_MAJOR}.${LAGraph_VERSION_MINOR}.${LAGraph_VERSION_SUB}
    SOVERSION ${LAGraph_VERSION_MAJOR}
    C_STANDARD_REQUIRED 11
    PUBLIC_HEADER "include/LAGraphX.h" )
set_property ( TARGET lagraphx PROPERTY C_STANDARD 11 )
target_link_libraries(lagraphx PUBLIC lagraph ${GRAPHBLAS_LIBRARIES})

#-------------------------------------------------------------------------------
# static lagraphx library properties
#-------------------------------------------------------------------------------

add_library ( lagraphx_static STATIC ${LAGRAPHX_LIB_SOURCES} )
SET_TARGET_PROPERTIES ( lagraphx_static PROPERTIES
    VERSION ${LAGraph_VERSION_MAJOR}.${LAGraph_VERSION_MINOR}.${LAGraph_VERSION_SUB}
    OUTPUT_NAME lagraphx
    POSITION_INDEPENDENT_CODE OFF
    SOVERSION ${LAGraph_VERSION_MAJOR}
    C_STANDARD_REQUIRED 11
    PUBLIC_HEADER "include/LAGraphX.h" )
set_property ( TARGET lagraphx_static PROPERTY C_STANDARD 11 )
target_link_libraries(lagraphx_static PUBLIC lagraph_static ${GRAPHBLAS_LIBRARIES})

#-------------------------------------------------------------------------------
# add OpenMP
#-------------------------------------------------------------------------------

if ( OPENMP_FOUND )
    target_link_libraries ( lagraphx PUBLIC lagraph ${OpenMP_C_LIBRARIES} )
    target_link_libraries ( lagraphx_static PUBLIC lagraph_static ${OpenMP_C_LIBRARIES} )
endif ( )

#-------------------------------------------------------------------------------
# applications
#-------------------------------------------------------------------------------

add_subdirectory(test)
add_subdirectory(benchmark)
#add_subdirectory(utility/test)
#add_subdirectory(algorithm/test)

#-------------------------------------------------------------------------------
# lagraph installation location
#-------------------------------------------------------------------------------

# install ( TARGETS lagraph lagraph_static
#     LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
#     PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
#     ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR} )
