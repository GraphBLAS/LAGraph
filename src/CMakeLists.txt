#-------------------------------------------------------------------------------
# LAGraph/src/CMakeLists.txt:  cmake script for LAGraph/src
#-------------------------------------------------------------------------------

# LAGraph, (c) 2019-2023 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
#
# For additional details (including references to third party source code and
# other files) see the LICENSE file or contact permission@sei.cmu.edu. See
# Contributors.txt for a full list of contributors. Created, in part, with
# funding and support from the U.S. Government (see Acknowledgments.txt file).
# DM22-0790

#-------------------------------------------------------------------------------

include_directories ( utility )

file ( GLOB LAGRAPH_LIB_SOURCES "utility/*.c" "algorithm/*.c" )

#-------------------------------------------------------------------------------
# JIT kernel header generation (one _jit.h per source file)
#-------------------------------------------------------------------------------

set ( LAGRAPH_JIT_HEADERS "" )
foreach ( SRC IN LISTS LAGRAPH_LIB_SOURCES )
    get_filename_component ( SRC_NAME "${SRC}" NAME_WE )
    set ( JIT_HEADER "${JIT_INCLUDE_DIR}/${SRC_NAME}_jit.h" )
    add_custom_command (
        OUTPUT  "${JIT_HEADER}"
        COMMAND Python3::Interpreter
                "${PROJECT_SOURCE_DIR}/jit_generator.py"
                --input "${SRC}" --output "${JIT_HEADER}"
        DEPENDS "${SRC}" "${PROJECT_SOURCE_DIR}/jit_generator.py"
        COMMENT "JIT: generating ${SRC_NAME}_jit.h"
    )
    list ( APPEND LAGRAPH_JIT_HEADERS "${JIT_HEADER}" )
endforeach ()

add_custom_target ( LAGraphJITGenerate DEPENDS ${LAGRAPH_JIT_HEADERS} )

#-------------------------------------------------------------------------------
# dynamic lagraph library properties
#-------------------------------------------------------------------------------

if ( BUILD_SHARED_LIBS )
    add_library ( LAGraph SHARED ${LAGRAPH_LIB_SOURCES} )
    set_target_properties ( LAGraph PROPERTIES
        VERSION ${LAGraph_VERSION_MAJOR}.${LAGraph_VERSION_MINOR}.${LAGraph_VERSION_SUB}
        SOVERSION ${LAGraph_VERSION_MAJOR}
        OUTPUT_NAME lagraph
        C_STANDARD_REQUIRED ON
        C_STANDARD 11
        PUBLIC_HEADER "include/LAGraph.h"
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/dlls )

    if ( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.25" )
        set_target_properties ( LAGraph PROPERTIES EXPORT_NO_SYSTEM ON )
    endif ( )

    target_link_libraries ( LAGraph PRIVATE GraphBLAS::GraphBLAS )
    if ( NOT NO_LIBM )
        target_link_libraries ( LAGraph PRIVATE "m" )
    endif ( )

    target_include_directories ( LAGraph PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${SUITESPARSE_INCLUDEDIR}> )

    target_compile_definitions ( LAGraph PRIVATE LG_LIBRARY )
    target_compile_definitions ( LAGraph PUBLIC LG_DLL )
    # install ( FILES ${LAGRAPH_INCLUDE} DESTINATION include)
endif ( )

#-------------------------------------------------------------------------------
# static lagraph library properties
#-------------------------------------------------------------------------------

if ( BUILD_STATIC_LIBS )
    add_library ( LAGraph_static STATIC ${LAGRAPH_LIB_SOURCES} )
    SET_TARGET_PROPERTIES ( LAGraph_static PROPERTIES
        OUTPUT_NAME lagraph
        POSITION_INDEPENDENT_CODE OFF
        C_STANDARD_REQUIRED ON
        C_STANDARD 11
        PUBLIC_HEADER "include/LAGraph.h" )

    if ( MSVC OR ("${CMAKE_C_SIMULATE_ID}" STREQUAL "MSVC") )
        set_target_properties ( LAGraph_static PROPERTIES
            OUTPUT_NAME lagraph_static )
    endif ( )

    if ( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.25" )
        set_target_properties ( LAGraph_static PROPERTIES EXPORT_NO_SYSTEM ON )
    endif ( )

    if ( TARGET GraphBLAS::GraphBLAS_static )
        target_link_libraries ( LAGraph_static PUBLIC GraphBLAS::GraphBLAS_static )
    else ( )
        target_link_libraries ( LAGraph_static PUBLIC GraphBLAS::GraphBLAS )
    endif ( )
    if ( NOT NO_LIBM )
        target_link_libraries ( LAGraph_static PUBLIC "m" )
    endif ( )

    target_include_directories ( LAGraph_static PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${SUITESPARSE_INCLUDEDIR}> )
endif ( )

#-------------------------------------------------------------------------------
# add OpenMP
#-------------------------------------------------------------------------------

if ( LAGRAPH_HAS_OPENMP )
    if ( BUILD_SHARED_LIBS )
        target_link_libraries ( LAGraph PRIVATE OpenMP::OpenMP_C )
    endif ( )
    if ( BUILD_STATIC_LIBS )
        target_link_libraries ( LAGraph_static PRIVATE OpenMP::OpenMP_C )
    endif ( )
endif ( )

#-------------------------------------------------------------------------------
# JIT header dependency: ensure headers are generated before compilation
#-------------------------------------------------------------------------------

if ( BUILD_SHARED_LIBS )
    add_dependencies ( LAGraph LAGraphJITGenerate )
endif ( )
if ( BUILD_STATIC_LIBS )
    add_dependencies ( LAGraph_static LAGraphJITGenerate )
endif ( )

#-------------------------------------------------------------------------------
# applications
#-------------------------------------------------------------------------------

if ( BUILD_TESTING )
    add_subdirectory ( test )
    add_subdirectory ( benchmark )
endif ( )
