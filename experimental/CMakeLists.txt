#-------------------------------------------------------------------------------
# LAGraph/experimental/CMakeLists.txt:  cmake script for LAGraph/experimental
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

include_directories (
    ${PROJECT_SOURCE_DIR}/src/utility
    ${PROJECT_SOURCE_DIR}/deps/json_h )

file ( GLOB LAGRAPHX_LIB_SOURCES "utility/*.c" "algorithm/*.c" )

#-------------------------------------------------------------------------------
# JIT kernel header generation (one _jit.h per source file)
#-------------------------------------------------------------------------------

set ( LAGRAPHX_JIT_HEADERS "" )
foreach ( SRC IN LISTS LAGRAPHX_LIB_SOURCES )
    get_filename_component ( SRC_NAME "${SRC}" NAME_WE )
    set ( JIT_HEADER "${JIT_INCLUDE_DIR}/${SRC_NAME}_jit.h" )
    add_custom_command (
        OUTPUT  "${JIT_HEADER}"
        COMMAND Python3::Interpreter
                "${PROJECT_SOURCE_DIR}/tools/jit_generator.py"
                --input "${SRC}" --output "${JIT_HEADER}"
        DEPENDS "${SRC}" "${PROJECT_SOURCE_DIR}/tools/jit_generator.py"
        COMMENT "JIT: generating ${SRC_NAME}_jit.h"
    )
    list ( APPEND LAGRAPHX_JIT_HEADERS "${JIT_HEADER}" )
endforeach ()

add_custom_target ( LAGraphXJITGenerate DEPENDS ${LAGRAPHX_JIT_HEADERS} )

#-------------------------------------------------------------------------------
# dynamic lagraphx library properties
#-------------------------------------------------------------------------------

if ( BUILD_SHARED_LIBS )
    add_library ( LAGraphX SHARED ${LAGRAPHX_LIB_SOURCES} )
    set_target_properties ( LAGraphX PROPERTIES
        VERSION ${LAGraph_VERSION_MAJOR}.${LAGraph_VERSION_MINOR}.${LAGraph_VERSION_SUB}
        SOVERSION ${LAGraph_VERSION_MAJOR}
        OUTPUT_NAME lagraphx
        C_STANDARD_REQUIRED ON
        C_STANDARD 11
        PUBLIC_HEADER "include/LAGraphX.h"
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/dlls )

    if ( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.25" )
        set_target_properties ( LAGraphX PROPERTIES EXPORT_NO_SYSTEM ON )
    endif ( )

    target_link_libraries ( LAGraphX PRIVATE LAGraph GraphBLAS::GraphBLAS )

    target_include_directories ( LAGraphX PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${SUITESPARSE_INCLUDEDIR}> )

    target_compile_definitions ( LAGraphX PRIVATE LGX_LIBRARY )
    target_compile_definitions ( LAGraphX PUBLIC LGX_DLL )
endif ( )

#-------------------------------------------------------------------------------
# static lagraphx library properties
#-------------------------------------------------------------------------------

if ( BUILD_STATIC_LIBS )
    add_library ( LAGraphX_static STATIC ${LAGRAPHX_LIB_SOURCES} )
    set_target_properties ( LAGraphX_static PROPERTIES
        OUTPUT_NAME lagraphx
        POSITION_INDEPENDENT_CODE OFF
        C_STANDARD_REQUIRED ON
        C_STANDARD 11
        PUBLIC_HEADER "include/LAGraphX.h" )

    if ( MSVC OR ("${CMAKE_C_SIMULATE_ID}" STREQUAL "MSVC") )
        set_target_properties ( LAGraphX_static PROPERTIES
            OUTPUT_NAME lagraphx_static )
    endif ( )

    if ( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.25" )
        set_target_properties ( LAGraphX_static PROPERTIES EXPORT_NO_SYSTEM ON )
    endif ( )

    if ( TARGET GraphBLAS::GraphBLAS_static )
        target_link_libraries ( LAGraphX_static PRIVATE GraphBLAS::GraphBLAS_static )
    else ( )
        target_link_libraries ( LAGraphX_static PRIVATE GraphBLAS::GraphBLAS )
    endif ( )

    target_include_directories ( LAGraphX_static PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${SUITESPARSE_INCLUDEDIR}> )
endif ( )

#-------------------------------------------------------------------------------
# add OpenMP
#-------------------------------------------------------------------------------

if ( LAGRAPH_HAS_OPENMP )
    if ( BUILD_SHARED_LIBS )
        target_link_libraries ( LAGraphX PRIVATE OpenMP::OpenMP_C )
    endif ( )
    if ( BUILD_STATIC_LIBS )
        target_link_libraries ( LAGraphX_static PRIVATE OpenMP::OpenMP_C )
    endif ( )
endif ( )

#-------------------------------------------------------------------------------
# libm
#-------------------------------------------------------------------------------

if ( NOT NO_LIBM )
    if ( BUILD_SHARED_LIBS )
        target_link_libraries ( LAGraphX PRIVATE m )
    endif ( )
    if ( BUILD_STATIC_LIBS )
        target_link_libraries ( LAGraphX_static PUBLIC m )
    endif ( )
endif ( )

#-------------------------------------------------------------------------------
# JIT header dependency: ensure headers are generated before compilation
#-------------------------------------------------------------------------------

if ( BUILD_SHARED_LIBS )
    add_dependencies ( LAGraphX LAGraphXJITGenerate )
endif ( )
if ( BUILD_STATIC_LIBS )
    add_dependencies ( LAGraphX_static LAGraphXJITGenerate )
endif ( )

#-------------------------------------------------------------------------------
# applications
#-------------------------------------------------------------------------------

if ( BUILD_TESTING )
    add_subdirectory ( test )
    add_subdirectory ( benchmark )
endif ( )
