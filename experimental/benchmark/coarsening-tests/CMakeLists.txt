cmake_minimum_required ( VERSION 3.19 )

# policies
# cmake_policy ( SET CMP0042 NEW )    # enable MACOSX_RPATH by default
cmake_policy ( SET CMP0048 NEW )    # VERSION variable policy
cmake_policy ( SET CMP0054 NEW )    # if ( expression ) handling policy
cmake_policy ( SET CMP0104 NEW )    # initialize CUDA architectures

set ( CMAKE_MACOSX_RPATH TRUE )
enable_language ( C CXX )
include ( GNUInstallDirs )

message (STATUS  "source: ${CMAKE_SOURCE_DIR}" )
message (STATUS  "path: ${CMAKE_MODULE_PATH}" )
set ( LAGRAPH_ROOT "../../..")
set ( LAGRAPHX_ROOT "../../..")
set ( GRAPHBLAS_ROOT "../../../../GraphBLAS")

# set the module path for all Find*.cmake files.
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake_modules
    ${CMAKE_SOURCE_DIR}/../cmake_modules
    ${CMAKE_SOURCE_DIR}/../../../cmake_modules
    ${CMAKE_SOURCE_DIR}/../LAGraph/cmake_modules )

project ( mytests LANGUAGES C CXX )

add_executable ( coarsen "coarsen.cpp" )

#-------------------------------------------------------------------------------
# find library dependencies
#-------------------------------------------------------------------------------

# look for all SuiteSparse packages:
find_package ( LAGraph REQUIRED )
find_package ( LAGraphX REQUIRED )
find_package ( GraphBLAS REQUIRED )

# look for all external libaries:
find_package ( OpenMP REQUIRED )

#-------------------------------------------------------------------------------
# add the library dependencies
#-------------------------------------------------------------------------------

# OpenMP:
message ( STATUS "OpenMP C libraries:      ${OpenMP_C_LIBRARIES} ")
message ( STATUS "OpenMP C include:        ${OpenMP_C_INCLUDE_DIRS} ")
message ( STATUS "OpenMP C flags:          ${OpenMP_C_FLAGS} ")
target_link_libraries ( coarsen PUBLIC ${OpenMP_C_LIBRARIES} )
set ( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS} " )
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} " )
include_directories ( ${OpenMP_C_INCLUDE_DIRS} )

# libm:
if ( NOT WIN32 )
    target_link_libraries ( coarsen PUBLIC m )
endif ( )

# GraphBLAS:
target_link_libraries ( coarsen PUBLIC ${GRAPHBLAS_LIBRARIES} )
include_directories ( ${GRAPHBLAS_INCLUDE_DIR} )

# LAGraph:
target_link_libraries ( coarsen PUBLIC ${LAGRAPH_LIBRARIES} )
include_directories ( ${LAGRAPH_INCLUDE_DIR} )

# LAGraphX:
target_link_libraries ( coarsen PUBLIC ${LAGRAPHX_LIBRARIES} )
include_directories ( ${LAGRAPHX_INCLUDE_DIR} )

