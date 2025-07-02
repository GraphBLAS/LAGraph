cmake_minimum_required ( VERSION 3.19 )

# policies
cmake_policy ( SET CMP0042 NEW )    # enable MACOSX_RPATH by default
cmake_policy ( SET CMP0048 NEW )    # VERSION variable policy
cmake_policy ( SET CMP0054 NEW )    # if ( expression ) handling policy
cmake_policy ( SET CMP0104 NEW )    # initialize CUDA architectures
cmake_policy ( SET CMP0144 NEW )    # uppercase find_package names

set ( CMAKE_MACOSX_RPATH TRUE )
enable_language ( C CXX )
include ( GNUInstallDirs )

set ( LAGRAPH_ROOT "../../..")
set ( LAGRAPHX_ROOT "../../..")
set ( GRAPHBLAS_ROOT "../../../../GraphBLAS")
set ( GraphBLAS_CUDA_DIR "../../../../GraphBLAS/build/CUDA")

# set the module path for all Find*.cmake files.
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake_modules
    ${CMAKE_SOURCE_DIR}/../../../cmake_modules
	${CMAKE_SOURCE_DIR}/../../../../GraphBLAS/cmake_modules ) 

include ( SuiteSparsePolicy )

project ( mytests LANGUAGES C CXX )

add_executable ( gen_bipartite "gen_bipartite.cpp" )
add_executable ( gen_general "gen_general.cpp")
add_executable ( verify_matching "verify_matching.cpp")

#-------------------------------------------------------------------------------
# find library dependencies
#-------------------------------------------------------------------------------

# look for all SuiteSparse packages:
find_package ( LAGraph 
    PATHS ${LAGRAPH_ROOT}/build NO_DEFAULT_PATH REQUIRED )

find_library ( LAGRAPHX lagraphx PATHS ../../../build/experimental REQUIRED ) 
set ( LAGRAPHX_INCLUDE_DIR ${LAGRAPHX_ROOT}/include )

find_package ( GraphBLAS REQUIRED
    PATHS ${GRAPHBLAS_ROOT}/build NO_DEFAULT_PATH REQUIRED )

# look for all external libaries:
find_package ( OpenMP REQUIRED )

#-------------------------------------------------------------------------------
# add the library dependencies
#-------------------------------------------------------------------------------

# OpenMP:
message ( STATUS "OpenMP C libraries:      ${OpenMP_C_LIBRARIES} ")
message ( STATUS "OpenMP C include:        ${OpenMP_C_INCLUDE_DIRS} ")
message ( STATUS "OpenMP C flags:          ${OpenMP_C_FLAGS} ")
target_link_libraries ( gen_bipartite PUBLIC ${OpenMP_C_LIBRARIES} )
target_link_libraries ( gen_general PUBLIC ${OpenMP_C_LIBRARIES} )
set ( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS} " )
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} " )
include_directories ( ${OpenMP_C_INCLUDE_DIRS} )

# libm:
if ( NOT WIN32 )
    target_link_libraries ( gen_bipartite PUBLIC m )
    target_link_libraries ( gen_general PUBLIC m )
    target_link_libraries ( verify_matching PUBLIC m )
endif ( )

# GraphBLAS:
target_link_libraries ( gen_bipartite PUBLIC ${GRAPHBLAS_LIBRARIES} )
target_link_libraries ( gen_general PUBLIC ${GRAPHBLAS_LIBRARIES} )
include_directories ( ${GRAPHBLAS_INCLUDE_DIR} )

# LAGraph:
target_link_libraries ( gen_bipartite PUBLIC ${LAGRAPH_LIBRARIES} )
target_link_libraries ( gen_general PUBLIC ${LAGRAPH_LIBRARIES} )
include_directories ( ${LAGRAPH_INCLUDE_DIR} )

# LAGraphX:
target_link_libraries ( gen_bipartite PUBLIC ${LAGRAPHX} )
target_link_libraries ( gen_general PUBLIC ${LAGRAPHX} )
include_directories ( ${LAGRAPHX_INCLUDE_DIR} )
