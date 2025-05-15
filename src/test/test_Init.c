//------------------------------------------------------------------------------
// LAGraph/src/test/test_Init.c:  test LAGraph_Init and LAGraph_Finalize
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LAGraph_test.h"
#include "LG_internal.h"

//------------------------------------------------------------------------------
// global variables
//------------------------------------------------------------------------------

char msg [LAGRAPH_MSG_LEN] ;

//------------------------------------------------------------------------------
// test_Init:  test LAGraph_Init
//------------------------------------------------------------------------------

void test_Init (void)
{

    int status = LAGraph_Init (msg) ;
    OK (status) ;
    int32_t ver [3] ;
    char library [256], date [256] ;

    OK (GrB_get (GrB_GLOBAL, library, GrB_NAME)) ;
    OK (GrB_get (GrB_GLOBAL, &(ver [0]), GrB_LIBRARY_VER_MAJOR)) ;
    OK (GrB_get (GrB_GLOBAL, &(ver [1]), GrB_LIBRARY_VER_MINOR)) ;
    OK (GrB_get (GrB_GLOBAL, &(ver [2]), GrB_LIBRARY_VER_PATCH)) ;
    date [0] = '\0' ;
    OK (LG_GET_LIBRARY_DATE (date)) ;

    printf ("\nlibrary: %s %d.%d.%d (%s)\n", library, ver [0], ver [1], ver [2],
        date) ;

    #if LAGRAPH_SUITESPARSE
    printf (  "include: %s %d.%d.%d (%s)\n", GxB_IMPLEMENTATION_NAME,
        GxB_IMPLEMENTATION_MAJOR, GxB_IMPLEMENTATION_MINOR,
        GxB_IMPLEMENTATION_SUB, GxB_IMPLEMENTATION_DATE) ;
    // make sure the SuiteSparse:GraphBLAS version and date match
    TEST_CHECK (ver [0] == GxB_IMPLEMENTATION_MAJOR) ;
    TEST_CHECK (ver [1] == GxB_IMPLEMENTATION_MINOR) ;
    TEST_CHECK (ver [2] == GxB_IMPLEMENTATION_SUB) ;
    OK (strcmp (date, GxB_IMPLEMENTATION_DATE)) ;
    char compiler [1024] ;
    int compiler_version [3] ;
    OK (GrB_Global_get_String (GrB_GLOBAL, compiler, GxB_COMPILER_NAME)) ;
    OK (GrB_Global_get_VOID   (GrB_GLOBAL, (void *) compiler_version,
        GxB_COMPILER_VERSION)) ;
    printf ("GraphBLAS compiled with: %s v%d.%d.%d\n", compiler,
        compiler_version [0], compiler_version [1], compiler_version [2]) ;
    #endif

    // check the LAGraph version using both LAGraph.h and LAGraph_Version
    printf ("LAGraph version %d.%d.%d (%s) from LAGraph.h\n",
        LAGRAPH_VERSION_MAJOR, LAGRAPH_VERSION_MINOR, LAGRAPH_VERSION_UPDATE,
        LAGRAPH_DATE) ;

    char version_date [LAGRAPH_MSG_LEN] ;
    status = LAGraph_Version (ver, version_date, msg) ;
    OK (status) ;

    printf ("LAGraph version %d.%d.%d (%s) from LAGraph_Version\n",
        ver [0], ver [1], ver [2], version_date) ;

    // make sure the LAGraph version and date match
    TEST_CHECK (ver [0] == LAGRAPH_VERSION_MAJOR) ;
    TEST_CHECK (ver [1] == LAGRAPH_VERSION_MINOR) ;
    TEST_CHECK (ver [2] == LAGRAPH_VERSION_UPDATE) ;
    OK (strcmp (version_date, LAGRAPH_DATE)) ;

    OK (LAGraph_Finalize (msg)) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "Init", test_Init },
    // no brutal test: see test_Xinit
    { NULL, NULL }
} ;
