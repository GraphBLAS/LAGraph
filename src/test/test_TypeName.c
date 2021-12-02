//------------------------------------------------------------------------------
// LAGraph/src/test/test_TypeName.c:  test LAGraph_TypeName
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#include "LAGraph_test.h"

//------------------------------------------------------------------------------
// global variables
//------------------------------------------------------------------------------

GrB_Type type = NULL ;
char *name = NULL ;
char msg [LAGRAPH_MSG_LEN] ;

typedef int myint ;

//------------------------------------------------------------------------------
// test_TypeName:  test LAGraph_TypeName
//------------------------------------------------------------------------------

void test_TypeName (void)
{
    OK (LAGraph_Init (msg)) ;

    OK (LAGraph_TypeName (&name, GrB_BOOL, msg)) ;
    OK (strcmp (name, "bool")) ;

    OK (LAGraph_TypeName (&name, GrB_INT8, msg)) ;
    OK (strcmp (name, "int8")) ;

    OK (LAGraph_TypeName (&name, GrB_INT16, msg)) ;
    OK (strcmp (name, "int16")) ;

    OK (LAGraph_TypeName (&name, GrB_INT32, msg)) ;
    OK (strcmp (name, "int32")) ;

    OK (LAGraph_TypeName (&name, GrB_INT64, msg)) ;
    OK (strcmp (name, "int64")) ;

    OK (LAGraph_TypeName (&name, GrB_UINT8, msg)) ;
    OK (strcmp (name, "uint8")) ;

    OK (LAGraph_TypeName (&name, GrB_UINT16, msg)) ;
    OK (strcmp (name, "uint16")) ;

    OK (LAGraph_TypeName (&name, GrB_UINT32, msg)) ;
    OK (strcmp (name, "uint32")) ;

    OK (LAGraph_TypeName (&name, GrB_UINT64, msg)) ;
    OK (strcmp (name, "uint64")) ;

    OK (LAGraph_TypeName (&name, GrB_FP32, msg)) ;
    OK (strcmp (name, "single")) ;

    OK (LAGraph_TypeName (&name, GrB_FP64, msg)) ;
    OK (strcmp (name, "double")) ;

    OK (GrB_Type_new (&type, sizeof (myint))) ;
    OK (LAGraph_TypeName (&name, type, msg)) ;
    OK (strcmp (name, "user-defined")) ;

    TEST_CHECK (LAGraph_TypeName (NULL, NULL, msg) == GrB_NULL_POINTER) ;
    printf ("\nmsg: %s\n", msg) ;

    TEST_CHECK (LAGraph_TypeName (&name, NULL, msg) == GrB_NULL_POINTER) ;
    printf ("msg: %s\n", msg) ;

    TEST_CHECK (LAGraph_TypeName (NULL, GrB_BOOL, msg) == GrB_NULL_POINTER) ;
    printf ("msg: %s\n", msg) ;

    OK (LAGraph_Finalize (msg)) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "TypeName", test_TypeName },
    // no brutal test needed
    { NULL, NULL }
} ;

