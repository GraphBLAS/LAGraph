//------------------------------------------------------------------------------
// LAGraph/src/test/test_NameOfType .c:  test LAGraph_NameOfType 
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
char name [LAGRAPH_MAX_NAME_LEN] ;
char msg [LAGRAPH_MSG_LEN] ;

typedef int myint ;

//------------------------------------------------------------------------------
// test_NameOfType :  test LAGraph_NameOfType 
//------------------------------------------------------------------------------

void test_NameOfType  (void)
{
    OK (LAGraph_Init (msg)) ;

    OK (LAGraph_NameOfType  (name, GrB_BOOL, msg)) ;
    OK (strcmp (name, "bool")) ;

    OK (LAGraph_NameOfType  (name, GrB_INT8, msg)) ;
    OK (strcmp (name, "int8_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_INT16, msg)) ;
    OK (strcmp (name, "int16_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_INT32, msg)) ;
    OK (strcmp (name, "int32_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_INT64, msg)) ;
    OK (strcmp (name, "int64_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_UINT8, msg)) ;
    OK (strcmp (name, "uint8_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_UINT16, msg)) ;
    OK (strcmp (name, "uint16_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_UINT32, msg)) ;
    OK (strcmp (name, "uint32_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_UINT64, msg)) ;
    OK (strcmp (name, "uint64_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_FP32, msg)) ;
    OK (strcmp (name, "float")) ;

    OK (LAGraph_NameOfType  (name, GrB_FP64, msg)) ;
    OK (strcmp (name, "double")) ;

    OK (GrB_Type_new (&type, sizeof (myint))) ;
    OK (LAGraph_NameOfType  (name, type, msg)) ;
    printf ("name: [%s]\n", name) ;
    OK (strcmp (name, "user_defined_type")) ;

    TEST_CHECK (LAGraph_NameOfType (NULL, NULL, msg) == GrB_NULL_POINTER) ;
    printf ("\nmsg: %s\n", msg) ;

    TEST_CHECK (LAGraph_NameOfType (name, NULL, msg) == GrB_NULL_POINTER) ;
    printf ("msg: %s\n", msg) ;

    TEST_CHECK (LAGraph_NameOfType (NULL, GrB_BOOL, msg) == GrB_NULL_POINTER) ;
    printf ("msg: %s\n", msg) ;

    OK (LAGraph_Finalize (msg)) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "TypeName", test_NameOfType  },
    // no brutal test needed
    { NULL, NULL }
} ;

