//------------------------------------------------------------------------------
// LAGraph/src/test/test_vector.c:  test LG_check_vector
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

GrB_Vector X = NULL ;
GrB_Index *x = NULL ;
GrB_Index n = 10 ;
int64_t missing = 42 ;
char msg [LAGRAPH_MSG_LEN] ;

//------------------------------------------------------------------------------
// test_vector
//------------------------------------------------------------------------------

void test_vector (void)
{
    OK (LAGraph_Init (msg)) ;
    OK (GrB_Vector_new (&X, GrB_INT64, n)) ;
    OK (GrB_assign (X, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    OK (GrB_apply (X, NULL, NULL, GrB_ROWINDEX_INT64, X, 0, NULL)) ;
    OK (GrB_Vector_removeElement (X, 3)) ;
    OK (LAGraph_Vector_print (X, 3, stdout, msg)) ;
    x = LAGraph_Malloc (n, sizeof (int64_t)) ;
    OK (LG_check_vector (x, X, n, missing)) ;
    for (int i = 0 ; i < n ; i++)
    {
        TEST_CHECK (x [i] == ((i == 3) ? missing : i)) ;
    }
    OK (GrB_free (&X)) ;
    LAGraph_Free ((void **) &x) ;
    OK (LAGraph_Finalize (msg)) ;
}

//------------------------------------------------------------------------------
// test_vector_brutal
//------------------------------------------------------------------------------

#if LG_SUITESPARSE
void test_vector_brutal (void)
{
    OK (LG_brutal_setup (msg)) ;
    printf ("\n") ;

    for (int nbrutal = 0 ; ; nbrutal++)
    {
        /* allow for only nbrutal mallocs before 'failing' */
        LG_brutal = nbrutal ;
        /* try the method with brutal malloc */
        int brutal_result = GrB_Vector_new (&X, GrB_INT64, n) ;
        if (brutal_result != GrB_SUCCESS) continue ;
        brutal_result = GrB_assign (X, NULL, NULL, 0, GrB_ALL, n, NULL) ;
        if (brutal_result >= 0)
        {
            /* the method finally succeeded */
            break ;
        }
        GrB_free (&X) ;
        if (nbrutal > 10000) { printf ("Infinite!\n") ; abort ( ) ; }
    }
    LG_brutal = -1 ;  /* turn off brutal mallocs */

    OK (GrB_apply (X, NULL, NULL, GrB_ROWINDEX_INT64, X, 0, NULL)) ;
    OK (GrB_Vector_removeElement (X, 3)) ;
    OK (LAGraph_Vector_print (X, 3, stdout, msg)) ;

    x = LAGraph_Malloc (n, sizeof (int64_t)) ;
    LG_BRUTAL (LG_check_vector (x, X, n, missing)) ;

    for (int i = 0 ; i < n ; i++)
    {
        TEST_CHECK (x [i] == ((i == 3) ? missing : i)) ;
    }

    OK (GrB_free (&X)) ;
    LAGraph_Free ((void **) &x) ;
    OK (LG_brutal_teardown (msg)) ;
}
#endif

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "vector", test_vector },
    { "vector_brutal", test_vector_brutal },
    { NULL, NULL }
} ;

