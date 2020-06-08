//------------------------------------------------------------------------------
// LAGraph/Test/VectorToDense/vector_to_dense_test.c: test program for
// LAGraph_Vector_to_dense
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// Contributed by Marton Elekes, BME

// Usage:
//
// vector_to_dense_test

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&v) ;                 \
    GrB_free (&v_dense) ;           \
    GrB_free (&v_dense_ref) ;       \
}

#define ASSERT_TRUE(expr)                                     \
{                                                             \
    if(!(expr)) {                                             \
        fprintf(stderr, "\nTEST FAILED: %s\nFile: %s:%d\n\n", \
                #expr, __FILE__, __LINE__);                   \
        LAGRAPH_FREE_ALL;                                     \
        exit(EXIT_FAILURE);                                   \
    }                                                         \
}


int main(void) {

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Vector v = NULL;
    GrB_Vector v_dense = NULL;
    GrB_Vector v_dense_ref = NULL;

    LAGRAPH_TRY_CATCH (LAGraph_init());
    GxB_set (GxB_BURBLE, true) ;
    char *date, *compile_date, *compile_time ;
    int version [3] ;
    GxB_get (GxB_LIBRARY_VERSION, version) ;
    GxB_get (GxB_LIBRARY_DATE, &date) ;
    GxB_get (GxB_LIBRARY_COMPILE_DATE, &compile_date) ;
    GxB_get (GxB_LIBRARY_COMPILE_TIME, &compile_time) ;
    printf ("Library version %d.%d.%d\n", version [0], version [1], version [2]) ;
    printf ("Library date: %s\n", date) ;
    printf ("Compiled at %s on %s\n", compile_time, compile_date) ;

    GrB_Index I[] = {4, 0, 1};
    uint64_t X[] = {1, 2, 3};
    GrB_Index n = 6;
    GrB_Index nvals = sizeof(I) / sizeof(I[0]);
    GrB_Type type = GrB_UINT64;
    GrB_BinaryOp dup = GrB_PLUS_UINT64;

    // prepare reference arrays for dense vector
    GrB_Index I_ref[n];
    uint64_t X_ref[n];
    for (GrB_Index i = 0; i < n; ++i) {
        I_ref[i] = i;
        X_ref[i] = 0;
    }
    for (GrB_Index i = 0; i < nvals; ++i) {
        X_ref[I[i]] = X[i];
    }

    LAGr_Vector_new(&v, type, n);
    LAGr_Vector_new(&v_dense, type, n);
    LAGr_Vector_new(&v_dense_ref, type, n);

    LAGr_Vector_build(v, I, X, nvals, dup);
    LAGr_Vector_build(v_dense_ref, I_ref, X_ref, n, dup);

    uint64_t zero = 0;
    LAGRAPH_TRY_CATCH(LAGraph_Vector_to_dense(&v_dense, v, &zero));

    LAGRAPH_TRY_CATCH(GxB_fprint(v, GxB_COMPLETE, stdout));
    LAGRAPH_TRY_CATCH(GxB_fprint(v_dense, GxB_COMPLETE, stdout));
    LAGRAPH_TRY_CATCH(GxB_fprint(v_dense_ref, GxB_COMPLETE, stdout));

    // test
    bool isequal = false;
    LAGRAPH_TRY_CATCH(LAGraph_Vector_isequal(&isequal, v_dense, v_dense_ref, GrB_NULL));
    ASSERT_TRUE(isequal);

    LAGRAPH_FREE_ALL;
    LAGraph_finalize();
}
