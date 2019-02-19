//------------------------------------------------------------------------------
// LAGraph_BF_basic.c: Bellman-Ford method for single source shortest paths
// without returning the parents in the paths
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// LAGraph_BF performs a Bellman-Ford to find out shortest from given source
// vertex s. The result is vector d where d(k) is the shortest distance from s
// to k.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL   \
{                          \
	GrB_free(&dtmp) ;  \
}


/*
* Given a n-by-n adjacency matrix A and a source vertex s.
* If there is no negative-weight cycle reachable from s, return the distances of
* shortest paths from s as vector d. Otherwise, return d=NULL if there is
* negative-weight cycle.
* pd is given pointer to a GrB_Vector
* A has zeros on diagonal and weights on corresponding entries of edges
* s is given index for source vertex
*/
GrB_Info LAGraph_BF_basic
(
    GrB_Vector *pd,             //the pointer to the vector of distance
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Index s           //given index of the source
)
{
    GrB_Info info;
    GrB_Index n;           // n = # of vertices in graph
    // tmp vector to store distance vector after n (i.e., V) loops
    GrB_Vector dtmp = NULL;
    
    if (A == NULL || pd == NULL)
    {
        // required argument is missing
        LAGRAPH_ERROR ("required arguments are NULL", GrB_NULL_POINTER) ;
    }

    LAGRAPH_OK (GrB_Matrix_nrows(&n, A));
    
    if (s >= n || s < 0)
    {
        LAGRAPH_ERROR ("invalid value for source vertex s", GrB_INVALID_VALUE) ;
    }

    // Initialize distance vector, change the d[s] to 0
    LAGRAPH_OK (GrB_Vector_new(pd, GrB_FP64, n));
    LAGRAPH_OK (GrB_Vector_setElement_FP64(*pd, 0, s));

    // copy d to dtmp in order to create a same size of vector
    LAGRAPH_OK (GrB_Vector_dup(&dtmp, *pd));
   
    int32_t count = 0;     //number of iterations
    bool same = false;     //variable indicating if d=dtmp

    // terminate when no new path is found or more than n-1 loops
    while (!same && count < n - 1)
    {
        // excute semiring on d and A, and save the result to d
        LAGRAPH_OK (GrB_mxv(dtmp, GrB_NULL, GrB_NULL, GxB_MIN_PLUS_FP64, A, 
            *pd, GrB_NULL));    
        LAGRAPH_OK (LAGraph_Vector_isequal(&same, dtmp, *pd, GrB_NULL));
	if (!same)
        {
            GrB_Vector ttmp = dtmp;
            dtmp = (*pd);
            (*pd) = ttmp;
        }
        count++;
    }

    // check for negative-weight cycle only when there was a new path in the 
    // last loop, otherwise, there can't be a negative-weight cycle.
    if (!same)
    {
        // excute semiring again to check for negative-weight cycle
        LAGRAPH_OK (GrB_mxv(dtmp, GrB_NULL, GrB_NULL, GxB_MIN_PLUS_FP64, A, 
            *pd, GrB_NULL));
        LAGRAPH_OK (LAGraph_Vector_isequal(&same, dtmp, *pd, GrB_NULL));

	// if d != dtmp, then there is a negative-weight cycle in the graph
        if (!same)
        {
            printf("A negative-weight cycle found. \n");
	    GrB_free (pd);
            LAGRAPH_FREE_ALL;
            return (GrB_SUCCESS) ;
        }
    }

    LAGRAPH_FREE_ALL;
    return (GrB_SUCCESS) ;
}
