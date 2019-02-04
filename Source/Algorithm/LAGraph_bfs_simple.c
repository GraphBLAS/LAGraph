//------------------------------------------------------------------------------
// LAGraph_bfs_simple:  simple BFS
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// By default in SuiteSparse:GraphBLAS, the matrix A is stored by row (CSR
// format).  The method below computes q*A via GrB_vxm, which results in a
// "push" BFS algorithm.  If instead A is stored by column (using an extension
// to SuiteSparse:GraphBLAS), the algorithm becomes a "pull" algorithm.

// A faster algorithm would use the push-pull strategy, which requires both
// A and AT = A' to be passed in.  If A is in CSR format, then vxm(q,A)
// does the "push" and mxv(AT,q) does the "pull".

// Reference: Carl Yang, Aydin Buluc, and John D. Owens. 2018. Implementing
// Push-Pull Efficiently in GraphBLAS. In Proceedings of the 47th International
// Conference on Parallel Processing (ICPP 2018). ACM, New York, NY, USA,
// Article 89, 11 pages. DOI: https://doi.org/10.1145/3225058.3225122

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL    \
{                           \
    GrB_free (&v) ;         \
    GrB_free (&q) ;         \
    GrB_free (&desc) ;      \
}

// Given an n-by-n adjacency matrix A and a source node s, performs a BFS
// traversal of the graph and sets v[i] to the level in which node i is
// visited (v[s] == 1).  If i is not reacheable from s, then v[i] = 0. (Vector
// v should be empty on input.)  The graph A need not be Boolean on input;
// if it isn't Boolean, the semiring will properly typecast it to Boolean.
// However, best performance is obtained if A has type GrB_BOOL.

// The matrix A can have explicit entries equal to zero.  These are safely
// ignored.

// TODO: allow for a set of source nodes (pass in q, not s)

GrB_Info LAGraph_bfs_simple
(
    GrB_Vector *v_output,   // v [i] is the BFS level of node i in the graph
    const GrB_Matrix A,     // input graph, treated as if boolean in semiring
    GrB_Index s,            // starting node of the BFS
    int64_t max_level       // max # of levels to search (<0: nothing,
                            // 1: just the source, 2: source and neighbors, etc)
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Info info ;

    GrB_Vector q = NULL ;           // nodes visited at each level
    GrB_Vector v = NULL ;           // result vector
    GrB_Descriptor desc = NULL ;    // descriptor for vxm

    if (v_output == NULL)
    {
        // required output argument is missing
        return (GrB_NULL_POINTER) ;
    }

    (*v_output) = NULL ;

    GrB_Index nrows, ncols ;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A)) ;
    if (nrows != ncols)
    {
        // A must be square
        return (GrB_INVALID_VALUE) ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Index n = nrows ;
    max_level = LAGRAPH_MIN (n, max_level) ;

    // create an empty vector v.  Assume int32 is sufficient
    LAGRAPH_OK (GrB_Vector_new (&v, GrB_INT32, n)) ;

    // create a boolean vector q, and set q[s] to true
    LAGRAPH_OK (GrB_Vector_new (&q, GrB_BOOL, n)) ;
    LAGRAPH_OK (GrB_Vector_setElement (q, true, s)) ;

    // descriptor: invert the mask for vxm, and clear output before assignment
    LAGRAPH_OK (GrB_Descriptor_new (&desc)) ;
    LAGRAPH_OK (GrB_Descriptor_set (desc, GrB_MASK, GrB_SCMP)) ;
    LAGRAPH_OK (GrB_Descriptor_set (desc, GrB_OUTP, GrB_REPLACE)) ;

    //--------------------------------------------------------------------------
    // BFS traversal and label the nodes
    //--------------------------------------------------------------------------

    bool successor = true ; // true when some successor found
    for (int64_t level = 1 ; successor && level <= max_level ; level++)
    {
        // v<q> = level, using vector assign with q as the mask
        LAGRAPH_OK (GrB_assign (v, q, NULL, level, GrB_ALL, n, NULL)) ;

        // q<!v> = q ||.&& A ; finds all the unvisited
        // successors from current q, using !v as the mask
        LAGRAPH_OK (GrB_vxm (q, v, NULL, LAGraph_LOR_LAND_BOOL, q, A, desc)) ;

        // note that if A has no explicit zeros, then this works faster:
        // GrB_Vector_nvals (&nvals, q) ; successor = (nvals > 0) ; 

        // successor = ||(q)
        LAGRAPH_OK (GrB_reduce (&successor, NULL, LAGraph_LOR_MONOID, q, NULL));
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*v_output) = v ;       // return result
    v = NULL ;              // set to NULL so LAGRAPH_FREE_ALL doesn't free it
    LAGRAPH_FREE_ALL ;      // free all workspace (except for result v)
    return (GrB_SUCCESS) ;
}

