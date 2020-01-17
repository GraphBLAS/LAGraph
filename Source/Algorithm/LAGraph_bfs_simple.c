//------------------------------------------------------------------------------
// LAGraph_bfs_simple:  simple breadth-first search
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

// LAGraph_bfs_simple:  based on the breadth-first search in the GraphBLAS C
// API Specification by Scott McMillan, CMU.  Modified by Tim Davis, Texas A&M.

// Perform a single-source BFS, starting at a source node.  Returns a dense
// vector v such that v(i) > 0 if node is reachable from the source node.
// v(source)=1 and v(i)=k if the path with the fewest edges from the source
// node to i has k-1 edges.  If i is not reachable from the source node, then
// v(i) is zero.

// This method is a simple one for illustration only, and works well in
// practice, except in the following cases:

// (1) It takes Omega(n) time.  If nvals(v) << n is expected, use
// LAGraph_bfs_pushpull instead, which is much faster if v is expected to be
// very sparse.

// (2) It assumes that vxm(q,A) is fast, and implemented in a 'push' fashion,
// using saxpy operations instead of dot products.  This requires that the rows
// A(i,:) are efficient to access, which is the case if A is in CSR format
// internally in GraphBLAS (or perhaps in both CSR and CSC formats).  This
// method will be *exceedlingly* slow if A is a MATLAB matrix in CSC format.

// See LAGraph_bfs_pushpull, which handles the above two cases.

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL    \
{                           \
    GrB_free (&v) ;         \
    GrB_free (&q) ;         \
}

GrB_Info LAGraph_bfs_simple     // push-only BFS
(
    GrB_Vector *v_output,   // v(i) is the BFS level of node i in the graph
    GrB_Matrix A,           // input graph, treated as if boolean in semiring
    GrB_Index source        // starting node of the BFS
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Vector q = NULL ;           // nodes visited at each level
    GrB_Vector v = NULL ;           // result vector
    if (v_output == NULL) LAGRAPH_ERROR ("argument missing", GrB_NULL_POINTER) ;
    GrB_Index n, nvals ;
    LAGr_Matrix_nrows (&n, A) ;

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS ) \
        && ( GxB_IMPLEMENTATION >= GxB_VERSION (3,2,0) )
    GrB_Descriptor desc_s  = GrB_DESC_S ;
    GrB_Descriptor desc_rc = GrB_DESC_RC ;
    GrB_Semiring semiring  = GxB_ANY_PAIR_BOOL ;
    #else
    GrB_Descriptor desc_s  = NULL ;
    GrB_Descriptor desc_rc = LAGraph_desc_oocr ;
    GrB_Semiring semiring  = LAGraph_LOR_FIRST_BOOL ;
    #endif

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    // create an empty vector v, and make it dense
    LAGr_Vector_new (&v, (n > INT32_MAX) ? GrB_INT64 : GrB_INT32, n) ;
    LAGr_assign (v, NULL, NULL, 0, GrB_ALL, n, NULL) ;

    // create a boolean vector q, and set q(source) to true
    LAGr_Vector_new (&q, GrB_BOOL, n) ;
    LAGr_Vector_setElement (q, true, source) ;

    //--------------------------------------------------------------------------
    // BFS traversal and label the nodes
    //--------------------------------------------------------------------------

    for (int64_t level = 1 ; level <= n ; level++)
    {
        // v<q> = level
        LAGr_assign (v, q, NULL, level, GrB_ALL, n, desc_s) ;

        // break if q is empty
        LAGr_Vector_nvals (&nvals, q) ;
        if (nvals == 0) break ;

        // q'<!v> = q'*A
        LAGr_vxm (q, v, NULL, semiring, q, A, desc_rc) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*v_output) = v ;       // return result
    v = NULL ;              // set to NULL so LAGRAPH_FREE_ALL doesn't free it
    LAGRAPH_FREE_ALL ;      // free all workspace (except for result v)
    return (GrB_SUCCESS) ;
}

