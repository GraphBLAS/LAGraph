//------------------------------------------------------------------------------
// LAGr_EdgeBetweennessCentrality: edge betweenness-centrality
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-Licene-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICEnE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Casey and Tim Davis, Texas A&M University;
// Adapted and revised from GraphBLAS C API Spec, Appendix B.4.

//------------------------------------------------------------------------------

// LAGr_EdgeBetweennessCentrality: Batch algorithm for computing
// betweeness centrality, using push-pull optimization.

// This is an Advanced algorithm (G->AT is required).


//------------------------------------------------------------------------------

#define LG_FREE_WORK                            \
{                                               \
    GrB_free (&frontier) ;                      \
    GrB_free (&paths) ;                         \
    GrB_free (&bc_update) ;                     \
    GrB_free (&W) ;                             \
    if (S != NULL)                              \
    {                                           \
        for (int64_t i = 0 ; i < n ; i++)       \
        {                                       \
            if (S [i] == NULL) break ;          \
            GrB_free (&(S [i])) ;               \
        }                                       \
        LAGraph_Free ((void **) &S, NULL) ;     \
    }                                           \
}

#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    GrB_free (centrality) ;         \
}

#include "LG_internal.h"
#include <LAGraphX.h>

GrB_Info add_one_divide(double *z, const double *x, const double *y) {
    *z = (1 + *x) / *y;  // Operation (1 + x) / y
    return GrB_SUCCESS;
}

//------------------------------------------------------------------------------
// LAGr_EdgeBetweennessCentrality: edge betweenness-centrality
//------------------------------------------------------------------------------

int LAGr_EdgeBetweennessCentrality
(
    // output:
    GrB_Vector *centrality,     // centrality(i): betweeness centrality of i
    // input:
    LAGraph_Graph G,            // input graph
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    // Array of BFS search matrices.
    // S [i] is a sparse matrix that stores the depth at which each vertex is
    // first seen thus far in each BFS at the current depth i. Each column
    // corresponds to a BFS traversal starting from a source node.
    GrB_Vector *S = NULL ;

    // Frontier matrix, a sparse matrix.
    // Stores # of shortest paths to vertices at current BFS depth
    GrB_Vector frontier = NULL ;

    // Paths matrix holds the number of shortest paths for each node and
    // starting node discovered so far.  A dense matrix that is updated with
    // sparse updates, and also used as a mask.
    GrB_Matrix paths = NULL ;

    // Update matrix for betweenness centrality, values for each node for
    // each starting node.  A dense matrix.
    GrB_Matrix bc_update = NULL ;

    // Temporary workspace matrix (sparse).
    GrB_Matrix W = NULL ;

    GrB_Index n = 0 ;                   // # nodes in the graph

    LG_ASSERT (centrality != NULL && sources != NULL, GrB_NULL_POINTER) ;
    (*centrality) = NULL ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;

    GrB_Matrix A = G->A ;
    GrB_Matrix AT ;
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
        G->is_symmetric_structure == LAGraph_TRUE)
    {
        // A and A' have the same structure
        AT = A ;
    }
    else
    {
        // A and A' differ
        AT = G->AT ;
        LG_ASSERT_MSG (AT != NULL, LAGRAPH_NOT_CACHED, "G->AT is required") ;
    }

    // =========================================================================
    // === initialization =====================================================
    // =========================================================================

    // Initialize paths and frontier with source notes
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_new (&paths,    GrB_FP64, n, n)) ;
    GRB_TRY (GrB_Vector_new (&frontier, GrB_FP64, n)) ; // todo: change to bool?
    GRB_TRY (LG_SET_FORMAT_HINT (paths, LG_BITMAP + LG_FULL)) ;

    // TOOO: change this to just set everything to 1 instead of looping
    for (GrB_Index i = 0 ; i < n ; i++)
    {
        // paths (i,s(i)) = 1
        // frontier (i,s(i)) = 1
        double one = 1 ;
        GrB_Index src = sources [i] ;
        GRB_TRY (GrB_Matrix_setElement (paths,    one, i, src)) ;
        GRB_TRY (GrB_Matrix_setElement (frontier, one, i, src)) ;
    }

    // Initial frontier: frontier<!paths>= frontier*A
    GRB_TRY (GrB_mxm (frontier, paths, NULL, LAGraph_plus_first_fp64,
        frontier, A, GrB_DESC_RSC)) ;

    // Allocate memory for the array of S matrices
    LG_TRY (LAGraph_Malloc ((void **) &S, n+1, sizeof (GrB_Vector), msg)) ;
    S [0] = NULL ;

    // =========================================================================
    // === Breadth-first search stage ==========================================
    // =========================================================================

    bool last_was_pull = false ;
    GrB_Index frontier_size, last_frontier_size = 0 ;
    GRB_TRY (GrB_Vector_nvals (&frontier_size, frontier)) ;

    int64_t depth ;
    for (depth = 0 ; frontier_size > 0 && depth < n ; depth++)
    {

        //----------------------------------------------------------------------
        // S [depth] = structure of frontier
        //----------------------------------------------------------------------

        S [depth+1] = NULL ;
        LG_TRY (LAGraph_Vector_Structure (&(S [depth]), frontier, msg)) ;

        //----------------------------------------------------------------------
        // Accumulate path counts: paths += frontier
        //----------------------------------------------------------------------

        GRB_TRY (GrB_assign (paths, NULL, GrB_PLUS_FP64, frontier, GrB_ALL, n,
            GrB_ALL, n, NULL)) ;

        //----------------------------------------------------------------------
        // Update frontier: frontier<!paths> = frontier*A
        //----------------------------------------------------------------------

        // pull if frontier is more than 10% dense,
        // or > 6% dense and last step was pull
        double frontier_denity = ((double) frontier_size) / (double) (n*n) ;
        bool do_pull = frontier_denity > (last_was_pull ? 0.06 : 0.10 ) ;

        if (do_pull)
        {
            // frontier<!paths> = frontier*AT'
            GRB_TRY (LG_SET_FORMAT_HINT (frontier, LG_BITMAP)) ;
            GRB_TRY (GrB_mxm (frontier, paths, NULL, LAGraph_plus_first_fp64,
                frontier, AT, GrB_DESC_RSCT1)) ;
        }
        else // push
        {
            // frontier<!paths> = frontier*A
            GRB_TRY (LG_SET_FORMAT_HINT (frontier, LG_SPARSE)) ;
            GRB_TRY (GrB_mxm (frontier, paths, NULL, LAGraph_plus_first_fp64,
                frontier, A, GrB_DESC_RSC)) ;
        }

        //----------------------------------------------------------------------
        // Get size of current frontier: frontier_size = nvals(frontier)
        //----------------------------------------------------------------------

        last_frontier_size = frontier_size ;
        last_was_pull = do_pull ;
        GRB_TRY (GrB_Matrix_nvals (&frontier_size, frontier)) ;
    }

    GRB_TRY (GrB_free (&frontier)) ;

    // =========================================================================
    // === Betweenness centrality computation phase ============================
    // =========================================================================

    // bc_update = ones (n, n) ; a full matrix (and stays full)
    GRB_TRY (GrB_Matrix_new (&bc_update, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_assign (bc_update, NULL, NULL, 1, GrB_ALL, n, GrB_ALL, n,
        NULL)) ;
    // W: empty n-by-n array, as workspace
    GRB_TRY (GrB_Matrix_new (&W, GrB_FP64, n, n)) ;

    // Backtrack through the BFS and compute centrality updates for each vertex
    for (int64_t d = depth-1 ; d > 0 ; d--)
    {
        GrB_Vector f_d = S[d] ;
        GrB_Vector f_d1 = S[d - 1] ;

        // 18 w = S(d, :) ÷ p × v + S(d, :)
        // 19 U = A .× w
        // 20 w = S(d − 1, :) × p
        // 21 U = w .× U

        // make J Matrix

        GrB_Vector J_vec;
        GRB_TRY (GrB_Vector_new(&J_vec, GrB_FP64, n));

        GrB_Vector bc_update_col, paths_col;
        GRB_TRY(GrB_Vector_new(&bc_update_col, GrB_FP64, rows));
        GRB_TRY(GrB_Vector_new(&paths_col, GrB_FP64, rows));
        GRB_TRY(GrB_extract(bc_update_col, f_d, GrB_NULL, bc_update, d));
        GRB_TRY(GrB_extract(paths_col, f_d, GrB_NULL, paths, d));

        GRB_TRY(GrB_eWiseAdd(J_vec, f_d, NULL, add_one_divide, bc_update_col, paths_col, GrB_DESC_R));

        GrB_Matrix J_matrix;
        GRB_TRY(GrB_Matrix_diag(&J_matrix, J_vec, 0));


        // make I matrix

        GrB_Vector I_vec;
        GRB_TRY (GrB_Vector_new (&I_vec, GrB_FP64, n)) ;
        GRB_TRY (GrB_extract (I_vec, f_d1, GrB_NULL, paths, d-1)) ;

        GrB_Matrix I_matrix;
        GRB_TRY(GrB_Matrix_diag(&I_matrix, I_vec, 0));


        // combine

        // intermediate matrix for Fd1 * A
        GrB_Matrix Fd1A;
        GrB_Matrix_new(&Fd1A, GrB_FP64, n, n);
        GRB_TRY(GrB_eWiseMult(Fd1A, NULL, NULL, GrB_TIMES_FP64, J_matrix, A, NULL));

        GrB_Matrix update;
        GrB_Matrix_new(update, GrB_FP64, n, n);
        GRB_TRY(GrB_eWiseMult(*update, NULL, NULL, GrB_TIMES_FP64, Fd1A, I_matrix, NULL));

        // free intermediate matrix
        GrB_Matrix_free(&Fd1A);


        // 22 B = B + U
        GRB_TRY(GrB_eWiseAdd(centrality, NULL, NULL, GrB_PLUS_FP64, centrality, update, NULL));
        
        // 23 v = U +.
        GrB_Vector temp_update; 
        GrB_Vector_new(&temp_update, GrB_FP64, n); // Create a temporary vector

        // Reduce "update" matrix to a vector (sum each column)
        GRB_TRY(GrB_reduce(temp_update, NULL, NULL, GrB_PLUS_MONOID_FP64, update, NULL));
        GRB_TRY(GrB_eWiseAdd(bc_update, NULL, NULL, GrB_PLUS_FP64, bc_update, temp_update, NULL));

        GrB_Vector_free(&temp_update);
        
        // 24 d = d − 1
    }

    // =========================================================================
    // === finalize the centrality =============================================
    // =========================================================================

    // Initialize the centrality array with -n to avoid counting
    // zero length paths
    GRB_TRY (GrB_Vector_new (centrality, GrB_FP64, n)) ;
    GRB_TRY (GrB_assign (*centrality, NULL, NULL, -n, GrB_ALL, n, NULL)) ;

    // centrality (i) += sum (bc_update (:,i)) for all nodes i
    GRB_TRY (GrB_reduce (*centrality, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64,
        bc_update, GrB_DESC_T0)) ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
