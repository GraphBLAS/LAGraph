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

// Contributed by Casey Pei and Tim Davis, Texas A&M University;
// Adapted and revised from GraphBLAS C API Spec, Appendix B.4.

//------------------------------------------------------------------------------

// LAGr_EdgeBetweennessCentrality: Exact algorithm for computing
// betweeness centrality.

// This is an Advanced algorithm (no self edges allowed)

//------------------------------------------------------------------------------

#define useAssign
// #define debug

#define LG_FREE_WORK                                \
{                                                   \
    GrB_free (&frontier) ;                          \
    GrB_free (&J_vec) ;                             \
    GrB_free (&I_vec) ;                             \
    GrB_free (&J_matrix) ;                          \
    GrB_free (&I_matrix) ;                          \
    GrB_free (&Fd1A) ;                              \
    GrB_free (&paths) ;                             \
    GrB_free (&bc_vertex_flow) ;                    \
    GrB_free (&temp_update) ;                       \
    GrB_free (&Add_One_Divide) ;                    \
    GrB_free (&Update) ;                            \
    GrB_free (&HalfUpdate) ;                        \
    GrB_free (&HalfUpdateT) ;                       \
    GrB_free (&SymmetricUpdate) ;                   \
    GrB_free (&internal_sources) ;                  \
    if (Search != NULL)                             \
    {                                               \
        for (int64_t i = 0 ; i <= n ; i++)          \
        {                                           \
            GrB_free (&(Search [i])) ;              \
        }                                           \
        LAGraph_Free ((void **) &Search, NULL) ;    \
    }                                               \
}

#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    GrB_free (centrality) ;         \
}

#include "LG_internal.h"
#include <LAGraphX.h>

//------------------------------------------------------------------------------
// (1+x)/y function for double: z = (1 + x) / y
//------------------------------------------------------------------------------

void LG_EBC_add_one_divide_function (double *z, const double *x, const double *y)
{
    double a = (*(x)) ;
    double b = (*(y)) ;
    (*(z)) = (1 + a) / b ;
}

#define ADD_ONE_DIVIDE_FUNCTION_DEFN                                           \
"void LG_EBC_add_one_divide_function (double *z, const double *x, const double *y)\n" \
"{                                                                         \n" \
"    double a = (*(x)) ;                                                   \n" \
"    double b = (*(y)) ;                                                   \n" \
"    (*(z)) = (1 + a) / b ;                                                \n" \
"}"

//------------------------------------------------------------------------------
// LAGr_EdgeBetweennessCentrality: edge betweenness-centrality
//------------------------------------------------------------------------------

int LAGr_EdgeBetweennessCentrality
(
    // output:
    GrB_Matrix *centrality,     // centrality(i): betweeness centrality of i
    // input:
    LAGraph_Graph G,            // input graph
    GrB_Vector sources,         // source vertices to compute shortest paths (if NULL or empty, use all vertices)
    char *msg
)
{

#if LAGRAPH_SUITESPARSE

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    // Array of BFS search matrices.
    // Search[i] is a sparse matrix that stores the depth at which each vertex is
    // first seen thus far in each BFS at the current depth i. Each column
    // corresponds to a BFS traversal starting from a source node.
    GrB_Vector *Search = NULL ;

    // Frontier vector, a sparse matrix.
    // Stores # of shortest paths to vertices at current BFS depth
    GrB_Vector frontier = NULL ;

    // Paths matrix holds the number of shortest paths for each node and
    // starting node discovered so far. A dense vector that is updated with
    // sparse updates, and also used as a mask.
    GrB_Vector paths = NULL ;

    // The betweenness centrality for each vertex. A dense vector that
    // accumulates flow values during backtracking.
    GrB_Vector bc_vertex_flow = NULL ;

    // Update matrix for betweenness centrality for each edge. A sparse matrix
    // that holds intermediate centrality updates.
    GrB_Matrix Update = NULL ;

    // Binary operator for computing (1+x)/y in centrality calculations
    GrB_BinaryOp Add_One_Divide = NULL ;

    // Temporary vectors and matrices for intermediate calculations
    // Diagonal values for J_matrix
    GrB_Vector J_vec = NULL ;      

    // Diagonal values for I_matrix
    GrB_Vector I_vec = NULL ;      

    // Matrix for previous level contributions
    GrB_Matrix I_matrix = NULL ;   

    // Matrix for current level contributions
    GrB_Matrix J_matrix = NULL ;  
    
    // Intermediate product matrix
    GrB_Matrix Fd1A = NULL ;       

    // Temporary vector for centrality updates
    GrB_Vector temp_update = NULL ;

    // Temporary matrices for doing updates on
    // approximate and undirected graphs
    GrB_Matrix HalfUpdate = NULL ;
    GrB_Matrix HalfUpdateT = NULL ;
    GrB_Matrix SymmetricUpdate = NULL ;

    // Source nodes vector (will be created if NULL is passed)
    GrB_Vector internal_sources = NULL;
    bool created_sources = false;

    GrB_Index n = 0 ;                   // # nodes in the graph

    double t1_total = 0;
    double t2_total = 0;
    double t3_total = 0;

    LG_ASSERT (centrality != NULL, GrB_NULL_POINTER) ;
    (*centrality) = NULL ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;

    GrB_Matrix A = G->A ;
    #if 0
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
    #endif

    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GrB_Index nsources ;
    if (sources == NULL)
    {
        nsources = n ;
    }
    else
    {
        GRB_TRY (GrB_Vector_nvals (&nsources, sources)) ;
    }
    LG_ASSERT (nsources > 0, GrB_INVALID_VALUE) ;

    // =========================================================================
    // === initialization =====================================================
    // =========================================================================

    GRB_TRY (GxB_BinaryOp_new (&Add_One_Divide,
        (GxB_binary_function) LG_EBC_add_one_divide_function,
        GrB_FP64, GrB_FP64, GrB_FP64,
        "LG_EBC_add_one_divide_function", ADD_ONE_DIVIDE_FUNCTION_DEFN)) ;

    // Initialize the frontier, paths, Update, and bc_vertex_flow
    GRB_TRY (GrB_Vector_new (&paths,    GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new (&frontier, GrB_FP64, n)) ;
    GRB_TRY (GrB_Matrix_new (&Update, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_Vector_new (&bc_vertex_flow, GrB_FP64, n)) ;

    
    // Initialize centrality matrix with zeros using A as structural mask
    LG_TRY (GrB_Matrix_new(centrality, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_assign (*centrality, A, NULL, 0.0, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;

    // Allocate memory for the array of S vectors
    LG_TRY (LAGraph_Calloc ((void **) &Search, n + 1, sizeof (GrB_Vector), msg)) ;

    // =========================================================================
    // === Process source nodes ================================================
    // =========================================================================

    // If sources is NULL, create a dense vector with all vertices
    if (sources == NULL)
    {
        // Create a vector with all nodes as sources
        GRB_TRY (GrB_Vector_new (&internal_sources, GrB_INT64, n)) ;

        // internal_sources (0:n-1) = 0
        GRB_TRY (GrB_assign (internal_sources, NULL, NULL, 0, GrB_ALL, n, NULL)) ;

        // internal_sources (0:n-1) = 0:n-1
        GRB_TRY (GrB_apply (internal_sources, NULL, NULL, GrB_ROWINDEX_INT64,
            internal_sources, 0, NULL)) ;

        /*
        int64_t ns = n;
        for (GrB_Index i = 0; i < ns; i++)
        {
            GRB_TRY (GrB_Vector_setElement_INT64 (internal_sources, i, i)) ;
        }
        */

        // Use this vector instead
        sources = internal_sources;
        created_sources = true;
    }

    // =========================================================================
    // === Breadth-first search stage ==========================================
    // =========================================================================

    GrB_Index frontier_size, last_frontier_size = 0 ;
    GRB_TRY (GrB_Vector_nvals (&frontier_size, frontier)) ;

    int64_t depth;
    GrB_Index root;

    GRB_TRY (GrB_Vector_new(&J_vec, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new (&I_vec, GrB_FP64, n)) ;
    GRB_TRY (GrB_Matrix_new (&Fd1A, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_Vector_new(&temp_update, GrB_FP64, n)) ; // Create a temporary vector

    GRB_TRY (GrB_Matrix_new(&HalfUpdate, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_Matrix_new(&HalfUpdateT, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_Matrix_new(&SymmetricUpdate, GrB_FP64, n, n)) ;

    // Iterate through source nodes
    for (GrB_Index i = 0; i < nsources; i++)
    {
        GRB_TRY (GrB_Vector_extractElement(&root, sources, i)) ;
        
        // Verify the root index is valid
        LG_ASSERT (root < n, GrB_INVALID_VALUE) ;

        depth = 0 ;

        // root frontier: Search [0](root) = true
        GrB_free (&(Search [0])) ;
        GRB_TRY (GrB_Vector_new(&(Search [0]), GrB_BOOL, n)) ;
        GRB_TRY (GrB_Vector_setElement_BOOL(Search [0], (bool) true, root)) ;

        // clear paths, and then set paths (root) = 1
        GRB_TRY (GrB_Vector_clear (paths)) ;
        GRB_TRY (GrB_Vector_setElement (paths, (double) 1.0, root)) ;

        GRB_TRY (GrB_Matrix_clear (Update)) ;

        // Extract row root from A into frontier vector: frontier = A(root,:)
        GRB_TRY (GrB_Col_extract (frontier, NULL, NULL, A, GrB_ALL, n, root,
            GrB_DESC_T0)) ;

        GRB_TRY (GrB_Vector_nvals (&frontier_size, frontier)) ;
        GRB_TRY (GrB_assign (frontier, frontier, NULL, 1.0, GrB_ALL, n, GrB_DESC_S)) ;

        while (frontier_size != 0)
        {
            depth++ ;

            //----------------------------------------------------------------------
            // paths += frontier
            // Accumulate path counts for vertices at current depth
            //----------------------------------------------------------------------

            GRB_TRY (GrB_assign (paths, NULL, GrB_PLUS_FP64, frontier, GrB_ALL, n,
                NULL)) ;

            //----------------------------------------------------------------------
            // Search[depth] = structure(frontier)
            // Record the frontier structure at current depth
            //----------------------------------------------------------------------

            GrB_free (&(Search [depth])) ;
            LG_TRY (LAGraph_Vector_Structure (&(Search [depth]), frontier, msg)) ;

            //----------------------------------------------------------------------
            // frontier<!paths> = frontier * A
            //----------------------------------------------------------------------
            
            GRB_TRY (LG_SET_FORMAT_HINT (frontier, LG_SPARSE)) ;
            GRB_TRY (GrB_vxm (frontier, paths, NULL, /* LAGraph_plus_first_fp64 */
                GxB_PLUS_FIRST_FP64, frontier, 
                A, GrB_DESC_RSC )) ;

            //----------------------------------------------------------------------
            // Get size of current frontier: frontier_size = nvals(frontier)
            //----------------------------------------------------------------------

            last_frontier_size = frontier_size ;
            GRB_TRY (GrB_Vector_nvals (&frontier_size, frontier)) ;
        }


        // =========================================================================
        // === Betweenness centrality computation phase ============================
        // =========================================================================

        // bc_vertex_flow = ones (n, n) ; a full matrix (and stays full)
        GRB_TRY (GrB_assign(bc_vertex_flow, NULL, NULL, 0.0, GrB_ALL, n, NULL)) ;

        GRB_TRY (GrB_Matrix_clear (HalfUpdate)) ;
        GRB_TRY (GrB_Matrix_clear (HalfUpdateT)) ;
        GRB_TRY (GrB_Matrix_clear (SymmetricUpdate)) ;
        GRB_TRY (GrB_Matrix_clear (Fd1A)) ;
        GRB_TRY (GrB_Vector_clear (J_vec)) ;
        GRB_TRY (GrB_Vector_clear (I_vec)) ;
        GRB_TRY (GrB_Vector_clear (temp_update)) ;




        // Backtrack through the BFS and compute centrality updates for each vertex
        // GrB_Index fd1_size;

        while (depth >= 1)
        {        
            GrB_Vector f_d = Search [depth] ;
            GrB_Vector f_d1 = Search [depth - 1] ;

            //----------------------------------------------------------------------
            // j<S(depth, :)> = (1 + v) / p
            // J = diag(j)
            // Compute weighted contributions from current level
            //----------------------------------------------------------------------

            GRB_TRY (GrB_eWiseMult(J_vec, f_d, NULL, Add_One_Divide, bc_vertex_flow, paths, GrB_DESC_RS)) ;
            GRB_TRY (GrB_Matrix_diag(&J_matrix, J_vec, 0)) ;

            //----------------------------------------------------------------------
            // i<S(depth-1, :)> = p
            // I = diag(i)
            // Compute weighted contributions from previous level
            //----------------------------------------------------------------------

            GRB_TRY (GrB_Vector_extract (I_vec, f_d1, NULL, paths, GrB_ALL, n, GrB_DESC_RS)) ;
            GRB_TRY (GrB_Matrix_diag(&I_matrix, I_vec, 0)) ;

            //----------------------------------------------------------------------
            // Update = I × A × J 
            // Compute edge updates based on current level weights
            //----------------------------------------------------------------------

            double t1 = LAGraph_WallClockTime();
            GRB_TRY(GrB_mxm(Fd1A, NULL, NULL, LAGraph_plus_first_fp64,
                I_matrix, A, NULL));
            t1 = LAGraph_WallClockTime() - t1;
            t1_total += t1;

            double t2 = LAGraph_WallClockTime();
            GRB_TRY(GrB_mxm(Update, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64,
                Fd1A, J_matrix, NULL));
            t2 = LAGraph_WallClockTime() - t2;
            t2_total += t2;
            GRB_TRY (GrB_free (&I_matrix)) ;
            GRB_TRY (GrB_free (&J_matrix)) ;
            //----------------------------------------------------------------------
            // centrality<A> += Update
            // Accumulate centrality values for edges
            //----------------------------------------------------------------------

            #ifdef useAssign
                // centrality{A} += Update, using assign
                double t3 = LAGraph_WallClockTime();
                
                if (G->kind == LAGraph_ADJACENCY_UNDIRECTED) {
                    // First divide the Update matrix by 2 for symmetric distribution
                    GrB_apply(HalfUpdate, NULL, NULL, GrB_DIV_FP64, Update, 2.0, NULL);

                    // Create a transposed version of the update
                    GrB_transpose(HalfUpdateT, NULL, NULL, HalfUpdate, NULL);

                    // Add the original and transposed matrices to create a symmetric update
                    GrB_eWiseAdd(SymmetricUpdate, NULL, NULL, GrB_PLUS_FP64, HalfUpdate, HalfUpdateT, NULL);

                    // Apply the symmetric update to the centrality
                    GRB_TRY(GrB_assign(*centrality, A, GrB_PLUS_FP64, SymmetricUpdate, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S));

                }
                else {
                    GRB_TRY (GrB_assign(*centrality, A, GrB_PLUS_FP64, Update, GrB_ALL, n, GrB_ALL, n, 
                        GrB_DESC_S));
                }

                t3 = LAGraph_WallClockTime() - t3;
                t3_total += t3;
            #else
                // TODO: Approx update using ewise add not implemented
                // centrality = centrality + Update using eWiseAdd
                double t3 = LAGraph_WallClockTime();
                GRB_TRY (GrB_eWiseAdd (*centrality, NULL, NULL, GrB_PLUS_FP64, *centrality, Update, NULL));
                t3 = LAGraph_WallClockTime() - t3;
                t3_total += t3;
            #endif


            //----------------------------------------------------------------------
            // v = Update +.
            // Reduce update matrix to vector for next iteration
            //----------------------------------------------------------------------

            GRB_TRY (GrB_reduce(temp_update, NULL, NULL, GrB_PLUS_MONOID_FP64, Update, NULL)) ;
            GRB_TRY (GrB_eWiseAdd(bc_vertex_flow, NULL, NULL, GrB_PLUS_FP64, bc_vertex_flow, temp_update, NULL)) ;

            // 24 d = d − 1
            depth-- ;
        }
        
    }

    #ifdef debug
        printf("  I*A time: %g\n", t1_total);

        printf("  (I*A)*J time: %g\n", t2_total);

        #ifdef useAssign
            printf("  Centrality update using assign time: %g\n", t3_total);
        #else
            printf("  Centrality update using eWiseAdd time: %g\n", t3_total);
        #endif

        GxB_print (*centrality, GxB_FULL) ;

    #endif


    // =========================================================================
    // === finalize the centrality =============================================
    // =========================================================================

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
