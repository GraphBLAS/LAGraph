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

#define LG_FREE_WORK                            \
{                                               \
    GrB_free (&frontier) ;                      \
    GrB_free (&J_vec) ;                         \
    GrB_free (&I_vec) ;                         \
    GrB_free (&J_matrix) ;                      \
    GrB_free (&I_matrix) ;                      \
    GrB_free (&Fd1A) ;                          \
    GrB_free (&paths) ;                         \
    GrB_free (&bc_vertex_flow) ;                     \
    GrB_free (&temp_update) ;                   \
    GrB_free (&Add_One_Divide) ;                \
    GrB_free (&Update) ;                     \
    if (Search != NULL)                              \
    {                                           \
        for (int64_t i = 0 ; i < n ; i++)       \
        {                                       \
            GrB_free (&(Search [i])) ;               \
        }                                       \
        LAGraph_Free ((void **) &Search, NULL) ;     \
    }                                           \
}

#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    GrB_free (centrality) ;         \
}

#include "LG_internal.h"
#include <LAGraphX.h>

#undef  LAGRAPH_CATCH
#define LAGRAPH_CATCH(status)                                           \
{                                                                       \
    print ("LAGraph failure (file %s, line %d): status: %d",     \
        __FILE__, __LINE__, status) ;                                   \
    LG_ERROR_MSG ("LAGraph failure (file %s, line %d): status: %d",     \
        __FILE__, __LINE__, status) ;                                   \
    LG_FREE_ALL ;                                                       \
    return (status) ;                                                   \
}

#undef GRB_CATCH
#define GRB_CATCH(info)                                                 \
{                                                                       \
    printf ("GraphBLAS failure (file %s, line %d): info: %d",     \
        __FILE__, __LINE__, info) ;                                     \
    LG_ERROR_MSG ("GraphBLAS failure (file %s, line %d): info: %d",     \
        __FILE__, __LINE__, info) ;                                     \
    LG_FREE_ALL ;                                                       \
    return (info) ;                                                     \
}

//------------------------------------------------------------------------------
// (1+x)/y function for double: z = (1 + x) / y
//------------------------------------------------------------------------------

void add_one_divide_function (double *z, const double *x, const double *y)
{
    double a = (*(x)) ;
    double b = (*(y)) ;
    (*(z)) = (1 + a) / b ;
}

#define ADD_ONE_DIVIDE_FUNCTION_DEFN                                           \
"void add_one_divide_function (double *z, const double *x, const double *y)\n" \
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
    char *msg
)
{

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

    GrB_Index n = 0 ;                   // # nodes in the graph

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

    // =========================================================================
    // === initialization =====================================================
    // =========================================================================

    GRB_TRY (GxB_BinaryOp_new (&Add_One_Divide,
        (GxB_binary_function) add_one_divide_function,
        GrB_FP64, GrB_FP64, GrB_FP64,
        "add_one_divide_function", ADD_ONE_DIVIDE_FUNCTION_DEFN)) ;

    // Initialize the frontier, paths, Update, and bc_vertex_flow
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Vector_new (&paths,    GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new (&frontier, GrB_FP64, n)) ;
    GRB_TRY (GrB_Matrix_new (&Update, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_Vector_new (&bc_vertex_flow, GrB_FP64, n)) ;

    
    // Initialize centrality matrix with zeros using A as structural mask
    LG_TRY (GrB_Matrix_new(centrality, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_assign (*centrality, A, NULL, 0.0, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;

    // Allocate memory for the array of S vectors
    LG_TRY (LAGraph_Calloc ((void **) &Search, n+1, sizeof (GrB_Vector), msg)) ;

    // =========================================================================
    // === Breadth-first search stage ==========================================
    // =========================================================================

    GrB_Index frontier_size, last_frontier_size = 0 ;
    GRB_TRY (GrB_Vector_nvals (&frontier_size, frontier)) ;

    int64_t depth, root ;
    for (root = 0 ; root < n ; root++)
    {
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
        GRB_TRY (GrB_Vector_new (&bc_vertex_flow, GrB_FP64, n)) ;
        GRB_TRY (GrB_assign(bc_vertex_flow, NULL, NULL, 0.0, GrB_ALL, n, NULL)) ;

        GRB_TRY (GrB_Vector_new(&J_vec, GrB_FP64, n)) ;
        GRB_TRY (GrB_Vector_new (&I_vec, GrB_FP64, n)) ;
        GRB_TRY (GrB_Matrix_new (&Fd1A, GrB_FP64, n, n)) ;
        GRB_TRY (GrB_Vector_new(&temp_update, GrB_FP64, n)) ; // Create a temporary vector

        // Backtrack through the BFS and compute centrality updates for each vertex
        // GrB_Index fd1_size;

        // printf ("\n----------------------------- backtrack:\n") ;

        while (depth >= 1)
        {        
            // printf ("\n----------------------------- backtrack depth : %" PRId64 "\n", depth) ;
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

            GRB_TRY(GrB_mxm(Fd1A, NULL, NULL, LAGraph_plus_first_fp64,
                I_matrix, A, NULL)) ;

            GRB_TRY(GrB_mxm(Update, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64,
                Fd1A, J_matrix, NULL)) ;

            //----------------------------------------------------------------------
            // centrality<A> += Update
            // Accumulate centrality values for edges
            //----------------------------------------------------------------------

            GRB_TRY (GrB_assign(*centrality, A, GrB_PLUS_FP64, Update, GrB_ALL, n, GrB_ALL, n, 
            GrB_DESC_S)) ;

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

    
    // =========================================================================
    // === finalize the centrality =============================================
    // =========================================================================

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
