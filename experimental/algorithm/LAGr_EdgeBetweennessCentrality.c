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

// This is an Advanced algorithm (G->AT is required).


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
    GrB_free (&bc_update) ;                     \
    GrB_free (&temp_update) ;                   \
    GrB_free (&Add_One_Divide) ;                \
    GrB_free (&v) ;                     \
    GrB_free (&U) ;                     \
    if (S != NULL)                              \
    {                                           \
        for (int64_t i = 0 ; i < n ; i++)       \
        {                                       \
            GrB_free (&(S [i])) ;               \
        }                                       \
        LAGraph_Free ((void **) &S, NULL) ;     \
    }                                           \
}

#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    GrB_free (&centrality_temp) ;               \
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

void add_one_divide_function (void *z, const void *x, const void *y)
{
    double a = (*((double *) x)) ;
    double b = (*((double *) y)) ;
    (*((double *) z)) = (1 + a) / b ;
}

#define ADD_ONE_DIVIDE_FUNCTION_DEFN                                           \
"void add_one_divide_function (void *z, const void *x, const void *y)      \n" \
"{                                                                         \n" \
"    double a = (*((double *) x)) ;                                        \n" \
"    double b = (*((double *) y)) ;                                        \n" \
"    (*((double *) z)) = (1 + a) / b ;                                     \n" \
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
    GrB_Vector paths = NULL ;

    // the delta vector for each node for each starting node.  A dense matrix.
    GrB_Vector bc_update = NULL ;

    // Update matrix for betweenness centrality, values for each node for
    // each starting node.  A dense matrix.
    GrB_Matrix U = NULL ;

    GrB_Vector v = NULL ;

    GrB_BinaryOp Add_One_Divide = NULL ;

    GrB_Matrix centrality_temp = NULL;    

    GrB_Vector J_vec = NULL ;
    GrB_Vector I_vec = NULL ;
    GrB_Matrix I_matrix = NULL ;
    GrB_Matrix J_matrix = NULL ;
    GrB_Matrix Fd1A = NULL ;
    GrB_Vector temp_update = NULL ; 

    // Temporary workspace matrix (sparse).
    // GrB_Matrix W = NULL ;

    GrB_Index n = 0 ;                   // # nodes in the graph

    LG_ASSERT (centrality != NULL, GrB_NULL_POINTER) ;
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

    printf ("G->A is:\n")  ;GxB_print (A, 5) ;
    printf ("G->AT is:\n") ;GxB_print (AT, 5) ;

    // =========================================================================
    // === initialization =====================================================
    // =========================================================================

    GRB_TRY (GxB_BinaryOp_new (&Add_One_Divide, add_one_divide_function,
        GrB_FP64, GrB_FP64, GrB_FP64,
        "add_one_divide_function", ADD_ONE_DIVIDE_FUNCTION_DEFN)) ;

    // Initialize paths and frontier with source notes
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Vector_new (&paths,    GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new (&frontier, GrB_FP64, n)) ;
    // GRB_TRY (LG_SET_FORMAT_HINT (paths, LG_BITMAP + LG_FULL)) ;

    GRB_TRY (GrB_Matrix_new (&U, GrB_FP64, n, n)) ;

    GRB_TRY (GrB_Vector_new (&v, GrB_FP64, n)) ;

    GRB_TRY (GrB_Vector_new (&bc_update, GrB_FP64, n)) ;


    // FIXME: rename this:
    // Initialize centrality_temp matrix with zeros using A as structural mask
    LG_TRY (GrB_Matrix_new(&centrality_temp, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_assign (centrality_temp, A, NULL, 0.0, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;

    // Initial frontier: frontier<!paths>= frontier*A
    // GRB_TRY (GrB_vxv (frontier, paths, NULL, LAGraph_plus_first_fp64,
    //     frontier, A, GrB_DESC_RSC)) ;

    // Allocate memory for the array of S vectors
    LG_TRY (LAGraph_Calloc ((void **) &S, n+1, sizeof (GrB_Vector), msg)) ;

    // =========================================================================
    // === Breadth-first search stage ==========================================
    // =========================================================================

    GrB_Index frontier_size, last_frontier_size = 0 ;
    GRB_TRY (GrB_Vector_nvals (&frontier_size, frontier)) ;

    int64_t depth, root ;
    for (root = 0 ; root < n ; root++)
    {
        printf("----\n");
        printf("root: %ld \n", root) ;
        depth = 0 ;
//      GrB_free (&(S [0])) ;
//      LG_TRY (LAGraph_Vector_Structure (&(S [0]), frontier, msg)) ;

        // root frontier: S [0](root) = true
        GrB_free (&(S [0])) ;
        GRB_TRY (GrB_Vector_new(&(S [0]), GrB_BOOL, n)) ;
        GRB_TRY (GrB_Vector_setElement_BOOL(S [0], (bool) true, root)) ;

        // clear paths, and then set paths (root) = 1
        GRB_TRY (GrB_Vector_clear (paths)) ;
        GRB_TRY (GrB_Vector_setElement (paths, (double) 1.0, root)) ;

        GRB_TRY (GrB_Matrix_clear (U)) ;

        GRB_TRY (GrB_Vector_clear (v)) ;

        // Extract row root from A into frontier vector: frontier = AT(root,:)
        // GRB_TRY (GrB_Col_extract (frontier, NULL, NULL, A, GrB_ALL, n, root,
        //     GrB_DESC_T0)) ;
        GRB_TRY (GrB_Col_extract (frontier, NULL, NULL, AT, GrB_ALL, n, root,
            NULL)) ;
        GRB_TRY (GrB_Vector_nvals (&frontier_size, frontier)) ;

        GxB_print(frontier, 5) ;

        while (frontier_size != 0)
        {
            depth++ ;

            //----------------------------------------------------------------------
            // Accumulate path counts: paths += frontier
            //----------------------------------------------------------------------

            GRB_TRY (GrB_assign (paths, NULL, GrB_PLUS_FP64, frontier, GrB_ALL, n,
                NULL)) ;

            //----------------------------------------------------------------------
            // Add frontier to S: S(depth, :) = frontier
            //----------------------------------------------------------------------

            GrB_free (&(S [depth])) ;
            LG_TRY (LAGraph_Vector_Structure (&(S [depth]), frontier, msg)) ;

            //----------------------------------------------------------------------
            // Update frontier: frontier<!paths> = frontier*A
            //----------------------------------------------------------------------
            
            GRB_TRY (LG_SET_FORMAT_HINT (frontier, LG_SPARSE)) ;
            GRB_TRY (GrB_vxm (frontier, paths, NULL, GxB_PLUS_FIRST_FP64, frontier, 
                A, GrB_DESC_RSC )) ;

            //----------------------------------------------------------------------
            // Get size of current frontier: frontier_size = nvals(frontier)
            //----------------------------------------------------------------------

            last_frontier_size = frontier_size ;
            GRB_TRY (GrB_Vector_nvals (&frontier_size, frontier)) ;
        }

        // printf ("depth: %ld\n", depth) ;
        // for (int64_t d = 1 ; d <= depth ; d++)
        // {
        //     printf ("------------------- S [%ld]\n", d) ;
        //     GxB_print (S [d], 5) ;
        // }
   
        // printf("  after:\n") ;
        // GxB_print(frontier, 5) ;

        // GRB_TRY (GrB_free (&frontier)) ;


        // =========================================================================
        // === Betweenness centrality computation phase ============================
        // =========================================================================

        // bc_update = ones (n, n) ; a full matrix (and stays full)
        GRB_TRY (GrB_Vector_new (&bc_update, GrB_FP64, n)) ;
        GRB_TRY (GrB_assign(bc_update, NULL, NULL, 0.0, GrB_ALL, n, NULL)) ;
        // // W: empty n-by-n array, as workspace
        // GRB_TRY (GrB_Matrix_new (&W, GrB_FP64, n, n)) ;

        // GxB_print (bc_update, 5) ;


        GRB_TRY (GrB_Vector_new(&J_vec, GrB_FP64, n)) ;
        GRB_TRY (GrB_Vector_new (&I_vec, GrB_FP64, n)) ;
        GRB_TRY (GrB_Matrix_new (&Fd1A, GrB_FP64, n, n)) ;
        GRB_TRY (GrB_Vector_new(&temp_update, GrB_FP64, n)) ; // Create a temporary vector

        // Backtrack through the BFS and compute centrality updates for each vertex
        // GrB_Index fd1_size;
        while (depth >= 1)
        {        
            printf ("backtrack depth %ld\n", depth) ;
            GrB_Vector f_d = S[depth] ;
            GrB_Vector f_d1 = S[depth - 1] ;

            printf("#####################################\n") ;
            printf("predecessors (frontier):\n") ;
            GxB_print (f_d, 5) ;
            GxB_print (f_d1, 5) ;

            // 18 w = S(d, :) ÷ p × v + S(d, :)
            // 19 U = A .× w
            // 20 w = S(d − 1, :) × p
            // 21 U = w .× U

            // make J Matrix
            GRB_TRY (GrB_eWiseMult(J_vec, f_d, NULL, Add_One_Divide, bc_update, paths, GrB_DESC_RS)) ;
            // GRB_TRY (GrB_eWiseMult(J_vec, f_d, NULL, GrB_PLUS_FP64, bc_update, paths, GrB_DESC_RS)) ;
            GRB_TRY (GrB_Matrix_diag(&J_matrix, J_vec, 0)) ;

            // make I matrix
            GRB_TRY (GrB_Vector_extract (I_vec, f_d1, NULL, paths, GrB_ALL, n, GrB_DESC_RS)) ;
            GRB_TRY (GrB_Matrix_diag(&I_matrix, I_vec, 0)) ;

            // combine

            // intermediate matrix for Fd1 * A
            // GRB_TRY (GrB_eWiseMult(Fd1A, NULL, NULL, GrB_TIMES_FP64, I_matrix, AT, NULL)) ;
            GRB_TRY(GrB_mxm(Fd1A, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64,
                I_matrix, A, NULL)) ;
            // GxB_print (Fd1A, 5) ;

            // GRB_TRY (GrB_eWiseMult(U, NULL, NULL, GrB_TIMES_FP64, Fd1A, J_matrix, NULL)) ;
            GRB_TRY(GrB_mxm(U, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64,
                Fd1A, J_matrix, NULL)) ;
            GxB_print (U, 5) ;
            

            // GxB_print (centrality_temp, 5) ;

            // 22 B = B + U
            // GRB_TRY (GrB_assign(centrality_temp, centrality_temp, GrB_PLUS_FP64, U, GrB_ALL, n, GrB_ALL, n, 
            // GrB_DESC_S)) ;
            GRB_TRY (GrB_eWiseAdd (centrality_temp, NULL, NULL, GrB_PLUS_FP64, centrality_temp, U, NULL)) ;


            // 23 v = U +.

            // Reduce "update" matrix to a vector (sum each column)
            GRB_TRY (GrB_reduce(temp_update, NULL, NULL, GrB_PLUS_MONOID_FP64, U, NULL)) ;
            // GRB_TRY (GrB_reduce(temp_update, NULL, NULL, GxB_ANY_FP64_MONOID, U, NULL)) ;
            GRB_TRY (GrB_eWiseAdd(bc_update, NULL, NULL, GrB_PLUS_FP64, bc_update, temp_update, NULL)) ;

            printf("#####################################\n sigma:\n") ;
            GxB_print (paths, 5) ;
            printf("delta:\n") ;
            GxB_print (bc_update, 5) ;
            printf("betweenness:\n") ;
            GxB_print (centrality_temp, 5) ;
            printf("#####################################\n") ;

            // Grb_reduce_monoid

            // 24 d = d − 1
            depth-- ;
        }

   
    }

    
    // =========================================================================
    // === finalize the centrality =============================================
    // =========================================================================

    GxB_print(centrality_temp, GxB_COMPLETE) ;
    
    *centrality = centrality_temp;

    // Initialize the centrality array with -n to avoid counting
    // zero length paths
    // GRB_TRY (GrB_Vector_new (centrality, GrB_FP64, n)) ;
    // GRB_TRY (GrB_assign (*centrality, NULL, NULL, -n, GrB_ALL, n, NULL)) ;

    // // centrality (i) += sum (bc_update (:,i)) for all nodes i
    // GRB_TRY (GrB_reduce (*centrality, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64,
    //     bc_update, GrB_DESC_T0)) ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
