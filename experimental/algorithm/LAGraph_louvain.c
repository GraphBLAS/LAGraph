//------------------------------------------------------------------------------
// LAGraph_louvain.c: Louvain method
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2026 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Roi Lipman and Gabriel Gomez, FalkorDB

//------------------------------------------------------------------------------

#include <LAGraph.h>
#include <LAGraphX.h>
#include <LG_internal.h>

#undef  LG_FREE_ALL
#define LG_FREE_ALL          \
{                            \
    GrB_free (&V) ;          \
    GrB_free (&cont) ;       \
}

// construct C, C[i,:] = [i]
// and I, where I[i] is i's community id (readonly handle to C's col-indices)
static GrB_Info _initCommunities
(
    GrB_Matrix *C,        // C [i] = i
    uint64_t **I,         // C's indices array
    GrB_Index node_count  // node count
)
{
    //--------------------------------------------------------------------------
    // initialize C
    //--------------------------------------------------------------------------

    char *msg = NULL ;
    GrB_Vector V = NULL ;
    GxB_Container cont = NULL ;

    GRB_TRY (GrB_Vector_new (&V, GrB_BOOL, node_count)) ;
    GRB_TRY (GrB_assign (V, NULL, NULL, true, GrB_ALL, node_count, NULL)) ;
    GRB_TRY (GrB_Matrix_diag (C, V, 0)) ;
    GRB_TRY (GrB_free (&V)) ;

    // C should be CSR with 64-bit indices
    GRB_TRY (GrB_set (*C, 64,           GxB_COLINDEX_INTEGER_HINT   )) ;
    GRB_TRY (GrB_set (*C, GxB_SPARSE,   GxB_SPARSITY_CONTROL        )) ;
    GRB_TRY (GrB_set (*C, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT)) ;

    //--------------------------------------------------------------------------
    // get a handle to C's indices array
    //--------------------------------------------------------------------------

    GRB_TRY (GxB_Container_new (&cont)) ;
    GRB_TRY (GxB_unload_Matrix_into_Container (*C, cont, NULL)) ;

    int handling ;
    GrB_Type type ;
    uint64_t n, X_memsize ;

    GRB_TRY (GxB_Vector_unload (cont->i, (void**) I, &type, &n, &X_memsize,
        &handling, NULL)) ;

    // load I and mark it read-only
    GRB_TRY (GxB_Vector_load (cont->i, (void**) I, type, n, X_memsize,
        GxB_IS_READONLY, NULL)) ;

    GRB_TRY (GxB_load_Matrix_from_Container (*C, cont, NULL)) ;
    GRB_TRY (GrB_free (&cont)) ;

    LG_FREE_ALL ;
    return GrB_SUCCESS ;
}

// computes the modularity contribution of attaching a node of degree i_degree,
// with kin edges into a community, to that community.
static inline double gain
(
    uint64_t i_degree,
    uint64_t kjn,
    uint64_t sigma_tot,
    double M2
)
{
    return (2.0 * (double) kjn) / M2 -
        ((double) sigma_tot * (double) i_degree) / (M2 * M2) ;
}

#undef  LG_FREE_WORK
#define LG_FREE_WORK                             \
{                                                \
    GrB_free (&desc) ;                           \
    GrB_free (&c_list) ;                         \
    GrB_free (&x) ;                              \
    GrB_free (&S_squished) ;                     \
    GrB_free (&S_new) ;                          \
    GrB_free (&A_squished) ;                     \
    GrB_free (&A_new) ;                          \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                              \
{                                                \
    if (free_A) GrB_free (&A) ;                  \
    LG_FREE_WORK ;                               \
}

// helper function: phase 2 of Louvain.
// Input C is node->community (possibly with empty community columns).
// Output S maps coarse nodes to previous level nodes, and A_new is coarsened.
static GrB_Info LG_Louvain_aggregate
(
    GrB_Matrix *S_handle,     // output: coarse->prev level map (n_new x n_old)
    GrB_Matrix *A_handle,     // output: coarsened adjacency (n_new x n_new)
    const GrB_Matrix C,       // input node->community map (n_old x n_old)
    GrB_Matrix A,             // input adjacency (n_old x n_old)
    bool free_A,              // if true, helper may free input A
    char *msg
)
{
    GrB_Vector c_list = NULL, x = NULL ;
    GrB_Matrix S_squished = NULL, S_new = NULL ;
    GrB_Matrix A_squished = NULL, A_new = NULL ;
    GrB_Descriptor desc = NULL ;

    uint64_t n, n_new ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;

    GrB_Type type_A ;
    GrB_Semiring plus_first, plus_second ;
    int code = 0 ;
    GRB_TRY (GrB_get (A, &code, GrB_EL_TYPE_CODE)) ;
    switch (code)
    {
        case GrB_BOOL_CODE:
        case GrB_INT8_CODE:
        case GrB_INT16_CODE:
        case GrB_INT32_CODE:
        case GrB_INT64_CODE:
        case GrB_UINT8_CODE:
        case GrB_UINT16_CODE:
        case GrB_UINT32_CODE:
        case GrB_UINT64_CODE:
            type_A = GrB_INT64 ;
            plus_first = GxB_PLUS_FIRST_INT64 ;
            plus_second = GxB_PLUS_SECOND_INT64 ;
            break ;

        case GrB_FP32_CODE:
        case GrB_FP64_CODE:
            type_A = GrB_FP64 ;
            plus_first = GxB_PLUS_FIRST_FP64 ;
            plus_second = GxB_PLUS_SECOND_FP64 ;
            break ;

        default:
            LG_ERROR_MSG ("LAGraph failed (file %s, line %d): unsupported "
                "adjacency matrix type for Louvain aggregation",
                __FILE__, __LINE__) ;
            LG_FREE_ALL ;
            return (GrB_NOT_IMPLEMENTED) ;
    }

    GRB_TRY (GrB_Descriptor_new (&desc)) ;
    GRB_TRY (GrB_set (desc, GxB_USE_INDICES, GxB_ROWINDEX_LIST)) ;
    GRB_TRY (GrB_set (desc, GxB_USE_INDICES, GxB_COLINDEX_LIST)) ;

    // Find active communities as column indices in C.
    GRB_TRY (GrB_Vector_new (&c_list, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new (&x, GrB_BOOL, n)) ;
    GRB_TRY (GrB_assign (x, NULL, NULL, (bool) 0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_vxm (c_list, NULL, NULL, GxB_ANY_PAIR_BOOL, x, C, NULL)) ;
    GRB_TRY (GrB_Vector_nvals (&n_new, c_list)) ;
    LG_ASSERT_MSG (n_new > 0, GrB_INVALID_VALUE, "No active communities") ;

    // S_squished: previous node -> active-community-id (contiguous columns)
    GRB_TRY (GrB_Matrix_new (&S_squished, GrB_BOOL, n, n_new)) ;
    GRB_TRY (GxB_Matrix_extract_Vector (
        S_squished, NULL, NULL, C, NULL, c_list, desc)) ;

    // A_squished: aggregate destination by community
    GRB_TRY (GrB_Matrix_new (&A_squished, type_A, n, n_new)) ;
    GRB_TRY (GrB_mxm (A_squished, NULL, NULL, plus_first, A, S_squished, NULL)) ;

    // S_new: active-community-id -> previous nodes
    GRB_TRY (GrB_Matrix_new (&S_new, GrB_BOOL, n_new, n)) ;
    GRB_TRY (GrB_transpose (S_new, NULL, NULL, S_squished, NULL)) ;

    // A_new: aggregate source by community
    GRB_TRY (GrB_Matrix_new (&A_new, type_A, n_new, n_new)) ;
    GRB_TRY (GrB_set (A_new, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT)) ;
    GRB_TRY (GrB_mxm (A_new, NULL, NULL, plus_second, S_new, A_squished, NULL)) ;

    if (free_A)
    {
        GRB_TRY (GrB_free (&A)) ;
    }

    *S_handle = S_new ; S_new = NULL ;
    *A_handle = A_new ; A_new = NULL ;

    LG_FREE_WORK ;
    return GrB_SUCCESS ;
}

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                 \
{                                                    \
    GrB_free (&C) ;                                  \
    GrB_free (&D) ;                                  \
    GrB_free (&it) ;                                 \
    GrB_free (&Ni) ;                                 \
    GrB_free (&desc) ;                               \
    GrB_free (&node_map) ;                           \
    GrB_free (&new_node_map) ;                       \
    GrB_free (&S_level) ;                            \
    GrB_free (&A_next) ;                             \
    GrB_free (&C_to_orig) ;                          \
    GrB_free (&active_comms) ;                       \
    GrB_free (&x_active) ;                           \
    LAGraph_Free ((void**) &degree, NULL) ;          \
    LAGraph_Free ((void**) &community_degree, NULL); \
    LAGraph_Free ((void**) &tuple_i, NULL) ;         \
    LAGraph_Free ((void**) &tuple_j, NULL) ;         \
    LAGraph_Free ((void**) &tuple_x, NULL) ;         \
    LAGraph_Free ((void**) &com_i, NULL) ;           \
    LAGraph_Free ((void**) &com_x, NULL) ;           \
    LAGraph_Free ((void**) &d_i, NULL) ;             \
    LAGraph_Free ((void**) &d_x, NULL) ;             \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                  \
{                                                    \
    if (free_A) GrB_free (&A) ;                      \
    LG_FREE_WORK ;                                   \
}

// compute clustering by running the Louvain algorithm against the graph's
// adjacency matrix A
GrB_Info LAGraph_louvain
(
    GrB_Vector *com,  // output communities
    LAGraph_Graph G,  // graph adjacency matrix
    int itermax,      // max number of modularity improvements sweeps per level
    int levelmax,     // max number of modularity improve and cluster condense
    float e,          // min change in modularity considered an improvement
    char *msg         // error message
)
{
    if (com == NULL || G == NULL || msg == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    GrB_Matrix C = NULL ;            // current level node->community map
    GrB_Vector D = NULL ;            // current level degree vector
    GxB_Iterator it = NULL ;
    GrB_Vector Ni = NULL ;           // node i's neighbors aggregated by comm
    GrB_Descriptor desc = NULL ;

    GrB_Matrix node_map = NULL ;     // current-level node -> original node map
    GrB_Matrix new_node_map = NULL ;
    GrB_Matrix S_level = NULL ;      // coarse->current node map from aggregate
    GrB_Matrix A_next = NULL ;       // coarsened adjacency
    GrB_Matrix C_to_orig = NULL ;
    GrB_Vector active_comms = NULL ;
    GrB_Vector x_active = NULL ;

    uint64_t *degree = NULL ;           // node degree for current level
    uint64_t *community_degree = NULL ; // sigma_tot per community id
    bool *tuple_x = NULL ;              // values for C_to_orig extraction
    GrB_Index *tuple_i = NULL, *tuple_j = NULL ;
    GrB_Index *com_i = NULL ;
    uint64_t *com_x = NULL ;
    GrB_Index *d_i = NULL ;
    uint64_t *d_x = NULL ;

    uint64_t *I = NULL ;      // readonly handle to C col-indices
    GrB_Index original_n = 0 ;
    GrB_Index nrows = 0, ncols = 0 ;
    GrB_Index final_nrows = 0 ;

    GrB_Matrix A = G->A ;     // current adjacency over levels
    bool free_A = false ;     // true once A becomes internally allocated

    // find out if graph is symmetric, compute cached values, and check loops
    LG_TRY (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;
    LG_TRY (LAGraph_Cached_OutDegree (G, msg)) ;
    LG_TRY (LAGraph_Cached_NSelfEdges (G, msg)) ;
    LG_ASSERT_MSG (G->nself_edges == 0, GrB_INVALID_VALUE,
        "G->nself_edges must be zero") ;

    GRB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    original_n = nrows ;

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT_MSG (nrows == ncols, LAGRAPH_INVALID_GRAPH,
        "adjacency matrix must be square") ;

    GrB_Type t ;
    GRB_TRY (GxB_Matrix_type (&t, A)) ;
    LG_ASSERT_MSG (t == GrB_BOOL, LAGRAPH_INVALID_GRAPH, "A must be boolean") ;

    bool improved = true ;

    for (int level = 0 ; level < levelmax ; level++)
    {
        GRB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
        GRB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
        LG_ASSERT_MSG (nrows == ncols, LAGRAPH_INVALID_GRAPH,
            "coarsened adjacency matrix must be square") ;
        final_nrows = nrows ;

        //----------------------------------------------------------------------
        // initialize degree for this level
        //----------------------------------------------------------------------

        GRB_TRY (GrB_Vector_new (&D, GrB_UINT64, nrows)) ;
        GRB_TRY (GrB_assign (D, NULL, NULL, (uint64_t) 0, GrB_ALL, nrows, NULL)) ;
        GRB_TRY (GrB_reduce (D, NULL, GrB_PLUS_UINT64, GrB_PLUS_MONOID_UINT64,
            A, NULL)) ;

        LG_TRY (LAGraph_Malloc ((void**) &degree, nrows, sizeof (uint64_t), msg)) ;
        memset (degree, 0, nrows * sizeof (uint64_t)) ;

        GrB_Index d_nvals = 0 ;
        GRB_TRY (GrB_Vector_nvals (&d_nvals, D)) ;
        LG_TRY (LAGraph_Malloc ((void**) &d_i, d_nvals, sizeof (GrB_Index), msg)) ;
        LG_TRY (LAGraph_Malloc ((void**) &d_x, d_nvals, sizeof (uint64_t), msg)) ;
        GRB_TRY (GrB_Vector_extractTuples_UINT64 (d_i, d_x, &d_nvals, D)) ;
        for (GrB_Index k = 0 ; k < d_nvals ; k++)
        {
            degree [d_i [k]] = d_x [k] ;
        }
        LG_TRY (LAGraph_Free ((void**) &d_i, msg)) ;
        LG_TRY (LAGraph_Free ((void**) &d_x, msg)) ;
        GRB_TRY (GrB_free (&D)) ;

        LG_TRY (LAGraph_Malloc ((void**) &community_degree, nrows,
            sizeof (uint64_t), msg)) ;
        memcpy (community_degree, degree, nrows * sizeof (uint64_t)) ;

        //----------------------------------------------------------------------
        // initialize communities C and readonly community id handle I
        //----------------------------------------------------------------------

        GRB_TRY (_initCommunities (&C, &I, nrows)) ;
        LG_ASSERT (I != NULL, GrB_NULL_POINTER) ;

        GRB_TRY (GrB_Vector_new (&Ni, GrB_UINT64, ncols)) ;
        GRB_TRY (GxB_Iterator_new (&it)) ;
        GRB_TRY (GrB_Descriptor_new (&desc)) ;
        GRB_TRY (GrB_set (desc, GrB_STRUCTURE, GrB_MASK_FIELD)) ;
        GRB_TRY (GrB_set (desc, GxB_USE_INDICES, GxB_ROWINDEX_LIST)) ;

        //----------------------------------------------------------------------
        // phase 1: local moving
        //----------------------------------------------------------------------

        double modularity_gain = 0 ;
        improved = true ;

        double M2 = 0.0 ;
        for (GrB_Index i = 0 ; i < nrows ; i++)
        {
            M2 += (double) degree [i] ;
        }
        LG_ASSERT_MSG (M2 > 0.0, GrB_INVALID_VALUE,
            "sum of degrees must be positive") ;

        for (int iter = 0 ; iter < itermax && improved ; iter++)
        {
            modularity_gain = 0 ;

            for (GrB_Index i = 0 ; i < nrows ; i++)
            {
                uint64_t ic = I [i] ;

                GRB_TRY (GrB_Col_extract (
                    Ni, NULL, NULL, A, GrB_ALL, ncols, i, GrB_DESC_T0)) ;

                // Aggregate node i's edge-weight to each candidate community.
                GRB_TRY (GrB_vxm (
                    Ni, NULL, NULL, GxB_PLUS_FIRST_UINT64, Ni, C, NULL)) ;

                uint64_t i_degree = degree [i] ;

                uint64_t kin = 0 ;
                GrB_Info kin_info = GrB_Vector_extractElement (&kin, Ni, ic) ;
                if (kin_info != GrB_SUCCESS && kin_info != GrB_NO_VALUE)
                {
                    GRB_TRY (kin_info) ;
                }
                if (kin_info == GrB_NO_VALUE) kin = 0 ;

                uint64_t sigma_tot = community_degree [ic] ;
                sigma_tot -= i_degree ;

                double base_gain = gain (i_degree, kin, sigma_tot, M2) ;
                double max_modularity = base_gain ;
                uint64_t best_community = ic ;

                GRB_TRY (GxB_Vector_Iterator_attach (it, Ni, NULL)) ;
                GrB_Info info = GxB_Vector_Iterator_seek (it, 0) ;
                while (info != GxB_EXHAUSTED)
                {
                    GrB_Index jc = GxB_Vector_Iterator_getIndex (it) ;
                    uint64_t kjn = GxB_Iterator_get_UINT64 (it) ;
                    info = GxB_Vector_Iterator_next (it) ;

                    if (jc == ic)
                    {
                        continue ;
                    }

                    sigma_tot = community_degree [jc] ;
                    double g = gain (i_degree, kjn, sigma_tot, M2) ;
                    if (g > max_modularity)
                    {
                        best_community = jc ;
                        max_modularity = g ;
                    }
                }

                if (best_community != ic)
                {
                    I [i] = best_community ;
                    community_degree [ic] -= i_degree ;
                    community_degree [best_community] += i_degree ;
                    modularity_gain += max_modularity - base_gain ;
                }
            }

            improved = (modularity_gain > (double) e) ;
        }

        //----------------------------------------------------------------------
        // decide whether to stop at this level before aggregation
        //----------------------------------------------------------------------

        if (!improved || (level + 1) >= levelmax || nrows <= 1)
        {
            break ;
        }

        // Count active communities at this level.
        GRB_TRY (GrB_Vector_new (&active_comms, GrB_BOOL, nrows)) ;
        GRB_TRY (GrB_Vector_new (&x_active, GrB_BOOL, nrows)) ;
        GRB_TRY (GrB_assign (x_active, NULL, NULL, (bool) 0, GrB_ALL, nrows, NULL)) ;
        GRB_TRY (GrB_vxm (active_comms, NULL, NULL, GxB_ANY_PAIR_BOOL,
            x_active, C, NULL)) ;
        GrB_Index n_communities = 0 ;
        GRB_TRY (GrB_Vector_nvals (&n_communities, active_comms)) ;
        GRB_TRY (GrB_free (&active_comms)) ;
        GRB_TRY (GrB_free (&x_active)) ;

        // No further coarsening is possible.
        if (n_communities == nrows || n_communities <= 1)
        {
            break ;
        }

        //----------------------------------------------------------------------
        // phase 2: aggregate
        //----------------------------------------------------------------------

        GRB_TRY (LG_Louvain_aggregate (&S_level, &A_next, C, A, free_A, msg)) ;

        // Compose level maps so node_map always maps current nodes -> original.
        if (node_map == NULL)
        {
            node_map = S_level ;
            S_level = NULL ;
        }
        else
        {
            GRB_TRY (GrB_Matrix_new (&new_node_map, GrB_BOOL,
                n_communities, original_n)) ;
            GRB_TRY (GrB_mxm (new_node_map, NULL, NULL, GxB_ANY_PAIR_BOOL,
                S_level, node_map, NULL)) ;
            GRB_TRY (GrB_free (&node_map)) ;
            GRB_TRY (GrB_free (&S_level)) ;
            node_map = new_node_map ;
            new_node_map = NULL ;
        }

        // Prepare next level.
        GRB_TRY (GrB_free (&C)) ;
        I = NULL ;
        A = A_next ;
        A_next = NULL ;
        free_A = true ;

        GRB_TRY (GrB_free (&it)) ;
        GRB_TRY (GrB_free (&Ni)) ;
        GRB_TRY (GrB_free (&desc)) ;
        LG_TRY (LAGraph_Free ((void**) &degree, msg)) ;
        LG_TRY (LAGraph_Free ((void**) &community_degree, msg)) ;
    }

    //----------------------------------------------------------------------
    // set output
    //----------------------------------------------------------------------

    LG_ASSERT (C != NULL, GrB_NULL_POINTER) ;
    GRB_TRY (GrB_Vector_new (com, GrB_UINT64, original_n)) ;

    if (node_map == NULL)
    {
        LG_ASSERT_MSG (final_nrows == original_n, GrB_INVALID_VALUE,
            "unexpected shape mismatch without node aggregation") ;
        GRB_TRY (GxB_Vector_load (*com, (void**) (&I), GrB_UINT64, original_n,
            sizeof (uint64_t) * original_n, GrB_DEFAULT, NULL)) ;
    }
    else
    {
        // C_to_orig = C' * node_map: community rows, original-node columns
        GRB_TRY (GrB_Matrix_new (&C_to_orig, GrB_BOOL, final_nrows, original_n)) ;
        GRB_TRY (GrB_mxm (C_to_orig, NULL, NULL, GxB_ANY_PAIR_BOOL,
            C, node_map, GrB_DESC_T0)) ;

        GrB_Index cnvals = 0 ;
        GRB_TRY (GrB_Matrix_nvals (&cnvals, C_to_orig)) ;
        LG_TRY (LAGraph_Malloc ((void**) &tuple_i, cnvals, sizeof (GrB_Index), msg)) ;
        LG_TRY (LAGraph_Malloc ((void**) &tuple_j, cnvals, sizeof (GrB_Index), msg)) ;
        LG_TRY (LAGraph_Malloc ((void**) &tuple_x, cnvals, sizeof (bool), msg)) ;
        LG_TRY (LAGraph_Malloc ((void**) &com_i, cnvals, sizeof (GrB_Index), msg)) ;
        LG_TRY (LAGraph_Malloc ((void**) &com_x, cnvals, sizeof (uint64_t), msg)) ;

        GRB_TRY (GrB_Matrix_extractTuples_BOOL (
            tuple_i, tuple_j, tuple_x, &cnvals, C_to_orig)) ;

        for (GrB_Index k = 0 ; k < cnvals ; k++)
        {
            com_i [k] = tuple_j [k] ;
            com_x [k] = (uint64_t) tuple_i [k] ;
        }

        GRB_TRY (GrB_Vector_build_UINT64 (
            *com, com_i, com_x, cnvals, GrB_SECOND_UINT64)) ;
    }

    LG_FREE_ALL ;
    return GrB_SUCCESS ;
}
