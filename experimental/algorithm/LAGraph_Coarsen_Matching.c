//------------------------------------------------------------------------------
// LAGraph_Coarsen_Matching: Coarsen an undirected graph using an edge matching
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Vidith Madhu, Texas A&M University

//------------------------------------------------------------------------------

// TODO: ready for src? need a vanilla non-GxB, and incidence graphs.

/*
This method is used to coarsen an undirected graph. The coarsening is based on a maximal matching,
which is handled by LAGraph_MaximalMatching. 

The coarsening step involves a reduction from a graph G to G', where we use a bijection f from
nodes in G to nodes in G'. We can consider f(u) to be the parent of node u. 
For each edge (u, v) in G, we add an edge (f(u), f(v)) to G' iff f(u) != f(v). In our case,
this bijection is given by the maximal matching, where for every matched edge, one of the endpoints of the edge is the
parent (representative) of both endpoints, and any node not part of a matched edge is its own parent.

This method performs a single coarsening step on the input graph.

The inputs to this algorithm are as follows in order:
1. an LAGraph_Graph containing the target graph to coarsen
2. the type of matching to perform (random, heavy, or light)
3. whether to retain the size of the graph when coarsening. If 1, then nodes that are eliminated by a coarsening step
are turned into singletons. If 0, the size of the graph is changed and nodes are explicitly relabeled.
4. whether edges that are combined during a coarsening step should have their edge weights summed (for an unweighted graph, this
counts the number of combined edges). If this option is false, then only the pattern of combined edges is retained.
6. Random seed used for maximal matching (same for every coarsening step)
7. msg for LAGraph error reporting

There are 4 outputs from the function:
1. A GrB_Matrix of the coarsened graph (if the input adjacency matrix is of type GrB_BOOL or GrB_UINT{8|16|32} or GrB_INT*, it will
have type GrB_INT64. If it is of type GrB_FP32, it will have type GrB_FP64. Else, it will have the same type as the input matrix.

2. A full GrB_Vector (parent_result) of length n where if parent_result[u] = v,
then node u has parent v. This parent mapping is derived from a maximal matching of the graph
and is used for the coarsening step (meaning node u collapses into node v).

3. A GrB_Vector (newlabels_result) of length n where if newlabels_result[u] = v,
then node u in G is relabeled as node v in G', where G' is the coarsened graph. In addition, newlabels_result[u] exists iff node u survives the i-th coarsening step.
If preserve_mapping = 1, then this result is returned as NULL since no relabeling occurs. This result is used to interpret the contents of parent_result, 
since the labels of nodes may change in an arbitrary fashion for each coarsening step.

4. A full GrB_Vector (inv_newlabels_result) of length n' (the number of vertices in the coarsened graph) where if inv_newlabels_result[u] = v,
then node u in G' had an original label as node v in G, where G' is the coarsened graph. In other words, this is simply the inverse of result (3).
If preserve_mapping = 1, this is returned as NULL.

NOTE: Results (2), (3), and (4) are only computed/returned if they are not passed as a NULL pointer on input.
Passing a NULL pointer denotes that the user does not want that particular result, and no return value is given.
This means that a NULL pointer may safely be passed on input.

This method requires O(n + e) space for an undirected graph with e edges and n nodes
*/

#include "LG_internal.h"
#include "LAGraphX.h"

// #define dbg
// #define burble

#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define LG_FREE_WORK                                        \
{                                                           \
    GrB_free(&parent_cpy) ;                                 \
    GrB_free(&one) ;                                        \
    GrB_free(&VALUEEQ_ROWINDEX_UINT64) ;                    \
    LAGraph_Free ((void**) &original_indices, NULL) ;       \
    LAGraph_Free ((void**) &original_values, NULL) ;        \
    LAGraph_Free ((void**) &preserved_values, NULL) ;       \
    LAGraph_Free ((void**) &S_rows, NULL) ;                 \
    LAGraph_Free ((void**) &S_cols, NULL) ;                 \
}                                                           \

#define LG_FREE_ALL                                 \
{                                                   \
    LG_FREE_WORK ;                                  \
    GrB_free(&S) ;                                  \
}                                                   \

#define F_INDEX_UNARY(f)  ((void (*)(void *, const void *, GrB_Index, GrB_Index, const void *)) f)

#if LAGRAPH_SUITESPARSE

void LG_CM_valueeq_index_func (bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const void *y) {
    (*z) = ((*x) == i) ;
}

static int LAGraph_Parent_to_S
(   
    // input/outputs:
    GrB_Matrix *result,             // resulting S matrix
    GrB_Vector *newlabels,          // The contents of some newlabels_result[i], where newlabels_result is as described at the top of the file.
                                    // if NULL on input, will not return anything. If preserve_mapping = 1, the return value is a NULL GrB_Vector.

    GrB_Vector *inv_newlabels,      // The contents of some inv_newlabels_result[i], where inv_newlabels_result is as described at the top of the file.
                                    // if NULL on input, will not return anything. If preserve_mapping = 1, the return value is a NULL GrB_Vector.

    // strictly inputs:
    GrB_Vector parent,              // dense integer vector of size n. parent[i] -> representative of node i
    int preserve_mapping,           // whether to preserve the original namespace of nodes, or to compress it down
    GrB_Type S_type,                // type of constructed S matrix
    char *msg
)
{
    GrB_Vector parent_cpy = NULL ;  // used so we don't modify the input parent vector. Also useful to have for computing newlabels
    GrB_Matrix S = NULL ;           // stores the resulting S matrix

    // ------------------------ BUILDING S MATRIX --------------------------
    // used to unpack parent_cpy and build S
    GrB_Index *S_cols = NULL, *S_rows = NULL, S_cols_size, S_rows_size ;
    GrB_Scalar one = NULL ;
    GrB_Index nvals ;
    // ---------------------------------------------------------------------
    

    // ----------------------- NODE COMPRESSION -----------------------------
    // No need to free this since it gets packed back
    uint64_t *ramp = NULL ;
    // used to unpack preserved nodes
    // note: No need to free the array since it is packed back
    GrB_Index *preserved_indices = NULL, preserved_indices_size, preserved_values_size ;
    // need to free this (not packed back)
    uint64_t *preserved_values = NULL ;
    // used to extractTuples for original parent
    GrB_Index *original_indices = NULL ;
    uint64_t *original_values = NULL ;

    GrB_Index num_preserved ;
    bool is_jumbled ;

    GrB_IndexUnaryOp VALUEEQ_ROWINDEX_UINT64 = NULL ;
    // ---------------------------------------------------------------------

    // number of nodes in original graph
    GrB_Index n ;

    GRB_TRY (GrB_Vector_nvals (&n, parent)) ;

    if (preserve_mapping) {
        // parent_cpy will be the same as parent
        GRB_TRY (GrB_Vector_dup (&parent_cpy, parent)) ;
    } else {
        // we want an empty vector for the compression step (see below code)
        GRB_TRY (GrB_Vector_new (&parent_cpy, GrB_UINT64, n)) ;
    }

    if (!preserve_mapping) {
        /*
        new code:
            - identify preserved nodes (GrB_select to parent_cpy)
            - unpack to get indices of preserved nodes
            - build ramp vector
            - pack back into parent_cpy with ramp as values, preserved node indices as indices (performs compression)
            - GrB_extract into parent_cpy from parent_cpy with row indices as values from original parent
                - This fills in the new parents for discarded nodes
        */
        
        GRB_TRY (GrB_IndexUnaryOp_new (&VALUEEQ_ROWINDEX_UINT64, F_INDEX_UNARY(LG_CM_valueeq_index_func), GrB_BOOL, GrB_UINT64, GrB_UINT64)) ;

        // identify preserved nodes
        GRB_TRY (GrB_select (parent_cpy, NULL, NULL, VALUEEQ_ROWINDEX_UINT64, parent, 0, NULL)) ;
        GRB_TRY (GrB_free (&VALUEEQ_ROWINDEX_UINT64)) ;
        
        // get indices of preserved nodes
        GRB_TRY (GxB_Vector_unpack_CSC (
            parent_cpy, 
            &preserved_indices, 
            (void**) &preserved_values, 
            &preserved_indices_size, 
            &preserved_values_size, 
            NULL, 
            &num_preserved,
            &is_jumbled,
            NULL
        )) ;

        LG_TRY (LAGraph_Free ((void**)(&preserved_values), msg)) ;
        
        // build ramp vector
        LG_TRY (LAGraph_Malloc ((void**) &ramp, num_preserved, sizeof(uint64_t), msg)) ;

        int64_t i;
        #pragma omp parallel for
        for (i = 0 ; i < num_preserved ; i++) {
            ramp [i] = i ;
        }

        if (inv_newlabels != NULL) {
            GRB_TRY (GrB_Vector_new (inv_newlabels, GrB_UINT64, num_preserved)) ;
            GRB_TRY (GrB_Vector_build (*inv_newlabels, ramp, preserved_indices, num_preserved, NULL)) ;
        }

        // pack back into parent_cpy (parent_cpy will now store the new labels of preserved nodes)
        GRB_TRY (GxB_Vector_pack_CSC (
            parent_cpy,
            &preserved_indices,
            (void**) &ramp,
            num_preserved * sizeof(GrB_Index),
            num_preserved * sizeof(uint64_t),
            false,
            num_preserved,
            is_jumbled,
            NULL
        )) ;

        if (newlabels != NULL) {
            GRB_TRY (GrB_Vector_dup (newlabels, parent_cpy)) ;
        }

        LG_TRY (LAGraph_Malloc ((void**) &original_indices, n, sizeof(GrB_Index), msg)) ;
        LG_TRY (LAGraph_Malloc ((void**) &original_values, n, sizeof(uint64_t), msg)) ;

        GRB_TRY (GrB_Vector_extractTuples (original_indices, original_values, &n, parent)) ;

        // fill in entries for discarded nodes
        GRB_TRY (GrB_Vector_extract (parent_cpy, NULL, NULL, parent_cpy, original_values, n, NULL)) ;

        GRB_TRY (GrB_Matrix_new (&S, S_type, num_preserved, n)) ;

    } else { 
        // result dim: n by n
        GRB_TRY (GrB_Matrix_new (&S, S_type, n, n)) ;
        // newlabels is the identity map, signified by null return values
        if (newlabels != NULL) {
            (*newlabels) = NULL ;
        }
        if (inv_newlabels != NULL) {
            (*inv_newlabels) = NULL ;
        }
    }

    GRB_TRY (GxB_Vector_unpack_CSC (parent_cpy, &S_cols, (void**) &S_rows, &S_cols_size, &S_rows_size, NULL, &nvals, NULL, NULL)) ;
    GRB_TRY (GrB_Scalar_new (&one, S_type)) ;
    GRB_TRY (GrB_Scalar_setElement (one, 1)) ;

    GRB_TRY (GxB_Matrix_build_Scalar (S, S_rows, S_cols, one, nvals)) ;    

    (*result) = S ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
#endif


#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define LG_FREE_WORK                                \
{                                                   \
    GrB_free(&E) ;                                  \
    GrB_free(&S) ;                                  \
    GrB_free(&matched_edges) ;                      \
    GrB_free(&E_t) ;                                \
    GrB_free(&S_t) ;                                \
    GrB_free(&edge_parent) ;                        \
    GrB_free(&node_parent) ;                        \
    GrB_free(&full) ;                               \
    LAGraph_Delete(&G_cpy, msg) ;                   \
}

#define LG_FREE_ALL                                 \
{                                                   \
    LG_FREE_WORK ;                                  \
}

#ifdef burble                                      
    #define CHKPT(msg){ printf("*** [CHKPT] *** %s\n", msg) ; }                                                            
#else                                                   
    #define CHKPT(msg){}
#endif

#define OPTIMIZE_PUSH_PULL

int LAGraph_Coarsen_Matching
(
    // outputs:
    GrB_Matrix *coarsened,                 // coarsened adjacency
    GrB_Vector *parent_result,             // refer to comments at top of file
    GrB_Vector *newlabel_result,           // refer to comments at top of file
    GrB_Vector *inv_newlabel_result,       // refer to comments at top of file

    // inputs:
    LAGraph_Graph G,                       // input graph
    LAGraph_Matching_kind matching_type,   // how to perform the coarsening
    bool preserve_mapping,                 // preserve original namespace of nodes
    bool combine_weights,                  // whether to sum edge weights or just keep the pattern
    uint64_t seed,                         // seed used for matching
    char *msg
)
{
#if LAGRAPH_SUITESPARSE

    LG_CLEAR_MSG ;

    LAGraph_Graph G_cpy = NULL ;            // used for the IncidenceMatrix function
    GrB_Matrix C = NULL ;                   // resulting adjacency matrix (used for output)
    GrB_Matrix E = NULL ;                   // incidence matrix
    GrB_Matrix E_t = NULL ;                 // transpose of incidence
    GrB_Matrix S = NULL ;                   // S matrix (S[i][j] -> node j maps to node i in coarsened graph)
    GrB_Matrix S_t = NULL ;                 // transpose of S
    GrB_Vector matched_edges = NULL ;       // result of maximal matching
    GrB_Vector edge_parent = NULL ;         // points to parent (representative) node for each edge
    GrB_Vector node_parent = NULL ;         // points to parent (representative) node for each node
    GrB_Vector full = NULL ;                // full vector

    GrB_Index nvals, nrows ;
    GrB_Type A_type ;

    // check properties (no self-loops, undirected
    LG_ASSERT_MSG (G->nself_edges == 0, LAGRAPH_NO_SELF_EDGES_ALLOWED, "G->nself_edges must be zero") ;

    LG_ASSERT (coarsened != NULL, GrB_NULL_POINTER) ;

    //------------------------------------------------------------------------------
    // check input graph, build local adjacency matrix to use for coarsening
    //------------------------------------------------------------------------------

    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED)
    {
        char typename[LAGRAPH_MAX_NAME_LEN] ;
        GrB_Type type ;
        LG_TRY (LAGraph_Matrix_TypeName (typename, G->A, msg)) ;
        LG_TRY (LAGraph_TypeFromName (&type, typename, msg)) ;

        if ((type == GrB_FP64) || (type == GrB_INT64 || type == GrB_UINT64)) {
            // output will keep the same type as input
            GRB_TRY (GrB_Matrix_dup (&C, G->A)) ;
            A_type = type ;
        } else {
            // output will become int64/fp64
            // reasoning: want to prevent overflow from combining edges and accomodate negative edge weights
            #ifdef dbg
                printf("Rebuilding with GrB_INT64/FP64, orig type was %s\n", typename);
            #endif

            bool is_float = (type == GrB_FP32) ;
            GRB_TRY (GrB_Matrix_nrows (&nrows, G->A)) ;
            A_type = (is_float ? GrB_FP64 : GrB_INT64) ;
            GRB_TRY (GrB_Matrix_new (&C, A_type, nrows, nrows)) ;
            GRB_TRY (GrB_assign (C, NULL, NULL, G->A, GrB_ALL, nrows, GrB_ALL, nrows, NULL)) ;
        }
    }
    else
    {
        // G is not undirected
        LG_ASSERT_MSG (false, LAGRAPH_INVALID_GRAPH, "G must be undirected") ;
    }
    CHKPT("Done with building C");

    // make new LAGraph_Graph to use for LAGraph_IncidenceMatrix and for useful functions (delete self-edges)
    LG_TRY (LAGraph_New (&G_cpy, &C, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    ASSERT (C == NULL) ;
    LG_TRY (LAGraph_Cached_NSelfEdges (G_cpy, msg)) ;

    GrB_Index num_nodes ;
    GrB_Index num_edges ;

    GRB_TRY (GrB_Matrix_nrows (&num_nodes, G_cpy->A)) ;
    GRB_TRY (GrB_Matrix_nvals (&num_edges, G_cpy->A)) ;
    num_edges /= 2 ; // since undirected

    CHKPT("Done building G_cpy");

    GRB_TRY (GrB_Matrix_new (&E_t, A_type, num_edges, num_nodes)) ;
    GRB_TRY (GrB_Vector_new (&edge_parent, GrB_UINT64, num_edges)) ;

    GRB_TRY (GrB_Vector_new (&node_parent, GrB_UINT64, num_nodes)) ;

    GRB_TRY (GrB_Vector_new (&full, GrB_BOOL, num_nodes)) ;
    GRB_TRY (GrB_assign (full, NULL, NULL, true, GrB_ALL, num_nodes, NULL)) ;

    // for push/pull optimization
    double sparsity_thresh = 
    #ifdef OPTIMIZE_PUSH_PULL
        0.04 ;
    #else
        1.0 ;
    #endif

    GrB_Index curr_level = 0 ;
    CHKPT("Starting coarsening step");

    //------------------------------------------------------------------------------
    // coarsening step
    //------------------------------------------------------------------------------

    // get incidence matrix
    LG_TRY (LAGraph_Incidence_Matrix (&E, G_cpy, msg)) ;
    CHKPT("Done with LAGraph_IncidenceMatrix");

    GRB_TRY (GrB_transpose (E_t, NULL, NULL, E, NULL)) ;
    // set to row major incase Incidence_Matrix gave col_major
    GRB_TRY (GrB_set(E, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT)) ;
    CHKPT("Starting maximal matching");
    // run maximal matching
    LG_TRY (LAGraph_MaximalMatching (&matched_edges, E, E_t, matching_type, seed, msg)) ;
    CHKPT("Done with maximal matching");

    // make edge_parent
    // want to do E_t * full and get the first entry for each edge (mask output with matched_edges)
    GRB_TRY (GrB_mxv (edge_parent, matched_edges, NULL, GxB_MIN_SECONDI_INT64, E_t, full, GrB_DESC_RS)) ;

    #ifdef dbg
        printf("Printing edge_parent for level (%lld)\n", curr_level) ;
        LG_TRY (LAGraph_Vector_Print (edge_parent, LAGraph_COMPLETE, stdout, msg)) ;
        printf("Printing E for level (%lld)\n", curr_level) ;
        LG_TRY (LAGraph_Matrix_Print (E, LAGraph_COMPLETE, stdout, msg)) ;
    #endif
    // now, we have edge_parent (each edge points to its parent node)
    // can do E * edge_parent with min_second to get node_parent
    GrB_Index num_matched ;
    GRB_TRY (GrB_Vector_nvals (&num_matched, edge_parent)) ;
    
    if (num_matched > sparsity_thresh * num_edges) {
        GRB_TRY (LG_SET_FORMAT_HINT (edge_parent, LG_BITMAP)) ;
        GRB_TRY (GrB_mxv (node_parent, NULL, NULL, GrB_MIN_SECOND_SEMIRING_UINT64, E, edge_parent, NULL)) ;
    } else {
        GRB_TRY (LG_SET_FORMAT_HINT (edge_parent, LG_SPARSE)) ;
        GRB_TRY (GrB_vxm (node_parent, NULL, NULL, GrB_MIN_FIRST_SEMIRING_UINT64, edge_parent, E_t, NULL)) ;
    }

    // populate non-existent entries in node_parent with their index
    // handles nodes that are not engaged in a matching
    GrB_apply (node_parent, node_parent, NULL, GrB_ROWINDEX_INT64, full, (uint64_t) 0, GrB_DESC_SC) ;

    if (parent_result != NULL) {
        // record a deep copy of the current node_parent for the output parent vector
        GRB_TRY (GrB_Vector_dup (parent_result, node_parent)) ;
    }

    #ifdef dbg
        printf("Printing node_parent for level (%lld)\n", curr_level) ;
        LG_TRY (LAGraph_Vector_Print (node_parent, LAGraph_COMPLETE, stdout, msg)) ;
        printf("Printing matched edges for level (%lld)\n", curr_level) ;
        LG_TRY (LAGraph_Vector_Print (matched_edges, LAGraph_COMPLETE, stdout, msg)) ;
    #endif

    // build the S matrix
    LG_TRY (LAGraph_Parent_to_S (
        &S,
        newlabel_result,
        inv_newlabel_result,
        node_parent,
        preserve_mapping, 
        A_type,
        msg
    )) ;

    #ifdef dbg
        printf("Printing S for level (%lld)\n", curr_level) ;
        LG_TRY (LAGraph_Matrix_Print (S, LAGraph_COMPLETE, stdout, msg)) ;
        printf("Printing C for level (%lld)\n", curr_level) ;
        LG_TRY (LAGraph_Matrix_Print (G_cpy->A, LAGraph_COMPLETE, stdout, msg)) ;
    #endif
    
    GrB_Index S_nrows, S_ncols ;

    // create S_t now that we know the dimensions of S
    GRB_TRY (GrB_Matrix_nrows (&S_nrows, S)) ;
    GRB_TRY (GrB_Matrix_ncols (&S_ncols, S)) ;
    
    GRB_TRY (GrB_Matrix_new (&S_t, A_type, S_ncols, S_nrows)) ;

    GRB_TRY (GrB_transpose (S_t, NULL, NULL, S, NULL)) ;

    // choose right semiring based on combine_weights and type of adjacency matrix
    GrB_Semiring combine_semiring = (A_type == GrB_FP64) ? GrB_PLUS_TIMES_SEMIRING_FP64 : GrB_PLUS_TIMES_SEMIRING_INT64 ;
    GrB_Semiring semiring = combine_weights ? combine_semiring : LAGraph_any_one_bool ;
    
    GRB_TRY (GrB_mxm (S, NULL, NULL, semiring, S, G_cpy->A, NULL)) ;

    #ifdef dbg
        printf("Printing S * C for level (%lld)\n", curr_level) ;
        LG_TRY (LAGraph_Matrix_Print (S, LAGraph_COMPLETE, stdout, msg)) ;
    #endif

    if (!preserve_mapping) {
        // re-instantiate adjacency matrix with new dimensions
        GRB_TRY (GrB_free (&G_cpy->A)) ;
        GRB_TRY (GrB_Matrix_new (&(G_cpy->A), A_type, S_nrows, S_nrows)) ;
    }
    GRB_TRY (GrB_mxm (G_cpy->A, NULL, NULL, semiring, S, S_t, NULL)) ;

    // make nself_edges unknown for delete self edges to work
    G_cpy->nself_edges = LAGRAPH_UNKNOWN ;
    // parent nodes for matched edges will form self-edges; need to delete
    LG_TRY (LAGraph_DeleteSelfEdges (G_cpy, msg)) ;

    //------------------------------------------------------------------------------
    // coarsening step done
    //------------------------------------------------------------------------------

    (*coarsened) = G_cpy->A ;
    G_cpy->A = NULL ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}

