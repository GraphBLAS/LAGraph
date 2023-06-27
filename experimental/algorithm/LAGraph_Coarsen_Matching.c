//------------------------------------------------------------------------------
// LAGraph_Coarsen_Matching: Coarsen an undirected graph using an edge matching
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Vidith Madhu, Texas A&M University

//------------------------------------------------------------------------------

/*
This method is used to coarsen an undirected graph. The coarsening is based on a maximal matching,
which is handled by LAGraph_MaximalMatching.
The inputs to this algorithm are as follows in order:
1. an LAGraph_Graph containing the target graph to coarsen
2. the type of matching to perform (random, heavy, or light)
3. whether to retain the size of the graph when coarsening. If 1, then nodes that are eliminated by a coarsening step
are turned into singletons. If 0, the size of the graph is changed and nodes are explicitly relabeled.
4. whether edges that are combined during a coarsening step should have their edge weights summed (for an unweighted graph, this
counts the number of combined edges)
5. How many coarsening steps to perform
6. Random seed used for maximal matching
7. msg for LAGraph error reporting

There are 3 outputs from the function:
1. A GrB_Matrix of the coarsened graph (if the input adjacency matrix is of type GrB_BOOL or GrB_UINT* or GrB_INT*, it will
have type GrB_INT64. Else, it will have the same type as the input matrix).
2. A list of GrB_Vectors (parent_result) of length nlevels, where if parent_result[i][u] = v,
then the parent of node u in G_{i} is node v in G_{i}, where G_0 is the initial graph. Note that this means 
the length of parent_result[i] is the number of nodes in G_{i}. If preserve_mapping = 1, then there is no need 
for such a result, and a NULL pointer is returned.
3. A list of GrB_Vectors (newlabels_result) of length nlevels, where if newlabels_result[i][u] = v,
then node u in G_{i} is relabeled as node v in G_{i + 1}. Again, the length of newlabels_result[i] is the number
of nodes in G_{i}; if preserve_mapping = 1, then this result is returned as NULL. This result is used to interpret
the contents of parent_result, since the naming of nodes may change in an arbitrary fashion for successive coarsenings.

More specifically, the coarsening step involves a reduction from a graph G to G', where we use a bijection f from
nodes in G to nodes in G'. We can consider f(u) to be the parent of node u. 
For each edge (u, v) in G, we add an edge (f(u), f(v)) to G' iff f(u) != f(v). In our case,
this bijection is given by the maximal matching, where for every matched edge, one of the endpoints of the edge is the
parent (representative) of both endpoints, and any node not part of a matched edge is its own parent.  
*/

#include "LG_internal.h"
#include "LAGraphX.h"

// #define dbg

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

void valueeq_index_func (bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const void *y) {
    (*z) = ((*x) == i) ;
}

static int LAGraph_Parent_to_S
(   
    // input/outputs:
    GrB_Matrix *result,             // resulting S matrix
    GrB_Vector *newlabels,          // if not NULL on input, and preserve_mapping = false, will return the new labels of nodes
                                    // more specifically, on output, mapping[i] exists iff node i survives in the coarsened graph, and has
                                    // value equal to the new label of node i.
                                    // if preserve_mapping = true, is returned as NULL on output

    GrB_Vector *inv_newlabels,      // same as above except the inverse (mapping of new nodes to old nodes)
                                    // more specifically, on output, mapping[i] is the node in the uncoarsened graph that maps to node i in the coarsened graph
                                    // again, if preserve_mapping = true, is NULL on output

    // strictly inputs:
    GrB_Vector parent,              // dense integer vector of size n. parent[i] -> representative of node i
    int preserve_mapping,           // whether to preserve the original namespace of nodes, or to compress it down
    char *msg
)
{
    GrB_Vector parent_cpy = NULL ;          // don't modify input parent. Also useful to have for compressing node labels.
    GrB_Matrix S = NULL ;

    // ------------------------ BUILDING S MATRIX --------------------------
    // used to unpack parent_cpy and build S
    GrB_Index *S_cols = NULL, *S_rows = NULL, S_cols_size, S_rows_size ;
    GrB_Scalar one = NULL ;
    GrB_Index nvals ;
    // ---------------------------------------------------------------------
    

    // ----------------------- NODE COMPRESSION -----------------------------
    // [0...(n - 1)]
    // No need to free this since it gets packed
    uint64_t *ramp = NULL ;
    // used to unpack preserved nodes
    // note: No need to free this array since they are packed back
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
            - unpack to get values of preserved
            - build ramp vector (full vector, GrB_apply w/ row index)
            - pack back into parent_cpy with ramp as values, preserved as indices
            - GrB_extract into parent_cpy from parent_cpy with row indices as values from original parent
                - This fills in the new parents for discarded nodes
        */
        
        GRB_TRY (GrB_IndexUnaryOp_new (&VALUEEQ_ROWINDEX_UINT64, F_INDEX_UNARY(valueeq_index_func), GrB_BOOL, GrB_UINT64, GrB_UINT64)) ;

        // identify preserved nodes
        GRB_TRY (GrB_select (parent_cpy, NULL, NULL, VALUEEQ_ROWINDEX_UINT64, parent, 0, NULL)) ;
        
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
        
        // build ramp vector
        // TODO: parallelize this
        LG_TRY (LAGraph_Malloc ((void**) &ramp, num_preserved, sizeof(uint64_t), msg)) ;
        for (GrB_Index i = 0 ; i < num_preserved ; i++) {
            ramp [i] = i ;
        }

        if (inv_newlabels != NULL) {
            GRB_TRY (GrB_Vector_new (inv_newlabels, GrB_UINT64, num_preserved)) ;
            GRB_TRY (GrB_Vector_build (*inv_newlabels, ramp, preserved_indices, num_preserved, GrB_SECOND_INT64)) ;
        }

        // pack back into parent_cpy
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
            // alternate method:
            // (*newlabels) = parent_cpy
            // and free parent_cpy in LG_FREE_ALL instead of LG_FREE_WORK 
            GRB_TRY (GrB_Vector_dup (newlabels, parent_cpy)) ;
        }

        LG_TRY (LAGraph_Malloc ((void**) &original_indices, n, sizeof(GrB_Index), msg)) ;
        LG_TRY (LAGraph_Malloc ((void**) &original_values, n, sizeof(uint64_t), msg)) ;

        // is extractTuples needed? Since this is a static function (and won't be exposed to users), 
        // is it OK to just unpack and rip the contents of the input out?
        GRB_TRY (GrB_Vector_extractTuples (original_indices, original_values, &n, parent)) ;

        // fill in entries for discarded nodes
        GRB_TRY (GrB_Vector_extract (parent_cpy, NULL, NULL, parent_cpy, original_values, n, NULL)) ;

        GRB_TRY (GrB_Matrix_new (&S, GrB_FP64, num_preserved, n)) ;

    } else { 
        // result dim: n by n
        GRB_TRY (GrB_Matrix_new (&S, GrB_FP64, n, n)) ;
        // newlabels is the identity map, signified by null return values
        if (newlabels != NULL) {
            (*newlabels) = NULL ;
        }
        if (inv_newlabels != NULL) {
            (*inv_newlabels) = NULL ;
        }
    }

    GRB_TRY (GxB_Vector_unpack_CSC (parent_cpy, &S_cols, (void**) &S_rows, &S_cols_size, &S_rows_size, NULL, &nvals, NULL, NULL)) ;
    GRB_TRY (GrB_Scalar_new (&one, GrB_FP64)) ;
    GRB_TRY (GrB_Scalar_setElement (one, 1)) ;

    GRB_TRY (GxB_Matrix_build_Scalar (S, S_rows, S_cols, one, nvals)) ;    

    (*result) = S ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

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
    LAGraph_Free ((void**)(&rows), msg) ;           \
    LAGraph_Free ((void**)(&cols), msg) ;           \
    LAGraph_Free ((void**)(&vals), msg) ;           \
}                                                   \

#define LG_FREE_ALL                                 \
{                                                   \
    LG_FREE_WORK ;                                  \
    LAGraph_Delete(&G_cpy, msg) ;                   \
    LAGraph_Free((void**) all_parents, msg) ;       \
}                                                   \

int LAGraph_Coarsen_Matching
(
    // outputs:
    GrB_Matrix *coarsened,                  // coarsened adjacency
    GrB_Vector **parent_result,             // array of parent mappings for each level; if preserve_mapping is true, is NULL
                                            // specifically, parent_result[i][u] = v if node u maps to node v in the i-th graph
    GrB_Vector **newlabels_result,          // array of mappings from old nodes to new nodes for each level; needed to interpret parent_result
                                            // if preserve_mapping is true, is NULL
                                            // refer to the description of the mapping result in Parent_to_S for the contents of newlabels_result[i]
    // inputs:
    LAGraph_Graph G,                        // input graph
    LAGraph_Matching_kind matching_type,    // how to perform the coarsening
    int preserve_mapping,                   // preserve original namespace of nodes
    int combine_weights,                    // whether to sum edge weights or just keep the pattern
    GrB_Index nlevels,                      // #of coarsening levels
    uint64_t seed,                          // used for matching
    char *msg
)
{
    LG_CLEAR_MSG ;
    LAGraph_Graph G_cpy = NULL ;            // used for the IncidenceMatrix function
    GrB_Matrix A = NULL ;                   // resulting adjacency matrix (used for output)
    GrB_Matrix E = NULL ;                   // incidence matrix
    GrB_Matrix E_t = NULL ;                 // transpose of incidence
    GrB_Matrix S = NULL ;                   // S matrix (S[i][j] -> node j maps to node i in coarsened graph)
    GrB_Matrix S_t = NULL ;                 // transpose of S
    GrB_Vector matched_edges = NULL ;       // result of maximal matching
    GrB_Vector edge_parent = NULL ;         // points to parent (representative) node for each edge
    GrB_Vector node_parent = NULL ;         // points to parent (representative) node for each node
    GrB_Vector full = NULL ;                // full vector

    GrB_Vector *all_parents = NULL ;        // resulting array of parents (used for output)
    GrB_Vector *all_newlabels = NULL ;       // resulting array of mappings (used for output)

    // used to build int64 A matrix if needed
    GrB_Index *rows = NULL ;
    GrB_Index *cols = NULL ;
    int64_t *vals = NULL ;
    GrB_Index nvals, nrows ;
    GrB_Type A_type ;
    // check properties (no self-loops, undirected)
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED)
    {
        char typename[LAGRAPH_MAX_NAME_LEN] ;
        GrB_Type type ;
        LG_TRY (LAGraph_Matrix_TypeName (typename, G->A, msg)) ;
        LG_TRY (LAGraph_TypeFromName (&type, typename, msg)) ;

        if ((type == GrB_FP32 || type == GrB_FP64) || (type == GrB_INT64 || type == GrB_UINT64)) {
            // output will keep the same type as input
            GRB_TRY (GrB_Matrix_dup (&A, G->A)) ;
            A_type = type ;
        } else {
            // output will become int64
            // reasoning: want to prevent overflow from combining edges and accomodate negative edge weights
            #ifdef dbg
                printf("Rebuilding A with GrB_INT64, orig type was %s\n", typename);
            #endif

            GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;
            GRB_TRY (GrB_Matrix_nrows (&nrows, G->A)) ;

            LG_TRY (LAGraph_Malloc ((void**)(&rows), nvals, sizeof(GrB_Index), msg)) ;
            LG_TRY (LAGraph_Malloc ((void**)(&cols), nvals, sizeof(GrB_Index), msg)) ;
            LG_TRY (LAGraph_Malloc ((void**)(&vals), nvals, sizeof(int64_t), msg)) ;
            // extractTuples casts all entries to int64
            GRB_TRY (GrB_Matrix_extractTuples (rows, cols, vals, &nvals, G->A)) ;

            GRB_TRY (GrB_Matrix_new (&A, GrB_INT64, nrows, nrows)) ;
            GRB_TRY (GrB_Matrix_build (A, rows, cols, vals, nvals, NULL)) ;

            LG_TRY (LAGraph_Free ((void**)(&rows), msg)) ;
            LG_TRY (LAGraph_Free ((void**)(&cols), msg)) ;
            LG_TRY (LAGraph_Free ((void**)(&vals), msg)) ;
            
            A_type = GrB_INT64 ;
        }
    }
    else
    {
        // G is not undirected
        LG_ASSERT_MSG (false, -105, "G must be undirected") ;
    }

    LG_ASSERT_MSG (G->nself_edges == 0, -107, "G->nself_edges must be zero") ;
    
    // make new LAGraph_Graph to use for building incidence matrix and for useful functions (delete self-edges)
    LG_TRY (LAGraph_New (&G_cpy, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    LG_TRY (LAGraph_Cached_NSelfEdges (G_cpy, msg)) ;

    // set A back (LAGraph_New sets A = NULL)
    A = G_cpy->A ;

    GrB_Index num_nodes ;
    GrB_Index num_edges ;

    GRB_TRY (GrB_Matrix_nrows (&num_nodes, A)) ;

    if (preserve_mapping) {
        GRB_TRY (GrB_Matrix_new (&S_t, GrB_UINT64, num_nodes, num_nodes)) ;
        GRB_TRY (GrB_Vector_new (&node_parent, GrB_UINT64, num_nodes)) ;
    }

    GRB_TRY (GrB_Vector_new (&full, GrB_BOOL, num_nodes)) ;

    GRB_TRY (GrB_assign (full, NULL, NULL, true, GrB_ALL, num_nodes, NULL)) ;

    LG_TRY (LAGraph_Malloc ((void**)(&all_parents), nlevels, sizeof(GrB_Vector), msg)) ;
    if (!preserve_mapping) {
        LG_TRY (LAGraph_Malloc ((void**)(&all_newlabels), nlevels, sizeof(GrB_Vector), msg)) ;
    }

    GrB_Index curr_level = 0 ;

    while (nlevels > 0) {
        // get E
        LG_TRY (LAGraph_Incidence_Matrix (&E, G_cpy, msg)) ;
        GRB_TRY (GrB_Matrix_nvals (&num_edges, A)) ;
        num_edges /= 2 ; // since undirected

        if (!preserve_mapping) {
            GRB_TRY (GrB_Matrix_nrows (&num_nodes, A)) ;

            // create node_parent for this level
            GRB_TRY (GrB_Vector_new (&node_parent, GrB_UINT64, num_nodes)) ;

            // ok to resize full since we are using its contents later
            GRB_TRY (GrB_Vector_resize (full, num_nodes)) ;
        }
        GRB_TRY (GrB_Matrix_new (&E_t, GrB_FP64, num_edges, num_nodes)) ;
        GRB_TRY (GrB_Vector_new (&edge_parent, GrB_UINT64, num_edges)) ;

        GRB_TRY (GrB_transpose (E_t, NULL, NULL, E, NULL)) ;

        // run maximal matching
        LG_TRY (LAGraph_MaximalMatching (&matched_edges, E, matching_type, seed, msg)) ;

        // TODO: Make single coarsening step a util function
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
        GRB_TRY (GrB_mxv (node_parent, NULL, NULL, GrB_MIN_SECOND_SEMIRING_UINT64, E, edge_parent, NULL)) ;

        // populate non-existent entries in node_parent with their index
        // handles nodes that are not engaged in a matching
        GrB_apply (node_parent, node_parent, NULL, GrB_ROWINDEX_INT64, full, (int64_t) 0, GrB_DESC_SC) ;

        // record a deep copy of the current node_parent for the current coarsening level
        GRB_TRY (GrB_Vector_dup (all_parents + curr_level, node_parent)) ;

        #ifdef dbg
            printf("Printing node_parent for level (%lld)\n", curr_level) ;
            LG_TRY (LAGraph_Vector_Print (node_parent, LAGraph_COMPLETE, stdout, msg)) ;
            printf("Printing matched edges for level (%lld)\n", curr_level) ;
            LG_TRY (LAGraph_Vector_Print (matched_edges, LAGraph_COMPLETE, stdout, msg)) ;
        #endif

        // build the S matrix
        LG_TRY (LAGraph_Parent_to_S (
            &S,
            (all_newlabels == NULL ? NULL : all_newlabels + curr_level), 
            NULL, 
            node_parent,
            preserve_mapping, 
            msg
        )) ;

        #ifdef dbg
            printf("Printing S for level (%lld)\n", curr_level) ;
            LG_TRY (LAGraph_Matrix_Print (S, LAGraph_COMPLETE, stdout, msg)) ;
            printf("Printing A for level (%lld)\n", curr_level) ;
            LG_TRY (LAGraph_Matrix_Print (A, LAGraph_COMPLETE, stdout, msg)) ;
        #endif

        GrB_Index S_rows, S_cols ;

        if (!preserve_mapping) {
            // need to create S_t for this level
            GRB_TRY (GrB_Matrix_nrows (&S_rows, S)) ;
            GRB_TRY (GrB_Matrix_ncols (&S_cols, S)) ;
            
            GRB_TRY (GrB_Matrix_new (&S_t, GrB_FP64, S_cols, S_rows)) ;
        }
        GRB_TRY (GrB_transpose (S_t, NULL, NULL, S, NULL)) ;

        GrB_Semiring semiring = combine_weights ? GrB_PLUS_TIMES_SEMIRING_FP64 : LAGraph_any_one_bool ;
        
        GRB_TRY (GrB_mxm (S, NULL, NULL, semiring, S, A, NULL)) ;

        #ifdef dbg
            printf("Printing S * A for level (%lld)\n", curr_level) ;
            LG_TRY (LAGraph_Matrix_Print (S, LAGraph_COMPLETE, stdout, msg)) ;
        #endif

        if (!preserve_mapping) {
            // resize result
            GRB_TRY (GrB_free (&A)) ;
            GRB_TRY (GrB_Matrix_new (&A, A_type, S_rows, S_rows)) ;
        }
        GRB_TRY (GrB_mxm (A, NULL, NULL, semiring, S, S_t, NULL)) ;

        G_cpy->A = A ;
        // make nself_edges unknown for delete self edges to work
        G_cpy->nself_edges = LAGRAPH_UNKNOWN ;
        // parent nodes for matched edges will form self-edges; need to delete
        LG_TRY (LAGraph_DeleteSelfEdges (G_cpy, msg)) ;
        A = G_cpy->A ;

        // want to free before we reassign what they point to
        GRB_TRY (GrB_free (&S)) ;
        GRB_TRY (GrB_free (&E)) ;
        GRB_TRY (GrB_free (&E_t)) ;
        GRB_TRY (GrB_free (&matched_edges)) ;
        GRB_TRY (GrB_free (&edge_parent)) ;

        if (!preserve_mapping){
            // also free node_parent and S_t for this level
            GRB_TRY (GrB_free (&node_parent)) ;
            GRB_TRY (GrB_free (&S_t)) ;
        }

        nlevels-- ;
        curr_level++ ;
    }
    (*coarsened) = A ;
    (*parent_result) = all_parents ;
    (*newlabels_result) = all_newlabels ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
