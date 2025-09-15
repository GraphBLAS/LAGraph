//------------------------------------------------------------------------------
// LAGraph_Incidence_Matrix: Given the adjacency matrix of an undirected graph with no 
// self-loops, builds its corresponding incidence matrix
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Vidith Madhu, Texas A&M University

// TODO: not ready for src; should handle all builtin types, with option to
// typecast to INT64, UINT64, or FP64 as is currently done.

// TODO: this method is required for MaximalMatching and CoarsenMatching

// Given an undirected graph G, construct the incidence matrix E.
/*
The incidence matrix E has size n-by-e where the
undirected graph G has n nodes and e edges.  If the kth edge of G is the edge
(i,j), then the column E(:,k) contains two entries:  E(i,k) and E(j,k), which
have the same value.  If the graph G is weighted, then both E(i,k) and E(j,k)
are equal to the weight of the (i,j) edge.  If G is unweighted, then both are
equal to 1 (and the matrix E is thus iso-valued).

The type of the result is always compatible with the type of the input graph, but
may be of a larger size; UINT8, UINT16, UINT32, INT8, INT16, INT32, INT64, and BOOL becomes
a UINT64. UINT64 remains as a UINT64. FP32 and FP64 both become FP64.

Note that complex types are NOT supported.
*/

#include "LG_internal.h"
#include "LAGraphX.h"

#define LOADTRICKIM
//#define dbg

#undef LG_FREE_ALL
#define LG_FREE_ALL                                           \
{                                                             \
   LAGraph_Free ((void**)(&row_indices), msg) ;               \
   LAGraph_Free ((void**)(&col_indices), msg) ;               \
   LAGraph_Free ((void**)(&values), msg) ;                    \
   LAGraph_Free ((void**)(&ramp), msg) ;                      \
   GrB_free (&con) ;                                          \
   GrB_free (&E_half) ;                                       \
   GrB_free (&A_tril) ;                                       \
   GrB_free (&x);                                             \
   GrB_free (&i);                                             \
   GrB_free (&j);                                             \
   GrB_free (&Ep);                                            \
   GrB_free (&Ei);                                            \
   GrB_free (&Ex);                                            \
   GrB_free (&fullx);                                         \
   GrB_free(&build_desc);                                     \
}                                                             


int LAGraph_Incidence_Matrix
(
    GrB_Matrix *result, // incidence
    LAGraph_Graph G, // must be undirected, no self-loops
    char *msg
)
{
#if LAGRAPH_SUITESPARSE
    
    GrB_Matrix E = NULL ;
    GrB_Matrix E_half = NULL ;
    GrB_Matrix A_tril = NULL ;
    GrB_Vector i = NULL, j = NULL, x = NULL, fullx = NULL;
    GrB_Vector Ep = NULL, Ei = NULL, Ex = NULL;
    GrB_Descriptor build_desc = NULL;
    GrB_Index *row_indices = NULL ;
    GrB_Index *col_indices = NULL ;
    void *values = NULL ;
    #if GxB_IMPLEMENTATION < GxB_VERSION (10,0,0)
    GrB_Vector con = NULL;  // unused, except by LG_FREE_ALL above
    #else
    GxB_Container con = NULL;
    #endif

    GrB_Index *ramp = NULL ;

    LG_ASSERT_MSG (
        G->kind == LAGraph_ADJACENCY_UNDIRECTED,
        LAGRAPH_INVALID_GRAPH, 
        "G must be undirected"
    ) ;

    LG_ASSERT_MSG (G->nself_edges == 0, LAGRAPH_NO_SELF_EDGES_ALLOWED, "G->nself_edges must be zero") ;

    GrB_Matrix A = G->A ;

    char typename[LAGRAPH_MAX_NAME_LEN] ;
    GrB_Type type ;
    #if GxB_IMPLEMENTATION < GxB_VERSION (10,0,0)
    LG_TRY (LAGraph_Matrix_TypeName (typename, A, msg)) ;
    LG_TRY (LAGraph_TypeFromName (&type, typename, msg)) ;
    #else 
    // No longer historical
    GRB_TRY (GxB_Matrix_type(&type, A));
    #endif

    GrB_Index nvals ;
    GrB_Index num_nodes, num_edges ;

    GRB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
    num_edges = nvals / 2 ;
    GRB_TRY (GrB_Matrix_nrows (&num_nodes, A)) ;

    GRB_TRY (GrB_Matrix_new (&A_tril, type, num_nodes, num_nodes)) ;
    GRB_TRY (GrB_Matrix_new (&E, type, num_nodes, num_edges)) ;
    GRB_TRY (GrB_Matrix_new (&E_half, type, num_nodes, num_edges)) ;

    // get just the lower triangular entries
    GRB_TRY (GrB_select (A_tril, NULL, NULL, GrB_TRIL, A, 0, NULL)) ;

    #if GxB_IMPLEMENTATION < GxB_VERSION (10,0,0)
    bool is_uint64 = (type == GrB_UINT64) ;
    bool is_float = ((type == GrB_FP32) || (type == GrB_FP64)) ;

    // which intermediate type to cast to
    // uint64_t : 0, double: 1, int64_t: 2
    // if the input matrix type is GrB_INT* or GrB_UINT{8|16|32}, it becomes a int64_t
    // if the input matrix type is GrB_UINT64, it becomes a uint64_t
    // if the input matrix type is GrB_FP{32|64}, it becomes a double
    int which_type = (is_uint64 ? 0 : (is_float ? 1 : 2)) ;

    size_t value_sizes[3] = {sizeof(uint64_t), sizeof(double), sizeof(int64_t)} ;

    // (*result) = E ;
    // return (GrB_SUCCESS) ;

    // Arrays to extract A into
    LG_TRY (LAGraph_Malloc ((void**)(&row_indices), num_edges, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&col_indices), num_edges, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&values), num_edges, value_sizes[which_type], msg)) ;

    switch (which_type){
        case 0:
            GRB_TRY (GrB_Matrix_extractTuples_UINT64 (row_indices, col_indices, values, &num_edges, A_tril)) ;
            break;
        case 1:
            GRB_TRY (GrB_Matrix_extractTuples_FP64 (row_indices, col_indices, values, &num_edges, A_tril)) ;
            break;
        case 2:
            GRB_TRY (GrB_Matrix_extractTuples_INT64 (row_indices, col_indices, values, &num_edges, A_tril)) ;
            break;
    }
    
    #ifdef dbg
        printf("Printing A_tril values\n");
        for (int i = 0; i < num_edges; i++){
            printf("%ld %ld %.5f\n", row_indices[i], col_indices[i], values[i]);
        }
    #endif
    LG_TRY (LAGraph_Malloc ((void**)(&ramp), num_edges, sizeof(GrB_Index), msg)) ;

    #pragma omp parallel for
    for (GrB_Index i = 0 ; i < num_edges ; i++) {
        ramp[i] = i ;
    }

    // build E_1 with (row_indices, ramp, values)
    // build E_2 with (col_indices, ramp, values)
    // E = E_1 + E_2
    switch (which_type) {
        case 0:
            GRB_TRY (GrB_Matrix_build_UINT64 (E_half, col_indices, ramp, values, num_edges, NULL)) ;
            GRB_TRY (GrB_Matrix_build_UINT64 (E, row_indices, ramp, values, num_edges, NULL)) ;
            break;
        case 1:
            GRB_TRY (GrB_Matrix_build_FP64 (E_half, col_indices, ramp, values, num_edges, NULL)) ;
            GRB_TRY (GrB_Matrix_build_FP64 (E, row_indices, ramp, values, num_edges, NULL)) ;
            break;
        case 2:
            GRB_TRY (GrB_Matrix_build_INT64 (E_half, col_indices, ramp, values, num_edges, NULL)) ;
            GRB_TRY (GrB_Matrix_build_INT64 (E, row_indices, ramp, values, num_edges, NULL)) ;
            break;
    }
    GRB_TRY (GrB_eWiseAdd (E, NULL, NULL, GrB_PLUS_FP64, E, E_half, NULL)) ;
    #else
    // The vector types and sizes get overriden by extract.
    GRB_TRY (GrB_Vector_new(&x, GrB_BOOL, 0)) ;
    GRB_TRY (GrB_Vector_new(&i, GrB_BOOL, 0)) ;
    GRB_TRY (GrB_Vector_new(&j, GrB_BOOL, 0)) ;
    GRB_TRY (GxB_Matrix_extractTuples_Vector(i, j, x, A_tril, NULL)) ;

    // FUTURE: x is always size nvals(A), and never iso-valued, even if A_tril
    // is iso-valued.  If A_tril is iso-valued, then Ex below could be
    // length 1 and con->iso could be false.

        // this load trick is quicker, but returns GrB_COLMAJOR
        #ifdef LOADTRICKIM
        GRB_TRY (GrB_Vector_new(&Ep, GrB_INT64, num_edges + 1)) ;
        GrB_Type ij_type = NULL;
        GRB_TRY (GxB_Vector_type(&ij_type, i));
        GRB_TRY (GrB_Vector_new(&Ei, ij_type, num_edges * 2)) ;
        GRB_TRY (GrB_assign(
            Ep, NULL, NULL, (int64_t) 1, GrB_ALL, 0, NULL));

        // Shuffle i and j into Ei.
        // Ei = [j[0], i[0], j[1], i[1], . . ., i[num_edges -1]]
        // Filling out Ei helps assign be much quicker.
        GRB_TRY (GrB_assign(
            Ei, NULL, NULL, (int64_t) 1, GrB_ALL, 0, NULL));
        GrB_Index stride[] = {(GrB_Index) 0, num_edges * 2 - 1, (GrB_Index) 2} ;
//      printf ("num_edges: %lu\n", num_edges) ;
//      printf ("stride: [%lu %lu %lu]\n", stride [0], stride [1], stride [2]) ;
        if (num_edges > 0)
        {
            GRB_TRY (GrB_Vector_assign(
                Ei, NULL, NULL, j, stride, GxB_STRIDE, NULL)) ;
            stride[GxB_BEGIN] = 1;
            GRB_TRY (GrB_Vector_assign(
                Ei, NULL, NULL, i, stride, GxB_STRIDE, NULL)) ;
        }

        // Ep = [0,2,4,...,2 * numedges]
        GRB_TRY (GrB_Vector_apply_IndexOp_INT64(
            Ep, NULL, NULL, GrB_ROWINDEX_INT64, Ep, (uint64_t) 0, NULL)) ;
        GRB_TRY (GrB_Vector_apply_BinaryOp2nd_INT64(
            Ep, NULL, NULL, GrB_TIMES_INT64, Ep, (int64_t) 2, NULL)) ;

        // Filling out Ex helps assign be much quicker.
        GRB_TRY (GrB_Vector_new(&Ex, type, num_edges * 2)) ;
        GRB_TRY (GrB_assign(
            Ex, NULL, NULL, (int64_t) 1, GrB_ALL, 0, NULL));
        stride[GxB_BEGIN] = 0;
        if (num_edges > 0)
        {
            GRB_TRY (GrB_Vector_assign(
                Ex, NULL, NULL, x, stride, GxB_STRIDE, NULL)) ;
            stride[GxB_BEGIN] = 1;
            GRB_TRY (GrB_Vector_assign(
                Ex, NULL, NULL, x, stride, GxB_STRIDE, NULL)) ;
        }

        //load up container
        GRB_TRY (GxB_Container_new(&con));
        GRB_TRY (GrB_free(&con->p));
        GRB_TRY (GrB_free(&con->i));
        GRB_TRY (GrB_free(&con->x));
        con->p = Ep; Ep = NULL ;
        con->i = Ei; Ei = NULL ;
        con->x = Ex; Ex = NULL ;
        con->format = GxB_SPARSE;
        con->orientation = GrB_COLMAJOR;
        con->nrows = num_nodes;
        con->ncols = num_edges;
        con->nvals = num_edges * 2;
        con->jumbled = false;
        con->iso = false;
        // Ep = [0,2,4,...,2 * numedges]
        // Ex = [x[0], x[0], x[1], x[1], . . ., x[num_edges -1]]
        // Ei = [j[0], i[0], j[1], i[1], . . ., i[num_edges -1]]
        // So each column k has two entries at j[k] and i[k] with values x[k]
        GRB_TRY (GxB_load_Matrix_from_Container(E, con, NULL));
        GRB_TRY (GrB_free(&con));
        #else
        GRB_TRY (GrB_Vector_new(
            &fullx, GrB_BOOL, num_edges)) ;
        GRB_TRY (GrB_assign(
            fullx, NULL, NULL, (bool) 1, GrB_ALL, 0, NULL));
        GRB_TRY (GrB_Descriptor_new(&build_desc));
        GRB_TRY (GrB_set(build_desc, GxB_USE_INDICES, GxB_COLINDEX_LIST));
        // fullx interpreted by index so is just a ramp.
        GRB_TRY (GxB_Matrix_build_Vector(E_half, j, fullx, x, NULL, build_desc));
        GRB_TRY (GxB_Matrix_build_Vector(E, i, fullx, x, NULL, build_desc));
        GRB_TRY (GrB_eWiseAdd (E, NULL, NULL, GrB_PLUS_FP64, E, E_half, NULL)) ;
        #endif
    #endif
    // LAGraph_Matrix_Print (E, LAGraph_COMPLETE, stdout, msg) ;
    LG_FREE_ALL ;
    
    (*result) = E ;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
