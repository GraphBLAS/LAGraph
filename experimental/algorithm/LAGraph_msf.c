//------------------------------------------------------------------------------
// LAGraph_msf.c
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Yongzhe Zhang (zyz915@gmail.com)
// Revised by Gabriel Gomez and Tim Davis

//------------------------------------------------------------------------------

/**
 * Code is based on Boruvka's minimum spanning forest algorithm
 * 
 * The algorithm calculates the minimum spanning forest of a graph represented 
 * by matrix A. It outputs a matrix containing lowest weight edges which connect 
 * the nodes in the forest, and a vector indicating the connected component 
 * "tree" which each node belongs to.
 */

// TODO: not ready for src but getting close.

// TODO: a "sanitize" input is fine for now in the experimental folder, but it
// doesn't fit with the standard LAGraph API.  It will need to be removed when
// this method is moved to the src folder.  The input will also become an
// LAGraph_Graph, not a plain GrB_Matrix A.

// FUTURE:  the idx type could be INT32 if the matrix has < 2 billion nodes,
// and the weight types could be expanded (int32 and float).

#include "LG_internal.h"
#include <LAGraph.h>
#include <LAGraphX.h>

//------------------------------------------------------------------------------
// tuple: a tuple containing (weight,index)
//------------------------------------------------------------------------------

// LG_MSF_tuple_int is used if the input graph uses integer weights of any type;
// LG_MSF_tuple_fp is used if the input graph is FP32 or FP64.  Likewise for the
// other *_int and *_fp types and operators.

typedef struct
{
    int64_t wInt;
    uint64_t idx;
} LG_MSF_tuple_int;

#define TUPLE_INT      \
"typedef struct    \n" \
"{                 \n" \
"    int64_t wInt; \n" \
"    uint64_t idx; \n" \
"} LG_MSF_tuple_int;"

typedef struct
{
    double wFp;
    uint64_t idx;
} LG_MSF_tuple_fp;

#define TUPLE_FP       \
"typedef struct    \n" \
"{                 \n" \
"    double wFp;   \n" \
"    uint64_t idx; \n" \
"} LG_MSF_tuple_fp;"

//------------------------------------------------------------------------------
// context_type: context for IndexUnaryOps (using the theta input)
//------------------------------------------------------------------------------

typedef struct
{
    uint64_t    *parent;   // parent of each vertex in the spanning forest
    struct
    {
        int64_t wInt;
        uint64_t idx;
    } *w_partner;          // partner vertex in the spanning forest
} LG_MSF_context_int;

#define LG_MSF_CONTEXT_INT    \
"typedef struct           \n" \
"{                        \n" \
"    uint64_t    *parent; \n" \
"    struct               \n" \
"    {                    \n" \
"        int64_t wInt;    \n" \
"        uint64_t idx;    \n" \
"    } *w_partner;        \n" \
"} LG_MSF_context_int;"

typedef struct
{
    uint64_t    *parent;   // parent of each vertex in the spanning forest
    struct
    {
        double wFp;
        uint64_t idx;
    } *w_partner;          // partner vertex in the spanning forest
} LG_MSF_context_fp;

#define LG_MSF_CONTEXT_FP     \
"typedef struct           \n" \
"{                        \n" \
"    uint64_t    *parent; \n" \
"    struct               \n" \
"    {                    \n" \
"        double wFp;      \n" \
"        uint64_t idx;    \n" \
"    } *w_partner;        \n" \
"} LG_MSF_context_fp;"

//------------------------------------------------------------------------------
// selectEdge: index-unary operator to select edges of min weight
//------------------------------------------------------------------------------

// generate solution:
// for each element A(i, j), it is selected if
//   1. weight[i] == A(i, j)    -- where weight[i] stores i's minimum edge weight
//   2. parent[j] == partner[i] -- j belongs to the specified connected component

void LG_MSF_selectEdge_int 
(
    bool *z, 
    const int64_t *x, 
    GrB_Index i, 
    GrB_Index j, 
    const LG_MSF_context_int *theta
)
{
    (*z) = (theta->w_partner[i].wInt == *x) && 
        (theta->parent[j] == theta->w_partner[i].idx);
}

#define SELECTEDGE_INT \
"void LG_MSF_selectEdge_int                            \n"\
"(                                                     \n"\
"    bool *z,                                          \n"\
"    const int64_t *x,                                 \n"\
"    GrB_Index i,                                      \n"\
"    GrB_Index j,                                      \n"\
"    const LG_MSF_context_int *theta                   \n"\
")                                                     \n"\
"{                                                     \n"\
"    (*z) = (theta->w_partner[i].wInt == *x) &&        \n"\
"        (theta->parent[j] == theta->w_partner[i].idx);\n"\
"}"

void LG_MSF_selectEdge_fp 
(
    bool *z, 
    const double *x, 
    GrB_Index i, 
    GrB_Index j, 
    const LG_MSF_context_fp *theta
)
{
    (*z) = (theta->w_partner[i].wFp == *x) && 
        (theta->parent[j] == theta->w_partner[i].idx);
}

#define SELECTEDGE_FP \
"void LG_MSF_selectEdge_fp                             \n"\
"(                                                     \n"\
"    bool *z,                                          \n"\
"    const double *x,                                  \n"\
"    GrB_Index i,                                      \n"\
"    GrB_Index j,                                      \n"\
"    const LG_MSF_context_fp *theta                    \n"\
")                                                     \n"\
"{                                                     \n"\
"    (*z) = (theta->w_partner[i].wFp == *x) &&         \n"\
"        (theta->parent[j] == theta->w_partner[i].idx);\n"\
"}"

//------------------------------------------------------------------------------
// removeEdge: remove edge (i,j) when i and j have the same parent
//------------------------------------------------------------------------------

// edge removal:
// A(i, j) is removed when parent[i] == parent[j]

void LG_MSF_removeEdge_int 
(
    bool *z, 
    const int64_t *x, 
    GrB_Index i, 
    GrB_Index j, 
    const LG_MSF_context_int *theta
)
{
    (*z) = (theta->parent[i] != theta->parent[j]);
}

#define REMOVEEDGE_INT \
"void LG_MSF_removeEdge_int                        \n"\
"(                                                 \n"\
"    bool *z,                                      \n"\
"    const int64_t *x,                             \n"\
"    GrB_Index i,                                  \n"\
"    GrB_Index j,                                  \n"\
"    const LG_MSF_context_int *theta               \n"\
")                                                 \n"\
"{                                                 \n"\
"    (*z) = (theta->parent[i] != theta->parent[j]);\n"\
"}"

void LG_MSF_removeEdge_fp 
(
    bool *z, 
    const double *x, 
    GrB_Index i, 
    GrB_Index j, 
    const LG_MSF_context_fp *theta
)
{
    (*z) = (theta->parent[i] != theta->parent[j]);
}

#define REMOVEEDGE_FP \
"void LG_MSF_removeEdge_fp                         \n"\
"(                                                 \n"\
"    bool *z,                                      \n"\
"    const double *x,                              \n"\
"    GrB_Index i,                                  \n"\
"    GrB_Index j,                                  \n"\
"    const LG_MSF_context_fp *theta                \n"\
")                                                 \n"\
"{                                                 \n"\
"    (*z) = (theta->parent[i] != theta->parent[j]);\n"\
"}"

//------------------------------------------------------------------------------
// combine: create a tuple from a weight and an index
//------------------------------------------------------------------------------

void LG_MSF_combine_int 
(
    LG_MSF_tuple_int *z, 
    const int64_t *x, 
    const uint64_t *y
)
{
    z->wInt = *x;
    z->idx = *y;
}

#define COMBINE_INT \
"void LG_MSF_combine_int  \n"\
"(                        \n"\
"    LG_MSF_tuple_int *z, \n"\
"    const int64_t *x,    \n"\
"    const uint64_t *y    \n"\
")                        \n"\
"{                        \n"\
"    z->wInt = *x;        \n"\
"    z->idx = *y;         \n"\
"}"

void LG_MSF_combine_fp 
(
    LG_MSF_tuple_fp *z, 
    const double *x, 
    const uint64_t *y
)
{
    z->wFp = *x;
    z->idx = *y;
}

#define COMBINE_FP \
"void LG_MSF_combine_fp  \n"\
"(                       \n"\
"    LG_MSF_tuple_fp *z, \n"\
"    const double *x,    \n"\
"    const uint64_t *y   \n"\
")                       \n"\
"{                       \n"\
"    z->wFp = *x;        \n"\
"    z->idx = *y;        \n"\
"}"

//------------------------------------------------------------------------------
// get_first:  get first item in a tuple (the weight)
//------------------------------------------------------------------------------

void LG_MSF_get_first_int (int64_t *y, const LG_MSF_tuple_int *x)
{
    *y = x->wInt;
}

#define GET_FIRST_INT \
"void LG_MSF_get_first_int (int64_t *y, const LG_MSF_tuple_int *x)  \n" \
"{                                                                  \n" \
"    *y = x->wInt;                                                  \n" \
"}"

void LG_MSF_get_first_fp (double *y, const LG_MSF_tuple_fp *x)
{
    *y = x->wFp;
}

#define GET_FIRST_FP \
"void LG_MSF_get_first_fp (double *y, const LG_MSF_tuple_fp *x)   \n" \
"{                                                                \n" \
"    *y = x->wFp;                                                 \n" \
"}"

//------------------------------------------------------------------------------
// get_second:  get second item in a tuple (the index)
//------------------------------------------------------------------------------

void LG_MSF_get_second_int (uint64_t *y, const LG_MSF_tuple_int *x)
{
    *y = x->idx;
}

#define GET_SECOND_INT \
"void LG_MSF_get_second_int (uint64_t *y, const LG_MSF_tuple_int *x)  \n" \
"{                                                                    \n" \
"    *y = x->idx;                                                     \n" \
"}"

void LG_MSF_get_second_fp (uint64_t *y, const LG_MSF_tuple_fp *x)
{
    *y = x->idx;
}

#define GET_SECOND_FP \
"void LG_MSF_get_second_fp (uint64_t *y, const LG_MSF_tuple_fp *x)    \n" \
"{                                                                    \n" \
"    *y = x->idx;                                                     \n" \
"}"

//------------------------------------------------------------------------------
// tupleMin: z = the min tuple of x and y
//------------------------------------------------------------------------------

void LG_MSF_tupleMin_int 
(
    LG_MSF_tuple_int *z, 
    const LG_MSF_tuple_int *x, 
    const LG_MSF_tuple_int *y
)
{
    bool xSmaller = x->wInt < y->wInt || 
        (x->wInt == y->wInt && x->idx < y->idx);
    z->wInt = (xSmaller)? x->wInt: y->wInt;
    z->idx = (xSmaller)? x->idx: y->idx;
}

#define TUPLEMIN_INT \
"void LG_MSF_tupleMin_int                        \n"\
"(                                               \n"\
"    LG_MSF_tuple_int *z,                        \n"\
"    const LG_MSF_tuple_int *x,                  \n"\
"    const LG_MSF_tuple_int *y                   \n"\
")                                               \n"\
"{                                               \n"\
"    bool xSmaller = x->wInt < y->wInt ||        \n"\
"        (x->wInt == y->wInt && x->idx < y->idx);\n"\
"    z->wInt = (xSmaller)? x->wInt: y->wInt;     \n"\
"    z->idx = (xSmaller)? x->idx: y->idx;        \n"\
"}"

void LG_MSF_tupleMin_fp 
(
    LG_MSF_tuple_fp *z, 
    const LG_MSF_tuple_fp *x, 
    const LG_MSF_tuple_fp *y
)
{
    bool xSmaller = x->wFp < y->wFp || (x->wFp == y->wFp && x->idx < y->idx);
    z->wFp = (xSmaller)? x->wFp: y->wFp;
    z->idx = (xSmaller)? x->idx: y->idx;
}

#define TUPLEMIN_FP \
"void LG_MSF_tupleMin_fp                                                      \n"\
"(                                                                            \n"\
"    LG_MSF_tuple_fp *z,                                                      \n"\
"    const LG_MSF_tuple_fp *x,                                                \n"\
"    const LG_MSF_tuple_fp *y                                                 \n"\
")                                                                            \n"\
"{                                                                            \n"\
"    bool xSmaller = x->wFp < y->wFp || (x->wFp == y->wFp && x->idx < y->idx);\n"\
"    z->wFp = (xSmaller)? x->wFp: y->wFp;                                     \n"\
"    z->idx = (xSmaller)? x->idx: y->idx;                                     \n"\
"}"

//------------------------------------------------------------------------------
// tuple2nd: z = y
//------------------------------------------------------------------------------

void LG_MSF_tuple2nd_int 
(
    LG_MSF_tuple_int *z, 
    const void *x, 
    const LG_MSF_tuple_int *y
)
{
    z->wInt = y->wInt;
    z->idx = y->idx;
}

#define TUPLE2ND_INT \
"void LG_MSF_tuple2nd_int     \n"\
"(                            \n"\
"    LG_MSF_tuple_int *z,     \n"\
"    const void *x,           \n"\
"    const LG_MSF_tuple_int *y\n"\
")                            \n"\
"{                            \n"\
"    z->wInt = y->wInt;       \n"\
"    z->idx = y->idx;         \n"\
"}"

void LG_MSF_tuple2nd_fp 
(
    LG_MSF_tuple_fp *z, 
    const void *x, 
    const LG_MSF_tuple_fp *y
)
{
    z->wFp = y->wFp;
    z->idx = y->idx;
}

#define TUPLE2ND_FP \
"void LG_MSF_tuple2nd_fp     \n"\
"(                           \n"\
"    LG_MSF_tuple_fp *z,     \n"\
"    const void *x,          \n"\
"    const LG_MSF_tuple_fp *y\n"\
")                           \n"\
"{                           \n"\
"    z->wFp = y->wFp;        \n"\
"    z->idx = y->idx;        \n"\
"}"

//------------------------------------------------------------------------------
// tupleEq: true if two tuples are equal
//------------------------------------------------------------------------------

void LG_MSF_tupleEq_int 
(
    bool *z, 
    const LG_MSF_tuple_int *x, 
    const LG_MSF_tuple_int *y
)
{
    *z = (x->wInt == y->wInt) && (x->idx == y->idx);
}

#define TUPLEEQ_INT \
"void LG_MSF_tupleEq_int                             \n"\
"(                                                   \n"\
"    bool *z,                                        \n"\
"    const LG_MSF_tuple_int *x,                      \n"\
"    const LG_MSF_tuple_int *y                       \n"\
")                                                   \n"\
"{                                                   \n"\
"    *z = (x->wInt == y->wInt) && (x->idx == y->idx);\n"\
"}"

void LG_MSF_tupleEq_fp 
(
    bool *z, 
    const LG_MSF_tuple_fp *x, 
    const LG_MSF_tuple_fp *y
)
{
    *z = (x->wFp == y->wFp) && (x->idx == y->idx);
}

#define TUPLEEQ_FP \
"void LG_MSF_tupleEq_fp                            \n"\
"(                                                 \n"\
"    bool *z,                                      \n"\
"    const LG_MSF_tuple_fp *x,                     \n"\
"    const LG_MSF_tuple_fp *y                      \n"\
")                                                 \n"\
"{                                                 \n"\
"    *z = (x->wFp == y->wFp) && (x->idx == y->idx);\n"\
"}"

//------------------------------------------------------------------------------

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                 \
{                                                   \
    GrB_free (&S);                                  \
    GrB_free (&T);                                  \
    LAGraph_Free ((void **) &SI, msg);              \
    LAGraph_Free ((void **) &SJ, msg);              \
    LAGraph_Free ((void **) &SX, msg);              \
    LAGraph_Free ((void **) &context_int.parent, msg);      \
    LAGraph_Free ((void **) &context_fp.parent, msg);       \
    GrB_free (&f);                                  \
    GrB_free (&w);                                  \
    GrB_free (&I);                                  \
    GrB_free (&t);                                  \
    GrB_free (&edge);                               \
    GrB_free (&cedge);                              \
    GrB_free (&tedge);                              \
    GrB_free (&mask);                               \
    GrB_free (&index_v);                            \
    GrB_free (&combine);                            \
    GrB_free (&minComb);                            \
    GrB_free (&get_first);                          \
    GrB_free (&get_second);                         \
    GrB_free (&selectEdge);                         \
    GrB_free (&removeEdge);                         \
    GrB_free (&context_type);                       \
    GrB_free (&parent_v);                           \
    GrB_free (&ramp);                               \
    GrB_free (&tupleMin);                           \
    GrB_free (&tuple2nd);                           \
    GrB_free (&tupleEq);                            \
    GrB_free (&tupleMin_monoid);                    \
    GrB_free (&tupleMin2nd);                        \
    GrB_free (&tuple);                              \
    GrB_free (&max_weight);                         \
}

#ifdef DBG
#undef  GRB_CATCH
#define GRB_CATCH(info)                                                 \
{                                                                       \
    printf ("GraphBLAS failure (file %s, line %d): info: %d",           \
        __FILE__, __LINE__, info) ;                                     \
    LG_ERROR_MSG ("GraphBLAS failure (file %s, line %d): info: %d",     \
        __FILE__, __LINE__, info) ;                                     \
    LG_FREE_ALL ;                                                       \
    return (info) ;                                                     \
}
#endif

//------------------------------------------------------------------------------
// dump_tuple_vector: debugging only
//------------------------------------------------------------------------------

// #define DBG

#ifdef DBG
static void dump_tuple_vector
(
    char *vname,
    GrB_Vector v,           // of type tuple (int or fp)
    GrB_Type weight_type
)
{
    GrB_Index n = 0 ;
    GrB_Info info ;
    GrB_Vector_size (&n, v) ;
    printf ("\ntuple vector %s, size %lu\n", vname, n) ;
    if (weight_type == GrB_INT64)
    {
        printf ("weight type: int64\n") ;
        LG_MSF_tuple_int e ;
        for (int i = 0 ; i < n ; i++)
        {
            info = GrB_Vector_extractElement_UDT (&e, v, i) ;
            if (info == GrB_SUCCESS)
            {
                printf ("   (%d) (%ld, %lu)\n", i, e.wInt, e.idx) ;
            }
        }
    }
    else
    {
        printf ("weight type: double\n") ;
        LG_MSF_tuple_fp e ;
        for (int i = 0 ; i < n ; i++)
        {
            info = GrB_Vector_extractElement_UDT (&e, v, i) ;
            if (info == GrB_SUCCESS)
            {
                printf ("   (%d) (%g, %lu)\n", i, e.wFp, e.idx) ;
            }
        }
    }
}
#endif

//------------------------------------------------------------------------------
// LAGraph_msf
//------------------------------------------------------------------------------

int LAGraph_msf
(
    GrB_Matrix *forest_edges, // output: an unsymmetrical matrix, containing
                        // the edges in the spanning forest
    GrB_Vector *componentId,  // output: The connected component of each node
                        // componentId[i] is the representative node of the 
                        // component set that i is in.
    GrB_Matrix A,       // input matrix
    bool sanitize,      // if true, ensure A is symmetric
    char *msg
)
{
    #if LG_SUITESPARSE_GRAPHBLAS_V10
    LG_CLEAR_MSG ;

    LG_MSF_context_int context_int = {.parent = NULL, .w_partner = NULL } ;
    LG_MSF_context_fp  context_fp  = {.parent = NULL, .w_partner = NULL } ;
    LG_MSF_tuple_int inf_int = {.wInt = INT64_MAX, .idx = UINT64_MAX};
    LG_MSF_tuple_fp  inf_fp  = {.wFp  = INFINITY , .idx = UINT64_MAX};

    GrB_Info info;
    GrB_Index n;
    GrB_Matrix S = NULL, T = NULL;
    GrB_Vector f = NULL, I = NULL, t = NULL, parent_v = NULL, tedge = NULL,
        edge = NULL, cedge = NULL, mask = NULL, index_v = NULL, ramp = NULL,
        w = NULL ;

    GrB_Index *SI = NULL, *SJ = NULL;
    void *SX = NULL;
    GrB_Type context_type = NULL, tuple = NULL, weight_type = NULL, 
        ignore = NULL ;
    GrB_BinaryOp combine = NULL, tupleMin = NULL, tuple2nd = NULL, 
        tupleEq = NULL;
    GrB_Monoid tupleMin_monoid = NULL;
    GrB_Semiring minComb = NULL, tupleMin2nd = NULL;
    GrB_UnaryOp get_first = NULL, get_second = NULL;
    GrB_Scalar max_weight = NULL;
    int edge_handling = GrB_DEFAULT;
    uint64_t edge_size = 0, edge_n = 0;
    GrB_IndexUnaryOp selectEdge = NULL, removeEdge = NULL;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------

    if (forest_edges == NULL || A == NULL) return (GrB_NULL_POINTER) ;
    GrB_Index ncols ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    LG_ASSERT(n == ncols, GrB_DIMENSION_MISMATCH) ;

    GrB_BinaryOp min_weight = NULL ;
    size_t sx_size = 0 ;

    GrB_Type_Code tcode = 0;
    GrB_Matrix_get_INT32(A, (int *) &(tcode), GrB_EL_TYPE_CODE);
    switch (tcode)
    {
        case GrB_INT8_CODE:
        case GrB_INT16_CODE:
        case GrB_INT32_CODE:
        case GrB_INT64_CODE:
        case GrB_BOOL_CODE:
        case GrB_UINT8_CODE:
        case GrB_UINT16_CODE:
        case GrB_UINT32_CODE:
        case GrB_UINT64_CODE:
            // integer edge weights: use INT64
            weight_type = GrB_INT64 ;
            min_weight = GrB_MIN_INT64 ;
            sx_size = sizeof (int64_t) ;
            break;

        case GrB_FP32_CODE:
        case GrB_FP64_CODE:
            // floating-point edge weights: use FP64
            weight_type = GrB_FP64 ;
            min_weight = GrB_MIN_FP64 ;
            sx_size = sizeof (double) ;
            break;

        default:
            // other types are not supported
            LG_ASSERT(false, GrB_DOMAIN_MISMATCH) ;
            break;
    }

    //--------------------------------------------------------------------------
    // create types and operators
    //--------------------------------------------------------------------------

    GRB_TRY (GxB_Scalar_new(&max_weight, weight_type)) ;
    void *inf = NULL ;

    if (weight_type == GrB_INT64)
    {

        //----------------------------------------------------------------------
        // types and ops for INT64 weights
        //----------------------------------------------------------------------

        GRB_TRY (GxB_Type_new (&tuple, sizeof (LG_MSF_tuple_int),
            "LG_MSF_tuple_int", TUPLE_INT)) ;

        GRB_TRY (GxB_BinaryOp_new (
            &combine, (GxB_binary_function) LG_MSF_combine_int,
            tuple, weight_type, GrB_UINT64,
            "LG_MSF_combine_int", COMBINE_INT)) ;

        GRB_TRY (GxB_Scalar_setElement_INT64(max_weight, INT64_MAX)) ;

        GRB_TRY (GxB_BinaryOp_new (
            &tupleMin, (GxB_binary_function) LG_MSF_tupleMin_int,
            tuple, tuple, tuple,
            "LG_MSF_tupleMin_int", TUPLEMIN_INT)) ;

        GRB_TRY (GxB_BinaryOp_new (
            &tuple2nd, (GxB_binary_function) LG_MSF_tuple2nd_int,
            tuple, GrB_BOOL, tuple,
            "LG_MSF_tuple2nd_int", TUPLE2ND_INT)) ;

        GRB_TRY (GxB_BinaryOp_new (
            &tupleEq, (GxB_binary_function) LG_MSF_tupleEq_int,
            GrB_BOOL, tuple, tuple,
            "LG_MSF_tupleEq_int", TUPLEEQ_INT)) ;

        inf = (void *) (&inf_int) ;

        GRB_TRY (GxB_UnaryOp_new (
            &get_first, (GxB_unary_function) LG_MSF_get_first_int, weight_type, 
            tuple, "LG_MSF_get_first_int", GET_FIRST_INT)) ;

        GRB_TRY (GxB_UnaryOp_new (
            &get_second, (GxB_unary_function) LG_MSF_get_second_int, GrB_UINT64, 
            tuple, "LG_MSF_get_second_int", GET_SECOND_INT)) ;

        // context type
        GRB_TRY (GxB_Type_new (
            &context_type, sizeof (LG_MSF_context_int),
            "LG_MSF_context_int", LG_MSF_CONTEXT_INT)) ;

        // ops for GrB_select
        GRB_TRY(GxB_IndexUnaryOp_new (
            &selectEdge, (GxB_index_unary_function) LG_MSF_selectEdge_int, 
            GrB_BOOL, weight_type, context_type, 
            "LG_MSF_selectEdge_int", SELECTEDGE_INT)) ;

        GRB_TRY(GxB_IndexUnaryOp_new (
            &removeEdge, (void *) LG_MSF_removeEdge_int, GrB_BOOL, weight_type, 
            context_type, "LG_MSF_removeEdge_int", REMOVEEDGE_INT)) ;

    }
    else
    {

        //----------------------------------------------------------------------
        // types and ops for FP64 weights
        //----------------------------------------------------------------------

        GRB_TRY (GxB_Type_new (&tuple, sizeof (LG_MSF_tuple_fp),
            "LG_MSF_tuple_fp", TUPLE_FP)) ;

        GRB_TRY (GxB_BinaryOp_new (
            &combine, (GxB_binary_function) LG_MSF_combine_fp,
            tuple, weight_type, GrB_UINT64,
            "LG_MSF_combine_fp", COMBINE_FP)) ;

        GRB_TRY (GxB_Scalar_setElement_FP64(max_weight, INFINITY)) ;

        GRB_TRY (GxB_BinaryOp_new (
            &tupleMin, (GxB_binary_function) LG_MSF_tupleMin_fp,
            tuple, tuple, tuple,
            "LG_MSF_tupleMin_fp", TUPLEMIN_FP)) ;

        GRB_TRY (GxB_BinaryOp_new (
            &tuple2nd, (GxB_binary_function) LG_MSF_tuple2nd_fp,
            tuple, GrB_BOOL, tuple,
            "LG_MSF_tuple2nd_fp", TUPLE2ND_FP)) ;

        GRB_TRY (GxB_BinaryOp_new (
            &tupleEq, (GxB_binary_function) LG_MSF_tupleEq_fp,
            GrB_BOOL, tuple, tuple,
            "LG_MSF_tupleEq_fp", TUPLEEQ_FP)) ;

        inf = (void *) (&inf_fp) ;

        GRB_TRY (GxB_UnaryOp_new (
            &get_first, (GxB_unary_function) LG_MSF_get_first_fp, weight_type, 
            tuple, "LG_MSF_get_first_fp", GET_FIRST_FP)) ;

        GRB_TRY (GxB_UnaryOp_new (
            &get_second, (GxB_unary_function) LG_MSF_get_second_fp, GrB_UINT64, 
            tuple, "LG_MSF_get_second_fp", GET_SECOND_FP)) ;

        GRB_TRY (GxB_Type_new (
            &context_type, sizeof (LG_MSF_context_fp),
            "LG_MSF_context_fp", LG_MSF_CONTEXT_FP)) ;

        // ops for GrB_select
        GRB_TRY(GxB_IndexUnaryOp_new (
            &selectEdge, (GxB_index_unary_function) LG_MSF_selectEdge_fp, 
            GrB_BOOL, weight_type, context_type,
            "LG_MSF_selectEdge_fp", SELECTEDGE_FP)) ;

        GRB_TRY(GxB_IndexUnaryOp_new (
            &removeEdge, (void *) LG_MSF_removeEdge_fp, GrB_BOOL, weight_type, 
            context_type, "LG_MSF_removeEdge_fp", REMOVEEDGE_FP)) ;
    }

    GRB_TRY (GrB_Monoid_new_UDT (&tupleMin_monoid, tupleMin, inf)) ;
    GRB_TRY (GrB_Semiring_new (&minComb, tupleMin_monoid, combine)) ;
    GRB_TRY (GrB_Semiring_new (&tupleMin2nd, tupleMin_monoid, tuple2nd)) ;

    //--------------------------------------------------------------------------
    // create matrices and vectors
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Matrix_new (&S, weight_type, n, n)) ;

    if (sanitize)
    {
        // S = A+A', and typecasting to weight_type
        GRB_TRY (GrB_Matrix_eWiseAdd_BinaryOp
            (S, NULL, NULL, min_weight, A, A, GrB_DESC_T1)) ;
    }
    else
    {
        // S = A, typecasting to weight_type, if necessary 
        GRB_TRY (GrB_Matrix_assign
                (S, NULL, NULL, A, GrB_ALL, n, GrB_ALL, n, NULL)) ;
    }

    GRB_TRY (GrB_Matrix_new (&T, weight_type, n, n)) ;
    GRB_TRY (GrB_Vector_new (&w, weight_type, n)) ;
    GRB_TRY (GrB_Vector_new (&t, GrB_UINT64, n)) ;
    GRB_TRY (GrB_Vector_new (&f, GrB_UINT64, n)) ;
    GRB_TRY (GrB_Vector_new (&ramp, GrB_INT64, n + 1)) ;
    GRB_TRY (GrB_Vector_new (&edge, tuple, n)) ;
    GRB_TRY (GrB_Vector_new (&cedge, tuple, n)) ;
    GRB_TRY (GrB_Vector_new (&tedge, tuple, n)) ;
    GRB_TRY (GrB_Vector_new (&mask, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new (&index_v, GrB_UINT64, n)) ;
    GRB_TRY (GrB_Vector_new (&parent_v, GrB_UINT64, n)) ;

    LG_TRY (LAGraph_Malloc  ((void **) &SI, 2*n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc  ((void **) &SJ, 2*n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc  (&SX, 2*n, sx_size, msg)) ;

    // prepare vectors
    GRB_TRY (GrB_Vector_assign_UINT64 (
        f, NULL, NULL, (uint64_t) 0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64 (
        f, NULL, NULL, GrB_ROWINDEX_INT64, f, (int64_t) 0, NULL)) ;
    GRB_TRY (GrB_Vector_dup (&I, f)) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (
        ramp, NULL, NULL, (uint64_t) 0, GrB_ALL, n + 1, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64 (
        ramp, NULL, NULL, GrB_ROWINDEX_INT64, ramp, (int64_t) 0, NULL)) ;

    // create context
    if (weight_type == GrB_INT64)
    {
        LG_TRY (LAGraph_Malloc
            ((void **) &context_int.parent, n, sizeof (uint64_t), msg)) ;
        GRB_TRY (GrB_Vector_extractTuples (NULL, context_int.parent, &n, f)) ;
        GRB_TRY (GxB_Vector_load(parent_v, (void **) &context_int.parent,
            GrB_UINT64, n, n * sizeof (uint64_t), GxB_IS_READONLY, NULL)) ;
    }
    else
    {
        LG_TRY (LAGraph_Malloc
            ((void **) &context_fp.parent, n, sizeof (double), msg)) ;
        GRB_TRY (GrB_Vector_extractTuples (NULL, context_fp.parent, &n, f)) ;
        GRB_TRY (GxB_Vector_load(parent_v, (void **) &context_fp.parent,
            GrB_UINT64, n,  n * sizeof (double), GxB_IS_READONLY, NULL)) ;
    }

    //--------------------------------------------------------------------------
    // the main computation
    //--------------------------------------------------------------------------

    GrB_Index nvals, ntuples = 0, num;
    bool diff = false;
    GRB_TRY (GrB_Matrix_nvals (&nvals, S)) ;
    for (int iters = 1; nvals > 0; iters++)
    {
        // every vertex points to a root vertex at the beginning
        // edge[u] = u's minimum edge (weight and index are encoded together)
        GRB_TRY (GrB_Vector_assign_UDT (
            edge, NULL, NULL, inf, GrB_ALL, 0, NULL)) ;

        // each edge looks at its adjacent edges and picks the one with the 
        // minimum weight. This weight is put into a tuple with the 
        // representative value of the connected componect the edge connects to
        GRB_TRY (GrB_mxv (edge, NULL, tupleMin, minComb, S, f, NULL)) ;

        // cedge[u] = children's minimum edge  | if u is a root
        //          = (max_weight, u)          | otherwise
        GRB_TRY (GrB_Vector_apply_BinaryOp1st_Scalar (
            cedge, NULL, NULL, combine, max_weight, I, NULL)) ;
        LG_TRY (LAGraph_FastAssign_Semiring(
            cedge, NULL, tupleMin, parent_v, edge, ramp, tupleMin2nd, NULL, msg
        )) ;

        // if (f[u] == u) f[u] := get_second(cedge[u])  -- the index
        GRB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_EQ_UINT64, f, I, NULL)) ;
        GRB_TRY (GrB_apply (
            f, mask, GrB_SECOND_UINT64, get_second, cedge, NULL)) ;

        // identify all the vertex pairs (u, v) where f[u] == v and f[v] == u
        // and then select the minimum of u, v as the new root;
        // if (f[f[i]] == i) f[i] = min(f[i], i)
        GRB_TRY (GxB_Vector_extract_Vector (t, NULL, NULL, f, f, NULL)) ;
        GRB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_EQ_UINT64, I, t, NULL)) ;
        GRB_TRY (GrB_assign (f, mask, GrB_MIN_UINT64, I, GrB_ALL, 0, NULL)) ;

        // five steps to generate the solution
        // 1. new roots (f[i] == i) revise their entries in cedge
        GRB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_EQ_UINT64, I, f, NULL)) ;
        GRB_TRY (GrB_assign (cedge, mask, NULL, inf, GrB_ALL, 0, NULL)) ;

        // 2. every vertex tries to know whether one of its edges is selected
        GRB_TRY (GxB_Vector_extract_Vector (
            tedge, NULL, NULL, cedge, parent_v, NULL)) ;
        GRB_TRY (GrB_eWiseMult (mask ,NULL, NULL, tupleEq, edge, tedge, NULL)) ;

        // 3. each root picks a vertex from its children to 
        // generate the solution
        GRB_TRY (GrB_assign (index_v, NULL, NULL, n, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_assign (index_v, mask, NULL, I, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_assign (t, NULL, NULL, n, GrB_ALL, 0, NULL)) ;
        LG_TRY (LAGraph_FastAssign_Semiring(
            t, NULL, GrB_MIN_UINT64, parent_v, index_v, ramp,
            GrB_MIN_SECOND_SEMIRING_UINT64, NULL, msg
        )) ;
        GRB_TRY (GxB_Vector_extract_Vector (
            index_v, NULL, NULL, t, parent_v, NULL)) ;
        GRB_TRY (GrB_eWiseMult (
            mask ,NULL, NULL, GrB_EQ_UINT64, I, index_v, NULL)) ;

        // 4. generate the select function
        if (weight_type == GrB_INT64)
        {
            GRB_TRY (GxB_Vector_unload(
                edge, (void **) &context_int.w_partner, &ignore, &edge_n,
                &edge_size, &edge_handling, NULL)) ;
            GRB_TRY (GrB_Matrix_select_UDT (T, NULL, NULL, selectEdge, S,
                &context_int, NULL)) ;
            GRB_TRY (GxB_Vector_load(
                edge, (void **) &context_int.w_partner, tuple, edge_n,
                edge_size, edge_handling, NULL)) ;
        }
        else
        {
            GRB_TRY (GxB_Vector_unload(
                edge, (void **) &context_fp.w_partner, &ignore, &edge_n, 
                &edge_size, &edge_handling, NULL)) ;
            GRB_TRY (GrB_Matrix_select_UDT (
                T, NULL, NULL, selectEdge, S, &context_fp, NULL)) ;
            GRB_TRY (GxB_Vector_load(
                edge, (void **) &context_fp.w_partner, tuple, edge_n, edge_size,
                edge_handling, NULL)) ;
        }

        GRB_TRY (GrB_Vector_clear (t)) ;

        // 5. the generated matrix may still have redundant edges remove the 
        //  duplicates by GrB_mxv and store them as tuples in (SI,SJ,SX)
        GRB_TRY (GrB_Vector_clear (edge)) ;
        GRB_TRY (GrB_mxv (edge, mask, tupleMin, minComb, T, I, NULL)) ;
        GRB_TRY (GrB_Vector_nvals (&num, edge)) ;
        GRB_TRY (GrB_apply (t, NULL, NULL, get_second, edge, NULL)) ;
        GRB_TRY (GrB_Vector_extractTuples (NULL, SJ + ntuples, &num, t)) ;
        GRB_TRY (GrB_apply (w, NULL, NULL, get_first, edge, NULL)) ;

        if (weight_type == GrB_INT64)
        {
            GRB_TRY (GrB_Vector_extractTuples_INT64 (
                SI + ntuples, ((int64_t *) SX) + ntuples, &num, w)) ;
        }
        else
        {
            GRB_TRY (GrB_Vector_extractTuples_FP64 (
                SI + ntuples, ((double *) SX) + ntuples, &num, w)) ;
        }

        ntuples += num;

        // path halving until every vertex points on a root
        do {
            GRB_TRY (GxB_Vector_extract_Vector (t, NULL, NULL, f, f, NULL)) ;
            GRB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_NE_UINT64, f, t, NULL)) ;
            GRB_TRY (GrB_Vector_reduce_BOOL (&diff, NULL, GrB_LOR_MONOID_BOOL, mask, NULL)) ;
            GrB_Vector temp = f; f = t; t = temp;
        } while (diff);

        // remove the edges in the same connected component
        if (weight_type == GrB_INT64)
        {
            GRB_TRY (GrB_Vector_extractTuples (NULL, context_int.parent, &n, f)) ;
            GRB_TRY (GrB_Matrix_select_UDT (S, NULL, NULL, removeEdge, S, &context_int, NULL)) ;
        }
        else
        {
            GRB_TRY (GrB_Vector_extractTuples (NULL, context_fp.parent, &n, f)) ;
            GRB_TRY (GrB_Matrix_select_UDT (S, NULL, NULL, removeEdge, S, &context_fp, NULL)) ;
        }

        GrB_Matrix_nvals (&nvals, S);
    }

    // create forest_edges
    GRB_TRY (GrB_Matrix_clear (T)) ;
    if (weight_type == GrB_INT64)
    {
        GRB_TRY (GrB_Matrix_build_INT64 (
            T, SI, SJ, (int64_t *) SX, ntuples, GxB_IGNORE_DUP)) ;
    }
    else
    {
        GRB_TRY (GrB_Matrix_build_FP64 (
            T, SI, SJ, (double *) SX, ntuples, GxB_IGNORE_DUP)) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    *forest_edges = T;
    T = NULL ;

    if(componentId != NULL)
    {
        *componentId = f;
        f = NULL;
    }

    LG_FREE_ALL;
    return (GrB_SUCCESS) ;
    #else
    return (GrB_NOT_IMPLEMENTED) ;
    #endif
}
