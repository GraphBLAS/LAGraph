//------------------------------------------------------------------------------
// LAGraph_LouvainSeq.c: Runs the Louvain Algorithm with Isolate Sets on a given Graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2024 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Olumayowa Olowomeye, Texas A&M University

//------------------------------------------------------------------------------

// Current Test File: experimental/test/test_louvain.c
#include "LG_internal.h"
#include <LAGraphX.h>
#include <LAGraph.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#if LG_SUITESPARSE_GRAPHBLAS_V10

// #define TIMING
#define DEBUG 0

#if DEBUG
#define dbg(x) \
    if (DEBUG) GxB_print(x, 5) \
#define err(x, info)                                    \
    if (!(info == GrB_SUCCESS || info == GrB_NO_VALUE)) \
    {                                                   \
        char **err;                                     \
        GrB_error(err, x);                              \
        printf("\ninfo: %d error: %s\n", info, err);    \
    }
#else
#define dbg(x)
#define err(x, info)
#endif

typedef struct Theta
{
    double *d;
    double m;
    uint64_t seed;
} Theta;
#define THETA_DEFN                    \
    "typedef struct Theta"            \
    "{"                               \
    "    double *d;"                  \
    "    double m;"                   \
    "    uint64_t seed;"              \
    "} Theta;"

typedef struct argmax_tup
{
    double score; /* change in modularity */
    int64_t comm; /* who */
    double tb;
} argmax_tup;

#define AM_TUP                                       \
    "typedef struct argmax_tup\n"                    \
    "{\n"                                            \
    "    double score; /* change in modularity */\n" \
    "    int64_t comm;    /* who */\n"               \
    "    double tb;\n"                               \
    "} argmax_tup;\n"

void make_argmax_tup(argmax_tup *z,
                     const double *x, GrB_Index ix, GrB_Index jx,
                     const double *y, GrB_Index iy, GrB_Index jy,
                     const void *theta)
{
    Theta *_theta = (Theta *)theta;
    double ki_to_C = *x;       // Edge weight from node i to community C
    double ki = _theta->d[ix]; // Degree of node i
    double kC = *y;            // Sum of degrees in community C
    double m = _theta->m;      // Total edge weight of graph (already divided by 2)

    // Compute modularity gain (delta Q)
    z->score = (ki_to_C / m) - ((ki * kC) / (m * m));

    // Generate tiebreaker (seed)
    uint64_t seed = _theta->seed + (*y + ix + iy + jy);
    seed ^= seed << 13;
    seed ^= seed >> 7;
    seed ^= seed << 17;

    z->comm = (int64_t)jx;
    z->tb = seed;
}

#define MAKE_AM_TUP                                                       \
    "void make_argmax_tup(argmax_tup *z,\n"                               \
    "                     const double *x, GrB_Index ix, GrB_Index jx,\n" \
    "                     const double *y, GrB_Index iy, GrB_Index jy,\n" \
    "                     const void *theta)\n"                           \
    "{\n"                                                                 \
    "    Theta *_theta = (Theta *)theta;\n"                               \
    "    double ki_to_C = *x;\n"                                          \
    "    double ki = _theta->d[ix];\n"                                    \
    "    double kC = *y;\n"                                               \
    "    double m = _theta->m;\n"                                         \
    "    z->score = (ki_to_C / m) - ((ki * kC) / (m * m));\n"             \
    "    uint64_t seed = _theta->seed + (*y + ix + iy + jy);\n"           \
    "    seed ^= seed << 13;\n"                                           \
    "    seed ^= seed >> 7;\n"                                            \
    "    seed ^= seed << 17;\n"                                           \
    "    z->comm = (int64_t)jx;\n"                                        \
    "    z->tb = seed;\n"                                                 \
    "}\n"

void argmax_op(argmax_tup *z, argmax_tup *x, argmax_tup *y)
{
    if (x->score > y->score)
    {
        z->score = x->score;
        z->comm = x->comm;
        z->tb = x->tb;
    }
    else if (x->score == y->score)
    {
        if (x->tb > y->tb)
        {
            z->score = x->score;
            z->comm = x->comm;
            z->tb = x->tb;
        }
        else
        {
            z->score = y->score;
            z->comm = y->comm;
            z->tb = y->tb;
        }
    }
    else
    {
        z->score = y->score;
        z->comm = y->comm;
        z->tb = y->tb;
    }
}
#define AM_OP                                                       \
    "void argmax_op(argmax_tup *z, argmax_tup *x, argmax_tup *y)\n" \
    "{\n"                                                           \
    "    if (x->score > y->score)\n"                                \
    "    {\n"                                                       \
    "        z->score = x->score;\n"                                \
    "        z->comm = x->comm;\n"                                  \
    "        z->tb = x->tb;\n"                                      \
    "    }\n"                                                       \
    "    else if (x->score == y->score)\n"                          \
    "    {\n"                                                       \
    "        if (x->tb > y->tb)\n"                                  \
    "        {\n"                                                   \
    "            z->score = x->score;\n"                            \
    "            z->comm = x->comm;\n"                              \
    "            z->tb = x->tb;\n"                                  \
    "        }\n"                                                   \
    "        else\n"                                                \
    "        {\n"                                                   \
    "            z->score = y->score;\n"                            \
    "            z->comm = y->comm;\n"                              \
    "            z->tb = y->tb;\n"                                  \
    "        }\n"                                                   \
    "    }\n"                                                       \
    "    else\n"                                                    \
    "    {\n"                                                       \
    "        z->score = y->score;\n"                                \
    "        z->comm = y->comm;\n"                                  \
    "        z->tb = y->tb;\n"                                      \
    "    }\n"                                                       \
    "}\n"
static GrB_Info build_argmax_operator(
    GrB_Type Theta_UDT, GrB_Type Tuple,
    double *d, double m, uint64_t seed,
    GrB_BinaryOp *MAKEAMTUP_Bop,
    GrB_BinaryOp *AM_Bop,
    GrB_Monoid *AM_mon,
    GrB_Semiring *AM_Semiring,
    GxB_IndexBinaryOp *MAKEAMTUP_op,
    GrB_Scalar *_0,
    char *msg)
{
    Theta theta_scalar;
    theta_scalar.d = d;
    theta_scalar.seed = seed;
    theta_scalar.m = m;

    GRB_TRY(GrB_Scalar_setElement_UDT(*_0, (void *)&theta_scalar));
    GRB_TRY(GxB_IndexBinaryOp_new(MAKEAMTUP_op,
                                  (GxB_index_binary_function)make_argmax_tup,
                                  Tuple, GrB_FP64, GrB_FP64, Theta_UDT,
                                  "make_argmax_tup", MAKE_AM_TUP));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(MAKEAMTUP_Bop, *MAKEAMTUP_op, *_0));

    argmax_tup id;
    memset(&id, 0, sizeof(argmax_tup));
    id.tb = 0;
    id.comm = INT64_MAX;
    id.score = (double)(-INFINITY);

    GRB_TRY(GxB_BinaryOp_new(AM_Bop,
                             (GxB_binary_function)argmax_op, Tuple, Tuple, Tuple,
                             "argmax_op", AM_OP));

    GRB_TRY(GrB_Monoid_new_UDT(AM_mon, *AM_Bop, &id));

    GRB_TRY(GrB_Semiring_new(AM_Semiring, *AM_mon, *MAKEAMTUP_Bop));

    return GrB_SUCCESS;
}

#define LG_FREE_EXTRACT                  \
    {                                    \
        GrB_free(&extract_k_if_gain_op); \
    }
void extract_k_if_gain(void *out, const void *in)
{
    const argmax_tup *a = in;
    int64_t *k_out = out;
    if (a->score > 0.0)
    {
        *k_out = a->comm;
    }
    else
    {
        *k_out = a->comm; 
    }
}
#define EXTRACT_K_IF_GAIN_SRC                                       \
    "void extract_k_if_gain(void *out, const void *in) {\n"         \
    "    const argmax_tup *a = in;\n"                               \
    "    int64_t *k_out = out;\n"                                   \
    "    if (a->score > 0.0) {\n"                                   \
    "        *k_out = a->comm;\n"                                   \
    "    } else {\n"                                                \
    "        *k_out = a->comm;\n" \
    "    }\n"                                                       \
    "}\n"
#undef LG_FREE_ALL
#define LG_FREE_ALL                          \
    {                                        \
        /* vectors */                        \
        GrB_free(&x);                        \
        GrB_free(&y);                        \
        GrB_free(&k);                        \
        GrB_free(&iset);                     \
        GrB_free(&Wy);                       \
        GrB_free(&k_values);                 \
        GrB_free(&Si_old);                   \
        GrB_free(&A);                        \
        GrB_free(&Miset);                    \
        GrB_free(&W);                        \
        GrB_free(&A_rows);                   \
        GxB_Container_free(&S_container);    \
        GxB_Container_free(&k_container);    \
        GrB_free(&ri);                       \
        GrB_free(&rv);                       \
        GrB_free(&AM_Semiring);              \
        GrB_free(&AM_mon);                   \
        GrB_free(&AM_Bop);                   \
        GrB_free(&MAKEAMTUP_Bop);            \
        GrB_free(&MAKEAMTUP_op);             \
        GrB_free(&argmax_0);                 \
        GrB_free(&Theta_UDT);                \
        GrB_Type_free(&Tuple);               \
        LAGraph_Free((void **)&d_copy, msg); \
        LAGraph_Free((void **)&c_copy, msg); \
    }
#endif

int LAGraph_LouvainIS(
    // output
    GrB_Matrix *S_result,
    uint64_t seed,
    // input
    LAGraph_Graph G,
    char *msg)
{
#if LG_SUITESPARSE_GRAPHBLAS_V10
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];
    GrB_set(GrB_GLOBAL, false, GxB_BURBLE);

    GrB_Descriptor ri = NULL;
    GrB_Descriptor rv = NULL;

    GrB_Vector iset = NULL;
    GrB_Vector k = NULL;
    GrB_Vector x = NULL;
    GrB_Vector y = NULL;
    GrB_Vector Wy = NULL;
    GrB_Vector k_values = NULL;
    GrB_Vector Si_old = NULL;
    GrB_Vector gain_values = NULL;
    GrB_Vector gain_mask = NULL;

    GrB_Matrix S = NULL;
    GrB_Matrix A = NULL;
    GrB_Matrix W = NULL;
    GrB_Matrix Miset = NULL;
    GrB_Matrix A_rows = NULL;

    GrB_Index n;
    GrB_Index niset;
    GrB_Index ncols;

    GxB_Container S_container = NULL;
    GxB_Container k_container = NULL;

    GrB_Scalar argmax_0;
    GrB_Type Tuple = NULL;
    GxB_IndexBinaryOp MAKEAMTUP_op = NULL;
    GrB_BinaryOp MAKEAMTUP_Bop = NULL, AM_Bop = NULL;
    GrB_Monoid AM_mon = NULL;
    GrB_Semiring AM_Semiring = NULL;
    GrB_Type Theta_UDT = NULL;

    double *d_copy = NULL;
    uint32_t *c_copy = NULL;

    GRB_TRY(GxB_Type_new(&Theta_UDT, sizeof(Theta), "Theta", THETA_DEFN));
    GRB_TRY(GrB_Scalar_new(&argmax_0, Theta_UDT));

    GRB_TRY(GxB_Type_new(&Tuple, sizeof(argmax_tup), "argmax_tup", AM_TUP));

    GrB_UnaryOp extract_k_if_gain_op = NULL;
    GRB_TRY(GxB_UnaryOp_new(&extract_k_if_gain_op, extract_k_if_gain, GrB_INT64, Tuple, "extract_k_if_gain", EXTRACT_K_IF_GAIN_SRC));

    LG_ASSERT(S_result != NULL, GrB_NULL_POINTER);
    GrB_Info info;
    dbg(A);
    GRB_TRY(GrB_Matrix_dup(&A, G->A));
    dbg(Tuple);

    GRB_TRY(GrB_Matrix_nrows(&n, A));
    GRB_TRY(GrB_Matrix_ncols(&ncols, A));
    // printf("n: %lu\n", n);

    GRB_TRY(GrB_Vector_new(&y, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&x, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_new(&k, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&Wy, Tuple, n));
    GRB_TRY(GrB_Vector_new(&iset, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_new(&gain_values, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&k_values, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&gain_mask, GrB_BOOL, n));

    GRB_TRY(GrB_Matrix_new(&W, GrB_FP64, n, n));
    GRB_TRY(GrB_Descriptor_new(&ri));
    GRB_TRY(GrB_set(ri, GxB_USE_INDICES, GxB_ROWINDEX_LIST));
    GRB_TRY(GrB_Descriptor_new(&rv));
    GRB_TRY(GrB_set(rv, GxB_USE_VALUES, GxB_ROWINDEX_LIST));

    GRB_TRY(GxB_Container_new(&S_container));
    GRB_TRY(GxB_Container_new(&k_container));

    double m = 0.0;

    void *f = NULL;
    uint64_t f_size, f_nvals = 0, f_nheld = 0;
    GrB_Type ftype = NULL;
    int f_handling;

    void *c = NULL;
    uint64_t c_size, c_nvals = 0, c_nheld = 0;
    GrB_Type ctype = NULL;
    int c_handling;

    // S <- I

    GRB_TRY(GrB_assign(x, NULL, NULL, true, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&S, x, 0));
    GRB_TRY(GrB_set(S, GxB_SPARSE, GxB_SPARSITY_CONTROL));
    dbg(S);
    GRB_TRY(GxB_unload_Matrix_into_Container(S, S_container, NULL));
    GRB_TRY(GrB_Vector_dup(&Si_old, S_container->i));
    dbg(Si_old);
    GRB_TRY(GxB_load_Matrix_from_Container(S, S_container, NULL));
    GRB_TRY(GrB_Matrix_reduce_Monoid(k, NULL, NULL, GrB_PLUS_MONOID_FP64, A, NULL));
    GRB_TRY(GrB_set(k, GxB_SPARSE, GxB_SPARSITY_CONTROL));
    dbg(k);
    GRB_TRY(GrB_Vector_reduce_FP64(&m, NULL, GrB_PLUS_MONOID_FP64, k, NULL));
    m /= 2;
    GRB_TRY(GxB_unload_Vector_into_Container(k, k_container, NULL));
    GRB_TRY(GxB_Vector_unload(k_container->x, &f, &ftype, &f_nheld, &f_size, &f_handling, NULL));
    GRB_TRY(GxB_unload_Matrix_into_Container(S, S_container, NULL));
    GRB_TRY(GxB_Vector_unload(S_container->i, &c, &ctype, &c_nheld, &c_size, &c_handling, NULL));

    GRB_TRY(build_argmax_operator(
        Theta_UDT, Tuple, f, m, seed,
        &MAKEAMTUP_Bop, &AM_Bop, &AM_mon, &AM_Semiring,
        &MAKEAMTUP_op, &argmax_0, msg));
    GRB_TRY(GxB_Vector_load(k_container->x, &f, ftype, f_nheld, f_size, f_handling, NULL));
    dbg(k_container->x);
    GRB_TRY(GxB_load_Vector_from_Container(k, k_container, NULL));
    dbg(k);
    GRB_TRY(GxB_Vector_load(S_container->i, &c, ctype, c_nheld, c_size, c_handling, NULL));
    GRB_TRY(GxB_load_Matrix_from_Container(S, S_container, NULL));
    dbg(S);
    bool changed = true;

    int iter = 0;
    double Q = 0;
    double gamma = 1;
    int max_iter = 20;

    #ifdef TIMING
    double tsimple = LAGraph_WallClockTime();
    #endif
    GRB_TRY(LAGraph_IsolateSets(&Miset, A, seed, msg));
    #ifdef TIMING
    tsimple = LAGraph_WallClockTime() - tsimple;
    printf("Isolate Set calc time: %10.10f\n", tsimple);
    #endif

    GrB_Index loop;
    GRB_TRY(GrB_Matrix_nrows(&loop, Miset));

    // allocate an array of GrB_Matrix handles
    GrB_Matrix *A_iset = malloc(loop * sizeof(GrB_Matrix));
    if (A_iset == NULL)
    {
        fprintf(stderr, "Out of memory allocating A_iset\n");
        exit(1);
    }
    GrB_Index n_k_values;
    #ifdef TIMING
    double tsimple = LAGraph_WallClockTime();
    #endif
    for (int i = 0; i < loop; i++)
    {
        // extract column into iset
        GRB_TRY(GrB_Col_extract(iset, NULL, NULL, Miset, GrB_ALL, n, i, GrB_DESC_T0));
        // dbg(iset);

        GrB_Index niset;
        GrB_Vector_nvals(&niset, iset);

        // temp matrix to hold extracted rows
        GrB_Matrix A_rows = NULL;
        GRB_TRY(GrB_Matrix_new(&A_rows, GrB_FP64, niset, ncols));

        // fill A_rows from A
        GRB_TRY(GxB_Matrix_extract_Vector(A_rows, NULL, NULL, A, iset, NULL, ri));
        // dbg(A_rows);

        // now create the ith matrix in A_iset
        GRB_TRY(GrB_Matrix_new(&A_iset[i], GrB_FP64, n, n));

        // assign A_rows into it
        GRB_TRY(GxB_Matrix_assign_Vector(A_iset[i], NULL, NULL, A_rows, iset, NULL, ri));

        GrB_free(&A_rows);
        A_rows = NULL;
    }

    #ifdef TIMING
    tsimple = LAGraph_WallClockTime() - tsimple;
    printf("A_ist array set calc time: %10.10f\n", tsimple);
    double tsimple2 = LAGraph_WallClockTime();
    #endif
    bool recompute = false;

    while (changed && iter < max_iter)
    {
        recompute = changed;
        changed = false;
        for (int i = 0; i < loop; i++)
        {
            double tsimple3 = LAGraph_WallClockTime();

            dbg(S);
            GRB_TRY(GrB_mxm(W, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, A_iset[i], S, NULL));
            dbg(W);
            if(recompute){
                GRB_TRY(GrB_vxm(y, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, k, S, GrB_DESC_T0));
                dbg(y);
            }
            GRB_TRY(GrB_mxv(Wy, NULL, NULL, AM_Semiring, W, y, NULL));
            dbg(Wy);

            GRB_TRY(GxB_unload_Matrix_into_Container(S, S_container, NULL));

            GRB_TRY(GrB_Vector_apply(k_values, NULL, NULL, extract_k_if_gain_op, Wy, NULL));
            GRB_TRY(GrB_Vector_nvals(&n_k_values, k_values));
            if (n_k_values > 0)
            {
                GRB_TRY(GrB_assign(S_container->i, k_values, NULL, k_values, GrB_ALL, n, GrB_DESC_S));
            }

            dbg(S_container->i);
            GRB_TRY(GxB_load_Matrix_from_Container(S, S_container, NULL));
            dbg(S);
        }

        // break;
        GRB_TRY(GxB_unload_Matrix_into_Container(S, S_container, NULL));
        GRB_TRY(LAGraph_Vector_IsEqual(&changed, Si_old, S_container->i, msg));
        changed = !changed;
        if (changed)
        { // skip if they are the same
            GRB_TRY(GrB_assign(Si_old, NULL, NULL, S_container->i, GrB_ALL, n, GrB_DESC_S));
        }


        GRB_TRY(GxB_load_Matrix_from_Container(S, S_container, NULL));
        dbg(Si_old);

        iter++;
    }
    #ifdef TIMING
    tsimple2 = LAGraph_WallClockTime() - tsimple2;
    printf("Main calc time: %10.10f\n", tsimple2);
    printf("Iterations: %d\n", iter);
    #endif

    for (int i = 0; i < loop; i++)
    {
        GrB_free(&A_iset[i]);
    }
    (*S_result) = S;
    S = NULL;

    LG_FREE_EXTRACT;
    LG_FREE_ALL;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED);
#endif	
}
