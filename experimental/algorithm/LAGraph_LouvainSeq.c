//------------------------------------------------------------------------------
// LAGraph_LouvainSeq.c: Runs the Louvain Algorithm on a given graph
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
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#undef LG_FREE_ALL
#define LG_FREE_ALL                           \
    {                                         \
        GrB_free(&k);                         \
        GrB_free(&x);                         \
        GrB_free(&v);                         \
        GrB_free(&srxq);                      \
        GrB_free(&sr);                        \
        GrB_free(&q1);                        \
        GrB_free(&t);                         \
        GrB_free(&t_q);                       \
        GrB_free(&dS);                        \
        GrB_free(&dSTk);                      \
        GrB_free(&vtS);                       \
        GrB_free(&temp);                      \
        GrB_free(&y_rand);                    \
        GrB_free(&max_q1);                    \
        GrB_free(&za);                        \
        GrB_free(&z_dSTk);                    \
        GrB_free(&AS);                        \
        GrB_free(&StAS);                      \
        GrB_free(&z);                         \
        GrB_Scalar_free(&Theta);              \
        GrB_Type_free(&Tuple);                \
        GrB_Semiring_free(&IBop_MAX);         \
        GrB_BinaryOp_free(&MAKEFP64_Bop);     \
        GxB_IndexBinaryOp_free(&MAKEFP64_op); \
        GrB_BinaryOp_free(&MAXFP64_op);       \
        GrB_Monoid_free(&MAXFP64_mon);        \
        GxB_Container_free(&S_container);     \
    }
#define DEBUG 0
#if DEBUG
#define check() printf("here")
#define dbg(x) \
    if (DEBUG) \
    GxB_print(x, 5)
#define err(x, info)                                    \
    if (!(info == GrB_SUCCESS || info == GrB_NO_VALUE)) \
    {                                                   \
        char **err;                                     \
        GrB_error(err, x);                              \
        printf("\ninfo: %d error: %s\n", info, err);    \
    }
#else
#define check()
#define dbg(x)
#define err(x, info)
#endif
typedef struct tuple_fp64
{
    int64_t k;
    double v;
    uint64_t tb;
} tuple_fp64;
#define FP64_K "typedef struct tuple_fp64 { int64_t k ; double v ; uint64_t  tb ;} tuple_fp64 ;"
void make_fp64(tuple_fp64 *z,
               const double *x, GrB_Index ix, GrB_Index jx,
               const uint64_t *y, GrB_Index iy, GrB_Index jy,
               const void *theta)
{
    z->k = (int64_t)jx;
    z->v = (*x);
    uint64_t seed = (*y + ix + iy + jy);
    seed ^= seed << 13;
    seed ^= seed >> 7;
    seed ^= seed << 17;
    z->tb = seed;
}
void max_fp64(tuple_fp64 *z, const tuple_fp64 *x, const tuple_fp64 *y)
{

    if (x->v > y->v)
    {
        z->k = x->k;
        z->v = x->v;
    }
    else if (x->v == y->v)
    {
        if (x->tb > y->tb)
        {
            z->k = y->k;
            z->v = y->v;
        }
        else
        {
            z->k = x->k;
            z->v = x->v;
        }
    }
    else
    {
        z->k = y->k;
        z->v = y->v;
    }
}
#define MAX_FP64                                                                             \
    "void max_fp64(tuple_fp64 *z, const tuple_fp64 *x, const tuple_fp64 *y){ \n"             \
    "   if (x->v > y->v)                                                     \n"             \
    "   {                                                                    \n"             \
    "       z->k = x->k;                                                     \n"             \
    "       z->v = x->v;                                                     \n"             \
    "   }else if(x->v == y->v){ \n"                                                          \
    "         if(x->tb > y->tb){z->k = y->k;z->v = y->v;} else {z->k = x->k;z->v = x->v;}\n" \
    "    }else{                                                               \n"            \
    "       z->k = y->k;                                                     \n"             \
    "       z->v = y->v;                                                     \n"             \
    "    }                                                                   \n"             \
    "}"
#define MAKE_FP64                                                         \
    "void make_fp64(tuple_fp64 *z,                               \n"      \
    "               const double *x, GrB_Index ix, GrB_Index jx, \n"      \
    "               const uint64_t *y, GrB_Index iy, GrB_Index jy,    \n" \
    "               const void *theta)                           \n"      \
    "{                                                           \n"      \
    "    z->k = (int64_t)jx;                                     \n"      \
    "    z->v = (*x);                   \n"                               \
"    uint64_t seed = (*y + ix + iy + jy); \n"                             \
    "    seed ^= seed << 13 ; \n"                                         \
    "    seed ^= seed >> 7 ;\n"                                           \
    "    seed ^= seed << 17 ;\n"                                          \
    "    z->tb = seed;      \n"                                           \
    "}"


int LAGraph_LouvainSeq(
    // output
    GrB_Matrix *S_result, 
    // input
    LAGraph_Graph G,
    uint64_t seed,
    char *msg)
{
#if LG_SUITESPARSE_GRAPHBLAS_V10
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];
    GrB_set(GrB_GLOBAL, false, GxB_BURBLE);


    GrB_Vector t_q = NULL, q1 = NULL, t = NULL, v = NULL;
    GrB_Vector k = NULL;
    GrB_Vector x = NULL;
    GrB_Vector z = NULL;
    GrB_Vector za = NULL;
    GrB_Vector z_dSTk = NULL;
    GrB_Index n, b;
    GrB_Matrix S = NULL;
    GrB_Matrix dS = NULL;
    GrB_Vector sr = NULL;
    GrB_Vector srxq = NULL;
    GrB_Index vals_srxq;
    // GrB_Matrix dS = NULL ;
    GrB_Vector dSTk = NULL, vtS = NULL;
    GrB_Vector temp = NULL;
    GrB_Vector y_rand = NULL;
    GrB_Vector max_q1 = NULL;
    GrB_Matrix AS = NULL;
    GrB_Matrix StAS = NULL;
    GrB_Semiring IBop_MAX = NULL;
    GxB_IndexBinaryOp MAKEFP64_op = NULL;
    GrB_BinaryOp MAKEFP64_Bop = NULL, MAXFP64_op = NULL;
    GrB_Scalar Theta = NULL;
    GrB_Type Tuple = NULL;
    GrB_Monoid MAXFP64_mon = NULL;

    GxB_Container S_container = NULL;
    GrB_Index q1_size;
    // GrB_Container
    tuple_fp64 o;
    double o1;
    double k_i = 0;
    LG_ASSERT(S_result != NULL, GrB_NULL_POINTER);

    GrB_Matrix A = G->A;
    // GxB_print(A,5);
    // printf("rand init");
    double test = 23;
    // index bin op definitions
    //------------------------------------------------------------------------------------------------------------------

    GRB_TRY(GrB_Scalar_new(&Theta, GrB_BOOL));
    GRB_TRY(GrB_Scalar_setElement_BOOL(Theta, 0));
    GRB_TRY(GxB_Type_new(&Tuple, sizeof(tuple_fp64), "tuple_fp64", FP64_K));
    GRB_TRY(GxB_IndexBinaryOp_new(&MAKEFP64_op, (GxB_index_binary_function)make_fp64, Tuple, GrB_FP64, GrB_UINT64, GrB_BOOL, "make_fp64", MAKE_FP64));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&MAKEFP64_Bop, MAKEFP64_op, Theta));
    tuple_fp64 id;
    memset(&id, 0, sizeof(tuple_fp64));
    id.k = INT64_MAX;
    id.v = (double)(-INFINITY);
    GRB_TRY(GxB_BinaryOp_new(&MAXFP64_op, (GxB_binary_function)max_fp64, Tuple, Tuple, Tuple, "max_fp64", MAX_FP64));
    GRB_TRY(GrB_Monoid_new_UDT(&MAXFP64_mon, MAXFP64_op, &id));
    GRB_TRY(GrB_Semiring_new(&IBop_MAX, MAXFP64_mon, MAKEFP64_Bop));
    //------------------------------------------------------------------------------------------------------------------

    GRB_TRY(GrB_Matrix_nrows(&n, A));
    GRB_TRY(GrB_Matrix_ncols(&b, A));

    // k = [+_j A(:,j)]
    GRB_TRY(GrB_Vector_new(&k, GrB_FP64, n));
    GRB_TRY(GrB_Matrix_reduce_Monoid(k, NULL, NULL, GrB_PLUS_MONOID_FP64, A, NULL));
    // GxB_print(k,5);

    // m = .5[+_i k(i)]
    double m;
    GRB_TRY(GrB_Vector_reduce_FP64(&m, NULL, GrB_PLUS_MONOID_FP64, k, NULL));
    m *= .5;
    if (m == 0)
    {
        *S_result = NULL;
        LG_FREE_ALL;
        return 0;
    }
    // printf("m= %f\n", m);

    // S <- I
    GRB_TRY(GrB_Vector_new(&x, GrB_BOOL, n));
    GRB_TRY(GrB_assign(x, NULL, NULL, true, GrB_ALL, n, NULL));
    // GxB_print(i,5);
    GRB_TRY(GrB_Matrix_diag(&S, x, 0));
    GRB_TRY(GrB_set(S, false, GxB_ISO));
    GrB_set(S, GxB_SPARSE, GxB_SPARSITY_CONTROL);
    // GxB_print(S,5);
    GRB_TRY(GrB_Vector_new(&v, GrB_FP64, n));
    // var used in for loop
    GRB_TRY(GrB_Vector_new(&t_q, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&q1, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&z, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&za, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&dSTk, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&vtS, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&z_dSTk, GrB_FP64, n));
    // temp used to  set dS to 0 matrix
    GRB_TRY(GrB_Vector_new(&temp, GrB_FP64, n));
    GRB_TRY(GrB_assign(temp, NULL, NULL, 0, GrB_ALL, n, NULL));
    GRB_TRY(GxB_Container_new(&S_container));
    GRB_TRY(GrB_Vector_new(&max_q1, Tuple, 1));
    GrB_Vector_new(&srxq, GrB_FP64, n);
    GRB_TRY(GrB_Matrix_diag(&dS, temp, 0));
    GRB_TRY(GrB_Vector_new(&y_rand, GrB_UINT64, n));
    GRB_TRY(GrB_Vector_new(&sr, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&srxq, GrB_FP64, n));
    GRB_TRY(GrB_Matrix_new(&AS, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&StAS, GrB_FP64, n, n));
    GRB_TRY(GrB_assign(y_rand, NULL, NULL, 1, GrB_ALL, n, NULL));

    bool changed = true;
    int max_iter = 20;
    int iter = 0;
    int aggr_iter = 0;
    // uint64_t seed = 1212312224;

    // GxB_print(y_rand,5);
    GRB_TRY(GrB_mxv(z, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, S, k, NULL));
    // GxB_print(z,5);
    while (aggr_iter < max_iter)
    {
        while (changed && iter < max_iter)
        {
            changed = false;
            for (int i = 0; i < n; i++)
            { // extract tuples
                // printf("%d",i);

                // v = A(i,:)
                GRB_TRY(GrB_Col_extract(v, NULL, NULL, A, GrB_ALL, b, i, GrB_DESC_T0));
                // GxB_print(v,5);

                // -- extract k_i
                GRB_TRY(GrB_Vector_extractElement_FP64(&k_i, k, i));
                // GxB_print(S,5);

                // t_q =v any.pair S   O(|v|)
                GRB_TRY(GrB_vxm(t_q, NULL, NULL, GxB_ANY_PAIR_FP64, v, S, GrB_DESC_T0));
                // GxB_print(t_q,5);

                // sr = S(i,:)
                GRB_TRY(GrB_Col_extract(sr, NULL, NULL, S, GrB_ALL, 1, i, GrB_DESC_T0));
                // GxB_print(sr,5);

                // S(i,:) = 0/false
                GRB_TRY(GxB_unload_Matrix_into_Container(S, S_container, NULL));
                GRB_TRY(GrB_Vector_setElement_BOOL(S_container->x, false, i));
                GRB_TRY(GxB_load_Matrix_from_Container(S, S_container, NULL));
                dbg(S);

                ////////////////////////////////////////////////////////////
                //-------------q1<t_q> = a(kTS)+vTS----------- -----------//
                //-------------q1<t_q> = z+vTS----------------------------//
                // double alpha_p = 1;
                double alpha = -k_i / m;

                GRB_TRY(GrB_Vector_apply_BinaryOp2nd_FP64(za, NULL, NULL, GrB_TIMES_FP64, z, alpha, GrB_DESC_T0));
                dbg(za);

                GRB_TRY(GrB_eWiseAdd(za, NULL, NULL, GrB_PLUS_FP64, za, v, NULL));
                // printf("z*alpha + v\n");
                dbg(za);
                GRB_TRY(GrB_vxm(q1, t_q, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, za, S, GrB_DESC_R));
                dbg(q1);
                ///////////////////////////////////////////////////////////

                ///////////////////////////////////////////////////////////
                //-------------Index Binary OP Rand_argminmax -----------//
                GRB_TRY(GrB_Vector_nvals(&q1_size, q1));
                // GxB_print(q1,5);
                // printf("Size of q1: %ld\n",q1_size);
                // GRB_TRY(GrB_Vector_setElement_UINT64(y_rand,seed,0));
                seed += 3;
                GRB_TRY(GrB_assign(y_rand, t_q, NULL, seed, GrB_ALL, n, GrB_DESC_S));

                // GxB_print(q1,5);
                GRB_TRY(GrB_mxv(max_q1, NULL, NULL, IBop_MAX, (GrB_Matrix)q1, y_rand, GrB_DESC_T0));
                // GxB_print(q1,5);

                GRB_TRY(GrB_Vector_extractElement_UDT((void *)&o, max_q1, 0));
                // printf("choice:%ld tb: %ld\n",(long)o.k, (long)o.tb);
                // dbg(S);

                GRB_TRY(GxB_unload_Matrix_into_Container(S, S_container, NULL));
                GRB_TRY(GrB_Vector_setElement(S_container->i, o.k, i));
                GRB_TRY(GrB_Vector_setElement_BOOL(S_container->x, true, i));
                GRB_TRY(GxB_load_Matrix_from_Container(S, S_container, NULL));
                // GxB_print(S,5);

                //////////////////////////////////////////////////////////

                // GxB_print(sr,5);
                GRB_TRY(GrB_Vector_eWiseMult_BinaryOp(srxq, NULL, NULL, GrB_TIMES_FP64, sr, q1, NULL));
                GRB_TRY(GrB_Vector_nvals(&vals_srxq, srxq));
                // GxB_print(srxq,5);

                if (vals_srxq == 0)
                {
                    changed = true;
                }
                // if(i==3)break;
            }
            iter++;
            // break;

        }

        aggr_iter++;
    }
    // printf("Iterations: %d\n", iter);
    (*S_result) = S;
    S = NULL;
    LG_FREE_ALL;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED);
#endif
