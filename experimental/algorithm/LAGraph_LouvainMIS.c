#include "LG_internal.h"
#include <LAGraphX.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#define DEBUG 1
#define dbg(x) \
    if (DEBUG) \
    GxB_print(x, 5)
#undef LG_FREE_ALL
#define LG_FREE_ALL                  \
    {                                \
        GrB_free(&x);                \
        GrB_free(&y);                \
        GrB_free(&iset);             \
        GrB_free(&k);                \
        GrB_free(&neighbours);       \
        GrB_free(&S_vector);         \
        GrB_free(&temp);             \
        GrB_free(&temp1);            \
        GrB_free(&S);                \
        GrB_free(&S_container);      \
        GrB_free(&iset_container);   \
        GrB_free(&A_iset_container); \
        GrB_free(&A);                \
        GrB_free(&W);                \
        GrB_free(&Wy);               \
        GrB_free(&A_iset);           \
        GrB_free(&AM_Semiring);      \
        GrB_free(&AM_mon);           \
    }

typedef struct Theta
{
    GrB_Vector* d;
    GrB_Vector* c; //c arrays
    double m;
    uint64_t seed;
} Theta;

typedef struct argmax_tup
{
    double score; // change in modularity
    int64_t k; //who
} argmax_tup;

#define AM_TUP "typedef struct argmax_tup{ double score; int64t k; } argmax_tup;"

void make_argmax_tup(argmax_tup *z,
                     const double *x, GrB_Index ix, GrB_Index jx,
                     const double *y, GrB_Index iy, GrB_Index jy,
                     const void *theta)
{
    Theta *_theta = (Theta *)theta;
    double _d;
    uint64_t _c;

    uint64_t seed = _theta->seed;
    seed ^= seed << 13;
    seed ^= seed >> 7;
    seed ^= seed << 17;
    z->k = (int64_t)jx;
    z->score = (*x) - (_d) / (2 * _theta->m) * (*y);
}
#define MAKE_AM_TUP                                                               \
    "void make_argmax_tup(argmax_tup *z,\n"                                       \
    "                     const double *x, GrB_Index ix, GrB_Index jx,\n"         \
    "                     const double *y, GrB_Index iy, GrB_Index jy,\n"         \
    "                     const void *theta)\n"                                   \
    "{\n"                                                                         \
    "    Theta *_theta = (Theta *)theta;\n"                                       \
    "    double _d;\n"                                                            \
    "    uint64_t _c;\n"                                                          \
    "    GrB_Info info = GrB_Vector_extractElement_UINT64(&_d, _theta->d, ix);\n" \
    "    if (!(info == GrB_SUCCESS || info == GrB_NO_VALUE))\n"                   \
    "    {\n"                                                                     \
    "        char *err;\n"                                                        \
    "        GrB_error(&err, _theta->d);\n"                                       \
    "        printf(\"\\ninfo: %d error: %s\\n\", info, err);\n"                  \
    "    }\n"                                                                     \
    "    info = GrB_Vector_extractElement_UINT64(&_c, _theta->c, ix);\n"          \
    "    if (!(info == GrB_SUCCESS || info == GrB_NO_VALUE))\n"                   \
    "    {\n"                                                                     \
    "        char *err;\n"                                                        \
    "        GrB_error(&err, _theta->c);\n"                                       \
    "        printf(\"\\ninfo: %d error: %s\\n\", info, err);\n"                  \
    "    }\n"                                                                     \
    "    uint64_t seed = (*y + ix + iy + jy);\n"                                  \
    "    seed ^= seed << 13;\n"                                                   \
    "    seed ^= seed >> 7;\n"                                                    \
    "    seed ^= seed << 17;\n"                                                   \
    "    z->tb = seed;\n"                                                         \
    "    z->k = (int64_t)jx;\n"                                                   \
    "    z->score = (*x) - (_d) / (2 * _theta->m) * (*y);\n"                      \
    "}\n"

void argmax_op(argmax_tup *z, argmax_tup *x, argmax_tup *y)
{
    if (x->score > y->score)
    {
        z->score = x->score;
        z->k = x->k;
        z->tb = x->tb;
    }
    else if (x->score == y->score)
    {
        if (x->tb > y->tb)
        {
            z->score = x->score;
            z->k = x->k;
            z->tb = x->tb;
        }
        else
        {
            z->score = y->score;
            z->k = y->k;
            z->tb = y->tb;
        }
    }
    else
    {
        z->score = y->score;
        z->k = y->k;
        z->tb = y->tb;
    }
}
#define AM_OP                                                       \
    "void argmax_op(argmax_tup *z, argmax_tup *x, argmax_tup *y)\n" \
    "{\n"                                                           \
    "    if (x->score > y->score)\n"                                \
    "    {\n"                                                       \
    "        z->score = x->score;\n"                                \
    "        z->k = x->k;\n"                                        \
    "        z->tb = x->tb;\n"                                      \
    "    }\n"                                                       \
    "    else if (x->score == y->score)\n"                          \
    "    {\n"                                                       \
    "        if (x->tb > y->tb)\n"                                  \
    "        {\n"                                                   \
    "            z->score = x->score;\n"                            \
    "            z->k = x->k;\n"                                    \
    "            z->tb = x->tb;\n"                                  \
    "        }\n"                                                   \
    "        else\n"                                                \
    "        {\n"                                                   \
    "            z->score = y->score;\n"                            \
    "            z->k = y->k;\n"                                    \
    "            z->tb = y->tb;\n"                                  \
    "        }\n"                                                   \
    "    }\n"                                                       \
    "    else\n"                                                    \
    "    {\n"                                                       \
    "        z->score = y->score;\n"                                \
    "        z->k = y->k;\n"                                        \
    "        z->tb = y->tb;\n"                                      \
    "    }\n"                                                       \
    "}\n"

int LAGraph_LouvainMIS(
    // output
    GrB_Matrix *S_result,
    // input
    LAGraph_Graph G,
    char *msg)
{
    #if LG_SUITESPARSE_GRAPHBLAS_V10
    LG_CLEAR_MSG;

    char MATRIX_TYPE[LAGRAPH_MSG_LEN];
    char *err;
    if (DEBUG)
        GrB_set(GrB_GLOBAL, true, GxB_BURBLE);
    // Shortened monoids, Binary ops, and Semirings
    GrB_Monoid plusmon = GrB_PLUS_MONOID_FP64;
    GrB_Monoid timesmon = GrB_TIMES_MONOID_FP64;
    GrB_BinaryOp plusf64 = GrB_PLUS_FP64;
    GrB_BinaryOp timesf64 = GrB_TIMES_FP64;
    GrB_BinaryOp divf64 = GrB_DIV_FP64;
    GrB_Semiring stdmxm = GrB_PLUS_TIMES_SEMIRING_FP64;
    GrB_BinaryOp UDT_AM;

    // Declarations
    GrB_Vector iset = NULL;
    GrB_Vector k = NULL;
    GrB_Vector x = NULL;
    GrB_Vector y = NULL;
    GrB_Vector neighbours = NULL;
    GrB_Vector S_vector = NULL;
    GrB_Vector temp = NULL;
    GrB_Vector temp1 = NULL;
    GrB_Matrix S = NULL;
    GxB_Container S_container = NULL;
    GxB_Container iset_container = NULL;
    GxB_Container A_iset_container = NULL;
    GrB_Matrix A = NULL;
    GrB_Matrix A_iset = NULL;
    GrB_Index n;
    GrB_Matrix W = NULL;
    GrB_Vector Wy = NULL;
    Theta theta;
    GrB_Scalar argmax_0;
    GrB_Type Tuple = NULL;
    GxB_IndexBinaryOp MAKEAMTUP_op = NULL;
    GrB_BinaryOp MAKEAMTUP_Bop = NULL, AM_Bop = NULL;
    GrB_Monoid AM_mon = NULL;
    GrB_Semiring AM_Semiring = NULL;

    // GrB_Matrix extractedIset= NULL;

    // Initializing
    LG_TRY(LAGraph_CheckGraph(G, msg));
    LG_ASSERT(S_result != NULL, GrB_NULL_POINTER);

    A = G->A;
    dbg(A);
    uint64_t seed= 1231245; 
    // printf("here");
    // -----------------------------Index Binary OP: AM--------------------------//

    GRB_TRY(GxB_Type_new(&Tuple, sizeof(argmax_tup), "argmax_tup", AM_TUP));

    // -------------------------------------------------------//
    // printf("here0");
    // k = [+_j A(:,j)]
    GRB_TRY(GrB_Matrix_nrows(&n, A));
    GRB_TRY(GrB_Vector_new(&k, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&y, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&x, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&neighbours, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&S_vector, GrB_FP64, n));
    // GRB_TRY(GrB_Vector_new(&temp,GrB_FP64,n));
    // GRB_TRY(GrB_Vector_new(&temp1,GrB_FP64,n));
    GRB_TRY(GrB_Matrix_new(&W, GrB_FP64, n, n));
    GRB_TRY(GrB_Vector_new(&Wy, Tuple, n));
    GRB_TRY(GrB_Matrix_new(&A_iset, GrB_FP64, n, n));
    // GRB_TRY(GrB_Matrix_new(&extractedIset,GrB_FP64,n,n));
    GRB_TRY(GxB_Container_new(&S_container));
    GRB_TRY(GxB_Container_new(&iset_container));
    GRB_TRY(GxB_Container_new(&A_iset_container));
    GRB_TRY(GrB_Matrix_reduce_Monoid(k, NULL, NULL, plusmon, A, NULL));
    dbg(k);
    double m;
    GRB_TRY(GrB_Vector_reduce_FP64(&m, NULL, plusmon, k, NULL));
    m *= 0.5;
    printf("Total edge weight (m): %f\n", m);

    // S <- I
    GRB_TRY(GrB_assign(x, NULL, NULL, 1, GrB_ALL, n, NULL));
    // GxB_print(i,5);
    dbg(x);

    GRB_TRY(GrB_Matrix_diag(&S, x, 0));
    // GRB_TRY(GrB_set(S, false, GxB_ISO));
    GRB_TRY(GrB_set(S, GxB_SPARSE, GxB_SPARSITY_CONTROL));
    dbg(S);
    GRB_TRY(GxB_unload_Matrix_into_Container(S, S_container, NULL));
    dbg(S_container->x);

    GRB_TRY(GxB_load_Matrix_from_Container(S, S_container, NULL));
    dbg(S);
    // Compute Isolate Set
    GRB_TRY(LAGraph_IsolateSets(&iset, G, 32121, msg));
    dbg(iset);

    // Compute max change in modularity for each node in the isolate set
    // GRB_TRY(GrB_set(iset, false, GxB_ISO));
    GRB_TRY(GrB_set(iset, GxB_SPARSE, GxB_SPARSITY_CONTROL));
    // GRB_TRY(GxB_unload_Vector_into_Container(iset,iset_container,NULL));
    dbg(iset);

    // dbg(S)
    // dbg(x);
    GrB_Index niset;
    GrB_Vector_nvals(&niset, iset);
    GrB_Index ncols;
    GrB_Matrix_ncols(&ncols, A);
    GrB_Index *iset_rindices = malloc(niset * sizeof(GrB_Index));
    void *values = malloc(niset * sizeof(double));
    GRB_TRY(GrB_Vector_extractTuples_FP64(iset_rindices, values, &n, iset));
    for (int i = 0; i < niset; i++)
    {
        printf("%ld,", iset_rindices[i]);
    }
    GrB_Matrix A_rows;
    GrB_Matrix_new(&A_rows, GrB_FP64, niset, ncols); // only selected rows

    // Extract rows from A into A_rows
    GrB_Matrix_extract(A_rows, NULL, NULL, A,
                       iset_rindices, niset,
                       GrB_ALL, ncols,
                       NULL);
    GrB_Info info = GxB_subassign(
        A_iset,               // destination (full size)
        NULL,                 // mask
        NULL,                 // accum
        A_rows,               // source (only rows needed)
        iset_rindices, niset, // destination rows
        GrB_ALL, ncols,       // all columns
        NULL                  // descriptor
    );
    dbg(A_iset);
    // GrB_Info info = GrB_Matrix_extract(A_iset, NULL, NULL, A, iset_rindices, n, GrB_ALL, n, NULL);
    if (!(info == GrB_SUCCESS || info == GrB_NO_VALUE))
    {
        char *err;
        GrB_error(&err, A_iset);
        printf("\ninfo: %d error: %s\n", info, err);
    }
    GRB_TRY(GrB_Matrix_wait(A_iset, GrB_MATERIALIZE));
    // GRB_TRY(GxB_Matrix_assign_Vector(A_iset,NULL,NULL,A,temp1,temp, NULL));
    // GRB_TRY(GxB_Matrix_extractTuples_Vector(iset_container->i,iset_container->p,iset_container->x,A));
    // printf("here");
    dbg(A_iset);
    dbg(S);
    GRB_TRY(GrB_mxm(W, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, A_iset, S, NULL));
    dbg(W);
    GRB_TRY(GrB_vxm(y, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, k, S, GrB_DESC_T0));
    dbg(y);

    theta->d = &k;
    // theta.c = k;
    theta.seed = seed;
    theta.m = m;

    GRB_TRY(GrB_Scalar_new(&argmax_0, Theta));
    //theta_scalar = {d = c = }
    GRB_TRY(GrB_Scalar_setElement_UDT(argmax_0, &theta));//(argmax_0,(void *)theta_scalar);

    GRB_TRY(GxB_IndexBinaryOp_new(&MAKEAMTUP_op, (GxB_index_binary_function)make_argmax_tup, Tuple, GrB_FP64, GrB_FP64, GrB_UDT_CODE, "make_argmax_tup", MAKE_AM_TUP));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&MAKEAMTUP_Bop, MAKEAMTUP_op, argmax_0));
    argmax_tup id;
    memset(&id, 0, sizeof(argmax_tup));
    id.k = INT64_MAX;
    id.score = (double)(-INFINITY);
    id.tb = UINT64_MAX;
    GRB_TRY(GxB_BinaryOp_new(&AM_Bop, (GxB_binary_function)argmax_op, Tuple, Tuple, Tuple, "argmax_op", AM_OP));
    GRB_TRY(GrB_Monoid_new_UDT(&AM_mon, AM_Bop, &id));
    GRB_TRY(GrB_Semiring_new(&AM_Semiring, AM_mon, MAKEAMTUP_Bop));

    info = GrB_mxv(Wy, NULL, NULL, AM_Semiring, W, y, NULL);
    if (!(info == GrB_SUCCESS || info == GrB_NO_VALUE))
    {
        char *err;
        GrB_error(&err, Wy);
        printf("\ninfo: %d error: %s\n", info, err);
    }
    dbg(Wy);

    // Aggregate Graph
    // Iterate till no change is detected

    // double Q;
    // double gamma = 1;
    // GRB_TRY(LAGr_Modularity2(&Q,gamma,A,S,msg));
    // // printf("Iterations: %d\n", iter);
    // printf("Q:%.15g\n",Q);
    if (S_result != NULL)
    {
        if (*S_result != NULL)
        {
            GrB_free(S_result);
        }
        *S_result = S;
        S = NULL;
    }
    LG_FREE_ALL;
#elif 

    return 0;

}