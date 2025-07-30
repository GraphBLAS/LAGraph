#include "LG_internal.h"
#include <LAGraphX.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#define DEBUG 1
#define dbg(x) \
    if (DEBUG) \
    GxB_print(x, 5)
#define err(x, info)                                    \
    if (!(info == GrB_SUCCESS || info == GrB_NO_VALUE)) \
    {                                                   \
        char **err;                                      \
        GrB_error(err, x);                             \
        printf("\ninfo: %lu error: %s\n", info, err);   \
    }
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
    double *d;
    uint32_t *c; /* c arrays */
    double m;
    uint64_t seed;
} Theta;

#define THETA_DEFN                      \
"typedef struct Theta"                  \
"{"                                     \
"    double *d;"                      \
"    uint32_t *c; /* c arrays */"       \
"    double m;"                         \
"    uint64_t seed;"                    \
"} Theta;"

typedef struct argmax_tup
{
    double score; /* change in modularity */
    int64_t k;    /* who */
    double tb;
} argmax_tup;

#define AM_TUP                                      \
"typedef struct argmax_tup"                         \
"{"                                                 \
"    double score; /* change in modularity */"      \
"    int64_t k;    /* who */"                       \
"    double tb;"                                    \
"} argmax_tup;"

void make_argmax_tup(argmax_tup *z,
                     const double *x, GrB_Index ix, GrB_Index jx,
                     const double *y, GrB_Index iy, GrB_Index jy,
                     const void *theta)
{
    printf("x: %f,%lu,%lu\n", *x, ix, jx);
    printf("y: %f,%lu,%lu\n", *y, iy, jy);
    Theta *_theta = (Theta *)theta;
    uint64_t seed = _theta->seed + (*y + ix + iy + jy);
    seed ^= seed << 13;
    seed ^= seed >> 7;
    seed ^= seed << 17;
    z->k = (int64_t)jx;
    printf("%f or %f\n", _theta->d[iy], _theta->d[ix]);
    z->score = (*x) - ((_theta->d[ix]) * (*y)) / (2*_theta->m);
    z->tb = seed;
    printf("z: %f,%lu,%f\n\n\n", z->score, z->k, z->tb);
}
#define MAKE_AM_TUP                                                       \
    "void make_argmax_tup(argmax_tup *z,\n"                               \
    "                     const double *x, GrB_Index ix, GrB_Index jx,\n" \
    "                     const double *y, GrB_Index iy, GrB_Index jy,\n" \
    "                     const void *theta)\n"                           \
    "{\n"                                                                 \
    "    printf(\"x: %f,%lu,%lu\\n\", *x, ix, jx);\n"                     \
    "    printf(\"y: %f,%lu,%lu\\n\", *y, iy, jy);\n"                     \
    "    Theta *_theta = (Theta*)theta;\n"                                \
    "    uint64_t seed = _theta->seed + (*y + ix + iy + jy);\n"           \
    "    seed ^= seed << 13;\n"                                           \
    "    seed ^= seed >> 7;\n"                                            \
    "    seed ^= seed << 17;\n"                                           \
    "    z->k = (int64_t)jx;\n"                                           \
    "    printf(\"%d or %d\\n\", _theta->d[iy], _theta->d[ix]);\n"        \
    "    z->score = (*x) - ((_theta->d[jx]) * (*y)) / (2 * _theta->m);\n" \
    "    z->tb = seed;\n"                                                 \
    "    printf(\"z: %f,%lu,%f\\n\\n\\n\", z->score, z->k, z->tb);\n"     \
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
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];
    char *err;
    // if (DEBUG)
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
    GrB_Scalar argmax_0;
    GrB_Type Tuple = NULL;
    GxB_IndexBinaryOp MAKEAMTUP_op = NULL;
    GrB_BinaryOp MAKEAMTUP_Bop = NULL, AM_Bop = NULL;
    GrB_Monoid AM_mon = NULL;
    GrB_Semiring AM_Semiring = NULL;

    // Initializing
    LG_TRY(LAGraph_CheckGraph(G, msg));
    LG_ASSERT(S_result != NULL, GrB_NULL_POINTER);

    A = G->A;
    dbg(A);
    uint64_t seed = 1231245;
    // printf("here");
    // -----------------------------Index Binary OP: AM--------------------------//

    GRB_TRY(GxB_Type_new(&Tuple, sizeof(argmax_tup), "argmax_tup", AM_TUP));

    // -------------------------------------------------------//

    GRB_TRY(GrB_Matrix_nrows(&n, A));
    GRB_TRY(GrB_Vector_new(&y, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&x, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&neighbours, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&S_vector, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&k, GrB_FP64, n));
    GRB_TRY(GrB_Matrix_new(&W, GrB_FP64, n, n));
    GRB_TRY(GrB_Vector_new(&Wy, Tuple, n));
    GRB_TRY(GrB_Matrix_new(&A_iset, GrB_FP64, n, n));
    GRB_TRY(GxB_Container_new(&S_container));
    GRB_TRY(GxB_Container_new(&iset_container));
    GRB_TRY(GxB_Container_new(&A_iset_container));

    // k = [+_j A(:,j)]

    GRB_TRY(GrB_Matrix_reduce_Monoid(k, NULL, NULL, plusmon, A, NULL));
    GRB_TRY(GrB_set(k, GxB_SPARSE, GxB_SPARSITY_CONTROL));
    dbg(k);

    double m;
    GRB_TRY(GrB_Vector_reduce_FP64(&m, NULL, plusmon, k, NULL));
    m *= 0.5;
    printf("Total edge weight (m): %f\n", m);

    // S <- I
    GRB_TRY(GrB_assign(x, NULL, NULL, 1, GrB_ALL, n, NULL));
    dbg(x);
    GRB_TRY(GrB_Matrix_diag(&S, x, 0));
    GRB_TRY(GrB_set(S, GxB_SPARSE, GxB_SPARSITY_CONTROL));
    dbg(S);

    GRB_TRY(GxB_unload_Matrix_into_Container(S, S_container, NULL));
    dbg(S_container->x);

    GRB_TRY(GxB_load_Matrix_from_Container(S, S_container, NULL));
    // dbg(x);
    GRB_TRY(LAGraph_IsolateSets(&iset, G, seed, msg));
    dbg(iset);
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
    // GRB_TRY(GrB_Matrix_extra)
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
    GRB_TRY(GrB_mxm(W, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, A_iset, S, NULL));
    dbg(W);
    GRB_TRY(GrB_vxm(y, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, k, S, GrB_DESC_T0));
    dbg(y);
    dbg(k);

    GxB_Container k_container = NULL;
    GRB_TRY(GxB_Container_new(&k_container));
    GRB_TRY(GxB_unload_Vector_into_Container(k, k_container, NULL));
    dbg(k_container->p);
    dbg(k_container->i);
    dbg(k_container->x);
    void *f = NULL;
    uint64_t nvals = 0, nheld = 0;
    GrB_Type xtype = NULL;
    int x_handling;
    uint64_t x_size;
    info = GxB_Vector_unload(k_container->x, &f, &xtype, &nheld, &x_size, &x_handling, NULL); 
    // GxB_load_Vector_from_Container(k,k_container,NULL);   
    // dbg(k);
    dbg(xtype);
    printf("\nsize of f%ld\n",nheld);

    GrB_Type Theta_UDT = NULL;
    GRB_TRY(GxB_Type_new(&Theta_UDT, sizeof(Theta), "Theta", THETA_DEFN));
    GRB_TRY(GrB_Scalar_new(&argmax_0, Theta_UDT));
    Theta theta_scalar;
    theta_scalar.d = (double*)f;
    for(int i =0;i<nheld;i++){
        printf("d[%d]:%f\n", i, theta_scalar.d[i]);
    }
    theta_scalar.c = NULL;
    theta_scalar.m = m;
    theta_scalar.seed = seed;
    // GRB_TRY(GrB_Scalar_setElement_UDT(argmax_0, (void *)&theta_scalar));
    info = GrB_Scalar_setElement_UDT(argmax_0, (void *)&theta_scalar);
    printf("info: %lu",info);
    GRB_TRY(GxB_IndexBinaryOp_new(&MAKEAMTUP_op, (GxB_index_binary_function)make_argmax_tup, Tuple, GrB_FP64, GrB_FP64, Theta_UDT, "make_argmax_tup", MAKE_AM_TUP));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&MAKEAMTUP_Bop, MAKEAMTUP_op, argmax_0));
    argmax_tup id;
    memset(&id, 0, sizeof(argmax_tup));
    id.k = INT64_MAX;
    id.score = (double)(-INFINITY);
    GRB_TRY(GxB_BinaryOp_new(&AM_Bop, (GxB_binary_function)argmax_op, Tuple, Tuple, Tuple, "argmax_op", AM_OP));
    GRB_TRY(GrB_Monoid_new_UDT(&AM_mon, AM_Bop, &id));
    GRB_TRY(GrB_Semiring_new(&AM_Semiring, AM_mon, MAKEAMTUP_Bop));

    info = GrB_mxv(Wy, NULL, NULL, AM_Semiring, W, y, NULL);
    dbg(Wy);
    // err(Wy, info);
    argmax_tup test;
    argmax_tup test2;
    GRB_TRY(GrB_Vector_extractElement_UDT((void*)&test, Wy, 17));
    printf("test score: %f, k: %ld, tb: %f\n", test.score, test.k, test.tb);

    GRB_TRY(GrB_Vector_extractElement_UDT((void*)&test2, Wy, 20));
    printf("test score: %f, k: %ld, tb: %f\n", test2.score, test2.k, test2.tb);

    // Aggregate Graph
    // Iterate till no change is detected

    // double Q;
    // double gamma = 1;
    // GRB_TRY(LAGr_Modularity2(&Q,gamma,A,S,msg));
    // // printf("Iterations: %lu\n", iter);
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
#else
    LG_ASSERT(false, GrB_NOT_IMPLEMENTED);
#endif
    return 0;
}
