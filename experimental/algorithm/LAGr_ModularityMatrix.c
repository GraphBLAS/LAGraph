//------------------------------------------------------------------------------
// LAGr_Modularity2.c: Calculates the modularity of a graph
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
// Q = (1/2m)*Trace(S^T*B*S)
// B the modularity matrix = A-(kk^t/2m)

// Given a symmetric graph A with no-self edges, LAGr_Modularity calculates the
// Modularity of that graph and a given community matrix
//

// Newman ME. Modularity and community structure in networks.
// Proc Natl Acad Sci U S A. 2006 Jun 6;103(23):8577-82. doi: 10.1073/pnas.0601602103.
// pub 2006 May 24. PMID: 16723398; PMCID: PMC1482622.

// Current Test File: experimental/test/test_modularity.c

#include "LG_internal.h"
#include <LAGraphX.h>
#include <stdio.h>

#define LG_FREE_MOD      \
    {                    \
        GrB_free(&k);    \
        GrB_free(&kk_);  \
        GrB_free(&B);    \
        GrB_free(&BS);   \
        GrB_free(&S_BS); \
        GrB_free(&Diag); \
    }
#undef LG_FREE_ALL
#define LG_FREE_ALL  \
    {                \
        LG_FREE_MOD; \
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

int LAGr_ModularityMatrix(
    double *Q,
    double gamma,
    GrB_Matrix A,
    GrB_Matrix S,
    char *msg)
{
#if LG_SUITESPARSE_GRAPHBLAS_V10_2
    LG_CLEAR_MSG;
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];

    // GrB_set(GrB_GLOBAL, false, GxB_BURBLE);

    GrB_Index n;
    GrB_Matrix k = NULL;
    GrB_Matrix kk_ = NULL;
    GrB_Matrix B = NULL;
    GrB_Matrix BS = NULL;
    GrB_Matrix S_BS = NULL;
    GrB_Matrix Diag = NULL;
    GrB_Matrix mask = NULL;

    GRB_TRY(GrB_Matrix_nrows(&n, A));

    GRB_TRY(GrB_Matrix_new(&B, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&k, GrB_FP64, n, 1));
    GRB_TRY(GrB_Matrix_new(&kk_, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&BS, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&S_BS, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&Diag, GrB_FP64, n, n));

    double m = 0.0;
    double Q_ = 0.0;

    GRB_TRY(GrB_Matrix_reduce_Monoid((GrB_Vector)k, NULL, NULL, GrB_PLUS_MONOID_FP64, A, NULL));

    GRB_TRY(GrB_Matrix_reduce_FP64(&m, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, A, NULL));
    m /= 2.0;
    if (m == 0.0)
    {
        *Q = 0.0;
        LG_FREE_ALL;
        return GrB_SUCCESS;
    }
    GRB_TRY(GrB_mxm(kk_, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, k, k, GrB_DESC_T1));
    double inv_m = -gamma / (2.0 * m);
    GRB_TRY(GrB_Matrix_apply_BinaryOp2nd_FP64(kk_, NULL, NULL, GrB_TIMES_FP64, kk_, inv_m, GrB_DESC_R));
    GRB_TRY(GrB_eWiseAdd(B, NULL, NULL, GrB_PLUS_FP64, A, kk_, NULL));
    GRB_TRY(GrB_mxm(BS, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, B, S, NULL));
    GRB_TRY(GrB_mxm(S_BS, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, S, BS, GrB_DESC_T0));

    GRB_TRY(GrB_select(Diag, NULL, NULL, GrB_DIAG, S_BS, 0, NULL));

    GRB_TRY(GrB_Matrix_reduce_FP64(&Q_, NULL, GrB_PLUS_MONOID_FP64, Diag, NULL));
    Q_ *= -inv_m;
    *Q = Q_;

    LG_FREE_ALL;
    return (GrB_SUCCESS);
#else
    return (GrB_NOT_IMPLEMENTED);
#endif
}
