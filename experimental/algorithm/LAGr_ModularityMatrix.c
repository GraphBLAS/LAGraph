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

#include "GraphBLAS.h"
#include "LG_internal.h"
#include <LAGraphX.h>

#undef LG_FREE_WORK
#define LG_FREE_WORK      \
    {                     \
        GrB_free (&k);    \
        GrB_free (&B);    \
        GrB_free (&S_t);  \
    }
#undef LG_FREE_ALL
#define LG_FREE_ALL  \
    {                \
        LG_FREE_WORK; \
    }

int LAGr_AdjModularity (
    double *Q,
    double gamma,
    const GrB_Matrix A,
    const GrB_Matrix S,
    char *msg)
{
#if LG_SUITESPARSE_GRAPHBLAS_V10_2
    LG_CLEAR_MSG;

    GrB_Index nrows, ncols, nrows_s, ncols_s;
    GrB_Vector k = NULL, com_deg = NULL;
    GrB_Matrix B = NULL;
    GrB_Matrix S_t = NULL;

    LG_ASSERT (A != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (S != NULL, GrB_NULL_POINTER) ;

    GxB_print(A, 2);
    GRB_TRY (GrB_Matrix_nrows (&nrows, A));
    GRB_TRY (GrB_Matrix_ncols (&ncols, A));
    GRB_TRY (GrB_Matrix_nrows (&nrows_s, S));
    GRB_TRY (GrB_Matrix_ncols (&ncols_s, S));

    LG_ASSERT (nrows == ncols,   GrB_DIMENSION_MISMATCH) ;
    LG_ASSERT (ncols == nrows_s, GrB_DIMENSION_MISMATCH) ;

    GRB_TRY (GrB_Matrix_new (&B, GrB_FP64, nrows, ncols));
    GRB_TRY (GrB_Matrix_new (&S_t, GrB_BOOL, ncols_s, nrows_s));
    GRB_TRY (GrB_Vector_new (&k, GrB_FP64, ncols));
    GRB_TRY (GrB_Vector_new (&com_deg, GrB_FP64, ncols_s));

    double m = 0.0;
    double Q_ = 0.0;

    GRB_TRY (GrB_reduce (k, NULL, NULL, GrB_PLUS_MONOID_FP64, A, NULL));
    GRB_TRY (GrB_reduce (&m, NULL, GrB_PLUS_MONOID_FP64, k, NULL));

    m /= 2.0;
    if (m == 0.0)
    {
        *Q = 0.0;
        LG_FREE_WORK;
        return GrB_SUCCESS;
    }
    double inv_m = -gamma / (2.0 * m);
    // sum degrees per community
    GRB_TRY (GrB_transpose (S_t, NULL, NULL, S, NULL)) ;
    GRB_TRY (GrB_mxv (com_deg, NULL, NULL, GxB_PLUS_SECOND_FP64, S_t, k, NULL)) ;
    // square
    GRB_TRY (GrB_apply (com_deg, NULL, NULL, GxB_POW_FP64, com_deg, 2.0, NULL)) ;

    // sum all of the products
    // this computed \sum_{i,j} (k_{i}k_{j}\delta(\sigma_i\sigma_j))
    GRB_TRY (GrB_reduce (&Q_, NULL, GrB_PLUS_MONOID_FP64, com_deg, NULL)) ;
    Q_ *= inv_m;
    GRB_TRY (GrB_free (&com_deg)) ;

    // Count the number of edges per node in the cluster originating at a node
    // in the cluster
    GRB_TRY (GrB_mxm (B, S, NULL, GxB_PLUS_FIRST_FP64, A, S_t, GrB_DESC_ST1));

    // sum the intra-cluster edges per node for all nodes
    GRB_TRY (GrB_reduce (&Q_, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, B, NULL)) ;

    Q_ *= -inv_m;
    *Q = Q_;

    LG_FREE_WORK;
    return (GrB_SUCCESS);
#else
    return (GrB_NOT_IMPLEMENTED);
#endif
}
