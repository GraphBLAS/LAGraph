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

// TDO: This only works with a given Matrix and not a Graph;

// Q = (1/2m)*Trace(S^T*B*S)
// B the modularity matrix = A-(kk^t/2m)

// Given a symmetric graph A with no-self edges, LAGr_Modularity calculates the
// Modularity of that graph and a given community matrix
// 

//Newman ME. Modularity and community structure in networks. 
//Proc Natl Acad Sci U S A. 2006 Jun 6;103(23):8577-82. doi: 10.1073/pnas.0601602103. 
//pub 2006 May 24. PMID: 16723398; PMCID: PMC1482622.

// Current Test File: experimental/test/test_modularity.c

#include "LG_internal.h"
#include <LAGraphX.h>
#include <stdio.h>

#define LG_FREE_MOD                     \
    {                                   \
        GrB_free (&k) ;                 \
        GrB_free (&k_) ;                \
        GrB_free (&kk_) ;               \
        GrB_free (&BS) ;                \
        GrB_free (&S_BS) ;              \
        GrB_free (&Diag) ;              \
    }
#undef LG_FREE_ALL
#define LG_FREE_ALL                     \
    {                                   \
        LG_FREE_MOD;                   \
    }

#define DEBUG 0

int LAGr_Modularity2(
    //output
    double *Q, //modularity Q
    // GrB_Matrix B,
    //input
    double gamma, //Optional resolution limit: default is 1
    // LAGraph_Graph G, TODO
    GrB_Matrix A, //adjacency matrix
    GrB_Matrix S, //community matrix
    char* msg
)
{
    LG_CLEAR_MSG ;

    char MATRIX_TYPE[LAGRAPH_MSG_LEN];
    if (DEBUG)
        GrB_set (GrB_GLOBAL, true, GxB_BURBLE);
    //------------------------------------------------------------------------------
    // Declare Monoids, Brinary Operations, Semirings,(for easier reference) and Matrices
    //------------------------------------------------------------------------------
    GrB_Monoid plusmon = GrB_PLUS_MONOID_FP64;

    GrB_BinaryOp plusf64 = GrB_PLUS_FP64;
    GrB_BinaryOp timesf64 = GrB_TIMES_FP64;

    GrB_Semiring stdmxm = GrB_PLUS_TIMES_SEMIRING_FP64;

    GrB_Index n;
    GrB_Matrix k =    NULL; // vector where the ith value is the degree of vertex i
    GrB_Matrix k_ =   NULL; // transpose of k vector
    GrB_Matrix kk_ =  NULL; // Outer Product of k and transpose(k)
    GrB_Matrix B =    NULL; // B = A - (kk^t/2m)
    GrB_Matrix BS =   NULL; // BS
    GrB_Matrix S_BS = NULL; // S_BS
    GrB_Matrix Diag = NULL;

    GRB_TRY(GrB_Matrix_nrows(&n, A));
    // GRB_TRY(GrB_Matrix_new(&S_, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&B, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&k, GrB_FP64, n,1));
    GRB_TRY(GrB_Matrix_new(&k_, GrB_FP64, n,1));
    GRB_TRY(GrB_Matrix_new(&kk_, GrB_FP64, n,n));

    // GRB_TRY(GrB_transpose(S_,NULL,NULL,S,NULL));
    GRB_TRY (GrB_Matrix_reduce_Monoid ((GrB_Vector)k, NULL, NULL,plusmon, A, NULL));
    // GxB_print(k,3);
    //------------------------------------------------------------------------------
    // Calculation of the adjacency Matrix B = A - kk_/2m
    //------------------------------------------------------------------------------
    double m;
    GRB_TRY(GrB_Matrix_reduce_FP64(&m,plusf64,plusmon,A,NULL));
    m/=2;
    // printf("m:%f\n",m);
    GxB_print(S,5);
    // GxB_print(A,5);
    GRB_TRY(GrB_Matrix_reduce_Monoid ((GrB_Vector)k_,NULL,NULL,plusmon,A, GrB_DESC_T0));
    // GxB_print(k_,3);
    GRB_TRY(GrB_mxm(kk_,NULL,NULL,stdmxm,k,k_,GrB_DESC_T1));
    double inv_m = -gamma/(2*m); 
    GRB_TRY(GrB_Matrix_apply_BinaryOp2nd_FP64(kk_,NULL,NULL,timesf64,kk_,inv_m,NULL));
    GRB_TRY(GrB_eWiseAdd(B,NULL,NULL,plusf64,A,kk_,NULL));
    // GxB_print (B,3); //DENSE

    //------------------------------------------------------------------------------
    // Calculation of S_BS
    //------------------------------------------------------------------------------
    GRB_TRY(GrB_Matrix_new(&BS, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&S_BS, GrB_FP64, n, n));
    GRB_TRY(GrB_mxm(BS, NULL, NULL, stdmxm, B, S, NULL));
    GRB_TRY(GrB_mxm(S_BS,NULL,NULL,stdmxm,S,BS,GrB_DESC_T0));  
    // GxB_print (S_BS,3);


    //------------------------------------------------------------------------------
    // Final Computation of Modularity Q
    //------------------------------------------------------------------------------
    GRB_TRY(GrB_Matrix_new(&Diag,GrB_FP64,n,n));
    GRB_TRY(GrB_select(Diag,NULL,NULL,GrB_DIAG,S_BS,0,NULL));
    // GxB_print(Diag,5);
    double Q_;
    printf("here");
    GRB_TRY(GrB_Matrix_reduce_FP64(&Q_,NULL,plusmon,Diag,NULL));
    Q_ *= -inv_m;
    *Q = Q_;
    LG_FREE_ALL;
    return 0;
}
