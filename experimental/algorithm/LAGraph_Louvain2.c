//------------------------------------------------------------------------------
// LAGraph_Louvain2.c: Runs the Louvain Algorithm on a given graph(under construction)
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


// Current Test File: experimental/test/test_louvain2.c

#include "LG_internal.h"
#include <LAGraphX.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
// #include <iostream>
#undef LG_FREE_ALL
#define LG_FREE_ALL                     \
    {                                   \
        GrB_free (&S) ;                 \
        GrB_free (&k) ;                 \
        GrB_free (&x) ;                 \
        GrB_free (&v) ;                 \
        GrB_free (&srxq) ;              \
        GrB_free (&sr) ;                \
        GrB_free (&q1) ;                \
        GrB_free (&t) ;                 \
        GrB_free (&t_q) ;               \
        GrB_free (&Theta) ;             \
        GrB_free (&Semiring) ;          \
        GrB_free (&Mon) ;               \
        GrB_free (&Tuple) ;             \
        GrB_free (&dS) ;                \
        GrB_free (&dSk) ;               \
        GrB_free (&vtS) ;               \
        GrB_free (&temp) ;              \
        GrB_free (&y_rand) ;            \
        GrB_free (&max_q1) ;            \
        LAGraph_Free ((void *) &Sp, NULL) ;     \
        LAGraph_Free ((void *) &Sj, NULL) ;     \
        LAGraph_Free ((void *) &dSp, NULL) ;     \
        LAGraph_Free ((void *) &dSj, NULL) ;     \
    }
#define DEBUG 0
// uint64_t seed = 213;
typedef struct tuple_fp64{
    int64_t k;
    double v;
}tuple_fp64;
#define FP64_K "typedef struct tuple_fp64 { int64_t k ; double v ; } tuple_fp64 ;"
void make_fp64(tuple_fp64 *z,
               const double *x, GrB_Index ix, GrB_Index jx,
               const uint64_t *y, GrB_Index iy, GrB_Index jy,
               const void *theta)
{
    z->k = (int64_t)jx;
    z->v = (*x) + (((*y) ^ (jy) << iy));
    printf("z->k:%ld, %f\t",z->k,z->v);
}
void max_fp64(tuple_fp64 *z, const tuple_fp64 *x, const tuple_fp64 *y){
    if (x->v > y->v ){
        z->k = x->k;
        z->v = x->v;
    }else{
        z->k = y->k;
        z->v = y->v;
    }
}
#define MAX_FP64                                                             \
"void max_fp64(tuple_fp64 *z, const tuple_fp64 *x, const tuple_fp64 *y){ \n" \
"   if (x->v > y->v)                                                     \n" \
"   {                                                                    \n" \
"       z->k = x->k;                                                     \n" \
"       z->v = x->v;                                                     \n" \
"   }else{                                                               \n" \
"       z->k = y->k;                                                     \n" \
"       z->v = y->v;                                                     \n" \
"    }                                                                   \n" \
"}"          
#define MAKE_FP64                                                \
"void make_fp64(tuple_fp64 *z,                               \n" \
"               const double *x, GrB_Index ix, GrB_Index jx, \n" \
"               const uint64_t *y, GrB_Index iy, GrB_Index jy,    \n" \
"               const void *theta)                           \n" \
"{                                                           \n" \
"    z->k = (int64_t)jx;                                     \n" \
"    z->v = (*x) + (((*y) ^ (jy) << iy));                    \n" \
"       printf(\"z->k:%ld, %f\t \",z->k,z->v);                \n" \   
"}"



int LAGraph_Louvain2(
    //output
    GrB_Matrix *S_result,   // TODO: make this a vector
    // input
    LAGraph_Graph G,
    char* msg
)
{
    LG_CLEAR_MSG ;

    char MATRIX_TYPE[LAGRAPH_MSG_LEN];
    if (DEBUG)
        GrB_set (GrB_GLOBAL, true, GxB_BURBLE);

    //assignment of monoids, bops, and semis   
    GrB_Monoid plusmon = GrB_PLUS_MONOID_FP64;
    GrB_Monoid maxmon = GrB_MAX_MONOID_FP64;

    GrB_BinaryOp plusf64 = GrB_PLUS_FP64;
    GrB_BinaryOp timesf64 = GrB_TIMES_FP64;
    GrB_BinaryOp minusf64 = GrB_MINUS_FP64;


    GrB_Semiring stdmxm = GrB_PLUS_TIMES_SEMIRING_FP64;
    GrB_Semiring anypB = GxB_ANY_PAIR_FP64 ;

    double *Sx = NULL ; //try S as double not bool
    GrB_Index *Sp = NULL , *Sj = NULL, Sp_size = 0, Sj_size = 0, Sx_size = 0 ;
    bool S_jumbled = false, S_iso = false;
    GrB_Vector t_q = NULL, q1=NULL, t=NULL,v=NULL;
    GrB_Vector k = NULL ;
    GrB_Vector x = NULL ;
    GrB_Vector z = NULL ;
    GrB_Index n,b;
    GrB_Matrix S = NULL;
    GrB_Vector sr = NULL;
    // add these to GB_FREE_ALL:
    GrB_Semiring Semiring = NULL ;
    GrB_Vector srxq = NULL;
    GrB_Index vals_srxq;
    GrB_Matrix dS = NULL ;
    GrB_Vector dSk = NULL, vtS = NULL ;
    GrB_Vector temp = NULL ;
    GrB_Vector y_rand = NULL ;
    GrB_Vector max_q1 = NULL ;
    GxB_IndexBinaryOp Iop = NULL ;
    GrB_BinaryOp Bop = NULL, MonOp = NULL ;
    GrB_Scalar Theta = NULL ;
    GrB_Type Tuple = NULL ;
    GrB_Monoid Mon = NULL ;
    double *dSx = NULL ; 
    GrB_Index *dSp = NULL, *dSj = NULL, dSp_size, dSj_size, dSx_size ;
    bool dS_jumbled = false, dS_iso = false ;
    GrB_Index q1_size;
    tuple_fp64 o;
    double o1;
    double k_i = 0;
    // FIXME: add check to see if S_result is NULL

    (*S_result) = NULL ;

    GrB_Matrix A = G->A;
    // GxB_print(A,5);

//index bin op definitions
//------------------------------------------------------------------------------------------------------------------

    GRB_TRY(GrB_Scalar_new(&Theta, GrB_BOOL));
    GRB_TRY(GrB_Scalar_setElement_BOOL(Theta, 0));
    GRB_TRY(GxB_Type_new(&Tuple, sizeof(tuple_fp64), "tuple_fp64", FP64_K));
    GRB_TRY(GxB_IndexBinaryOp_new(&Iop,(GxB_index_binary_function)make_fp64, Tuple, GrB_FP64, GrB_UINT64, GrB_BOOL,"make_fp64", MAKE_FP64));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&Bop, Iop, Theta));
    tuple_fp64 id;
    memset(&id, 0, sizeof(tuple_fp64));
    id.k = INT64_MAX;
    id.v = (double)(-INFINITY);
    GRB_TRY(GxB_BinaryOp_new(&MonOp,(GxB_binary_function)max_fp64, Tuple, Tuple, Tuple, "max_fp64", MAX_FP64));
    GRB_TRY(GrB_Monoid_new_UDT(&Mon, MonOp, &id));
    GRB_TRY(GrB_Semiring_new(&Semiring, Mon, Bop));
//------------------------------------------------------------------------------------------------------------------

    GRB_TRY(LAGraph_Random_Init(msg));

    GRB_TRY(GrB_Matrix_nrows(&n,A));
    GRB_TRY(GrB_Matrix_ncols(&b,A));

    //k = [+_j A(:,j)]
    GRB_TRY(GrB_Vector_new(&k, GrB_FP64, n));
    GRB_TRY (GrB_Matrix_reduce_Monoid(k, NULL, NULL,plusmon, A, NULL));
    // GxB_print(k,5);

    //m = .5[+_i k(i)]
    double m;
    GRB_TRY(GrB_Vector_reduce_FP64(&m,NULL,plusmon,k,NULL));
    m*=.5;
    // printf("m= %f\n", m);

    // S <- I
    GRB_TRY(GrB_Vector_new(&x,GrB_FP64,n));
    GRB_TRY(GrB_assign (x, NULL, NULL, 1, GrB_ALL, n, NULL)) ;
    // GxB_print(i,5);
    GRB_TRY(GrB_Matrix_diag(&S,x,0));
    // GxB_print(S,5);
    GRB_TRY(GrB_Vector_new(&v,GrB_FP64,n));
    //var used in for loop
    GRB_TRY(GrB_Vector_new(&t_q, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&q1, GrB_FP64, n)); 
    GRB_TRY(GrB_Vector_new(&z,GrB_FP64,n));
    GRB_TRY(GrB_Vector_new(&dSk,GrB_FP64,n));
    GRB_TRY(GrB_Vector_new(&vtS,GrB_FP64,n));
    // temp used to  set dS to 0 matrix
    GRB_TRY(GrB_Vector_new(&temp,GrB_FP64,n));
    GRB_TRY(GrB_assign(temp, NULL, NULL, 0, GrB_ALL,n,NULL));

    GRB_TRY(GrB_Vector_new(&max_q1,Tuple,1));
    GRB_TRY(GrB_Vector_new(&srxq,GrB_BOOL,n));
    GRB_TRY(GrB_Matrix_diag(&dS,temp,0));  
    GRB_TRY(GrB_Vector_new(&y_rand, GrB_UINT64,n));
    GRB_TRY(GrB_Vector_new(&sr, GrB_FP64,n));
    GRB_TRY(GrB_Vector_new(&srxq, GrB_FP64,n));
    GRB_TRY (GrB_assign (y_rand, NULL, NULL, 1, GrB_ALL, n, NULL));

    bool changed = true;
    int max_iter = 20;
    int iter =0;
    uint64_t seed = 210;
    GRB_TRY(GrB_mxv(z,NULL,NULL,stdmxm,S,k,NULL));
    while(changed && iter < max_iter){
        changed = false;
        for(int i=0;i<n;i++){//extract tuples
            // v = A(i,:)
            GRB_TRY (GrB_Col_extract (v, NULL, NULL, A, GrB_ALL, b, i,GrB_DESC_T0));
            // GxB_print(v,5);
            // -- extract k_i
            GRB_TRY(GrB_Vector_extractElement_FP64(&k_i,k,i));
            
            //t_q =v any.pair S   O(|v|)
            GRB_TRY(GrB_vxm(t_q,NULL,NULL,anypB,v,S,NULL));
            // GxB_print(t_q,5);

            //sr = S(i,:)
            GRB_TRY(GrB_Col_extract(sr,NULL,NULL,S,GrB_ALL,1,i,GrB_DESC_T0));
            //S(i,:) = empty
            GRB_TRY (GxB_Matrix_unpack_CSR (S, &Sp, &Sj, (void ** )&Sx,
                &Sp_size, &Sj_size, &Sx_size, NULL, NULL, NULL)) ;
            Sx[i] = false;
            GRB_TRY (GxB_Matrix_pack_CSR (S, &Sp, &Sj, (void**)&Sx,
                Sp_size, Sj_size, Sx_size, false, false, NULL));

////////////////////////////////////////////////////////////
//-------------q1<t_q> = a(kTS)+vTS----------- -----------//

            double alpha = -k_i/m;
            //compute dS
            //if version 10 
             GRB_TRY (GxB_Matrix_unpack_CSR (dS, &dSp, &dSj, (void ** )&dSx,
                &dSp_size, &dSj_size, &dSx_size, NULL, NULL, NULL)) ;
            dSx[i] = -1;
            if(i>0) {dSx[i-1] =0;}
            GRB_TRY (GxB_Matrix_pack_CSR (dS, &dSp, &dSj, (void**)&dSx,
                dSp_size, dSj_size, dSx_size, false, false, NULL));
            // GxB_print(dS,5);

            //compute z
            GRB_TRY(GrB_mxv(dSk,NULL,NULL,stdmxm,dS,k,GrB_DESC_T0));
            // GxB_print(dSk,5);
            GRB_TRY(GrB_Vector_eWiseAdd_BinaryOp(z,NULL,NULL,plusf64,z,dSk,NULL));
            // GxB_print(z,5);
            GRB_TRY(GrB_Vector_apply_BinaryOp2nd_FP64(z,NULL,NULL,timesf64,dSk,alpha,NULL));
            // GxB_print(z,5);
            
            // vtS
            GRB_TRY(GrB_vxm(vtS,NULL,NULL,stdmxm,v,S,GrB_DESC_T0));
            // GxB_print(vtS,5);
            GRB_TRY(GrB_Vector_eWiseAdd_BinaryOp(q1,t_q,NULL,plusf64,z,vtS,GrB_DESC_RT0));
            // GxB_print(q1,5);
///////////////////////////////////////////////////////////
            
///////////////////////////////////////////////////////////
//-------------Index Binary OP Rand_argminmax -----------//
            GRB_TRY(GrB_Vector_nvals(&q1_size,q1));
            GxB_print(q1,5);
            printf("Size of q1: %ld\n",q1_size);

            seed++;
            // GRB_TRY(LAGraph_Random_Seed(y_rand,seed,msg));
            // GRB_TRY(GrB_Vector_setElement_UINT64(y_rand,seed,0));
            GRB_TRY (GrB_assign (y_rand, NULL, NULL, seed, GrB_ALL, n, NULL));
            // GxB_print(y_rand,5);
            // GxB_print(q1,5); 
            GRB_TRY(GrB_mxv(max_q1,NULL,NULL,Semiring,(GrB_Matrix)q1,y_rand,GrB_DESC_T0));

            GRB_TRY(GrB_Vector_extractElement_UDT((void*)&o,max_q1,0));
            printf("choice:%ld\n",(long)o.k);
            GRB_TRY (GxB_Matrix_unpack_CSR (S, &Sp, &Sj, (void ** )&Sx,
                &Sp_size, &Sj_size, &Sx_size, NULL, &S_jumbled, NULL)) ;
            Sj[i] = o.k;
            Sx[i] = true;
            GRB_TRY (GxB_Matrix_pack_CSR (S, &Sp, &Sj, (void**)&Sx,
                Sp_size, Sj_size, Sx_size, false, S_jumbled, NULL));
            // }else{
            //     GRB_TRY(GrB_Vector_extractElement_FP64(&o1,q1,0));
            //     printf("Monochoice:%f\n",o1);
            //     GRB_TRY (GxB_Matrix_unpack_CSR (S, &Sp, &Sj, (void ** )&Sx,
            //         &Sp_size, &Sj_size, &Sx_size, NULL, &S_jumbled, NULL)) ;
            //     Sj[i] = o1;
            //     Sx[i] = true;
            //     GRB_TRY (GxB_Matrix_pack_CSR (S, &Sp, &Sj, (void**)&Sx,
            //         Sp_size, Sj_size, Sx_size, false, S_jumbled, NULL));
            // }
//////////////////////////////////////////////////////////

            // GxB_print(sr,5);
            GRB_TRY(GrB_Vector_eWiseMult_BinaryOp(srxq,NULL,NULL,timesf64,sr,q1,NULL));
            GRB_TRY(GrB_Vector_nvals(&vals_srxq,srxq));
            // GxB_print(srxq,5);
            // printf("%ld",vals_srxq);
            // printf("%d",vals_srxq==0);
            if(vals_srxq==0){
                changed  = true;
            }
        }
        iter++;
    }
    // GxB_print(S,5);
    double Q;
    double gamma = 1;
    GRB_TRY(LAGr_Modularity2(&Q,gamma,A,S,msg));
    printf("Iterations: %d\n", iter);
    printf("Q:%.15g\n",Q);
    (*S_result) = S ;
    S = NULL;
    LG_FREE_ALL;
    return 0;
}
