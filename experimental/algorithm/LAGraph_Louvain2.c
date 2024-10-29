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

#define LG_FREE_LOUV2                   \
    {                                   \
        GrB_free (&A) ;                 \
        GrB_free (&k) ;                 \
        GrB_free (&x) ;                 \
        GrB_free (&v) ;                 \
        GrB_free (&sr) ;                \
        GrB_free (&q) ;                 \
        GrB_free (&q1) ;                \
        GrB_free (&t) ;                 \
        GrB_free (&p) ;                 \
        GrB_free (&srxt) ;              \
        GrB_free (&t_q) ;               \
    }
#define LG_FREE_ALL                     \
    {                                   \
        LG_FREE_LOUV2;                  \
    }
#define DEBUG 0
int LAGraph_Louvain2(
    //output
    GrB_Matrix S,
    // input
    LAGraph_Graph G,
    char* msg
)
{
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

    double *Sx; //try S as double not bool
    GrB_Index *Sp, *Sj, Sp_size, Sj_size, Sx_size ;
    bool S_jumbled, S_iso;
    GrB_Vector t_q = NULL, sr = NULL, q = NULL, q1=NULL, t=NULL, p=NULL,v=NULL;
    GrB_Vector srxt = NULL;
    GrB_Vector k;
    GrB_Vector x;
    GrB_Index nvals_srxt,nvals_t;
    GrB_Vector z;
    GrB_Index *coor = NULL;
    bool * vals = NULL;
    GrB_Index *p_cs=NULL;
    double * p_vals;


    GrB_Matrix A = G->A;
    // GxB_print(A,5);
    GrB_Index n,b;
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

    //var used in for loop
    GrB_Index vertices_changed;
    GRB_TRY(GrB_Vector_nvals(&vertices_changed,k)); 
    GRB_TRY(GrB_Vector_new(&v, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&srxt, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&t_q, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&sr, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&q, GrB_FP64, n));  
    GRB_TRY(GrB_Vector_new(&q1, GrB_FP64, n)); 
    GRB_TRY(GrB_Vector_new(&t, GrB_BOOL, n)); 
    GRB_TRY(GrB_Vector_new(&p,GrB_FP64,n));
    GRB_TRY(GrB_Vector_new(&z,GrB_FP64,n));

    GrB_Matrix dS;
    GrB_Vector dSk,vtS;
    GRB_TRY(GrB_Vector_new(&dSk,GrB_FP64,n));
    GRB_TRY(GrB_Vector_new(&vtS,GrB_FP64,n));
    // temp used to  set dS to 0 matrix
    GrB_Vector temp;
    GRB_TRY(GrB_Vector_new(&temp,GrB_FP64,n));
    GRB_TRY(GrB_assign(temp, NULL, NULL, 0, GrB_ALL,n, NULL));
    double *dSx; 
    GrB_Index *dSp, *dSj, dSp_size, dSj_size, dSx_size ;
    bool dS_jumbled, dS_iso;
    GRB_TRY(GrB_Matrix_diag(&dS,temp,0));



    bool changed = true;
    int max_iter = 20;
    int iter =0;

    while(changed && iter < max_iter){
        changed = false;
        double k_i;
        GRB_TRY(GrB_mxv(z,NULL,NULL,stdmxm,S,k,NULL));
        for(int i=0;i<n;i++){//extract tuples
            // v = A(i,:)
            GRB_TRY (GrB_Col_extract (v, NULL, NULL, A, GrB_ALL, b, i,GrB_DESC_T0));
            // GxB_print(v,5);

            // -- extract k_i
            GRB_TRY(GrB_Vector_extractElement_FP64(&k_i,k,i));
            
            //t_q =v any.pair S   O(|v|)
            GRB_TRY(GrB_vxm(t_q,NULL,NULL,anypB,v,S,NULL));
            // GxB_print(t_q,5);

            // sr = S(i,:)
            GRB_TRY(GrB_Col_extract(sr,NULL,NULL,S,GrB_ALL,1,i,GrB_DESC_T0));
            // GxB_print(sr,5);

            //S(i,:) = empty
            GRB_TRY (GxB_Matrix_unpack_CSR (S, &Sp, &Sj, (void ** )&Sx,
                &Sp_size, &Sj_size, &Sx_size, NULL, &S_jumbled, NULL)) ;
            Sx[i] = false;
            GRB_TRY (GxB_Matrix_pack_CSR (S, &Sp, &Sj, (void**)&Sx,
                Sp_size, Sj_size, Sx_size, NULL, S_jumbled, NULL));
////////////////////////////////////////////////////////////

            double alpha = -k_i/m;

            //z += dS^t*k
            //compute dS
             GRB_TRY (GxB_Matrix_unpack_CSR (dS, &dSp, &dSj, (void ** )&dSx,
                &dSp_size, &dSj_size, &dSx_size, NULL, &dS_jumbled, NULL)) ;
            dSx[i] = -1;
            if(i>0) {dSx[i-1] =0;}
            GRB_TRY (GxB_Matrix_pack_CSR (dS, &dSp, &dSj, (void**)&dSx,
                dSp_size, dSj_size, dSx_size, NULL, dS_jumbled, NULL));
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

            //Compute q1
            //q1<t_q> = alpha (ktS) + (vtS)
            // GxB_print(t_q,5);
            GRB_TRY(GrB_Vector_eWiseAdd_BinaryOp(q1,t_q,NULL,plusf64,z,vtS,GrB_DESC_RT0));
            // GxB_print(q1,5);
///////////////////////////////////////////////////////////
            

            //t = (q1 == [max_i q_1(i)])
            double max_q1;
            GRB_TRY(GrB_Vector_reduce_FP64(&max_q1,NULL,maxmon,q1,NULL));
            GRB_TRY(GrB_Vector_select_FP64(t,NULL,NULL,GrB_VALUEEQ_FP64,q1,max_q1,NULL));
            // GxB_print(t,5);

            
            GRB_TRY(GrB_Vector_nvals(&nvals_t,t));
            while(nvals_t>1){
                // printf("%ld\n",nvals_t);
                GRB_TRY(LAGraph_Malloc ((void **) &p_cs, nvals_t, sizeof (GrB_Index), msg)); //free p_cs and P-vals
                GRB_TRY(LAGraph_Malloc ((void **) &p_vals, nvals_t, sizeof (double), msg)) ;
                GRB_TRY(GrB_Vector_extractTuples_FP64(p_cs,p_vals,&nvals_t,t));

                // p = random() x t
                for(int j = 0;j<nvals_t;j++){
                    double y = rd();
                    GRB_TRY(GrB_Vector_setElement_FP64(p,y*p_vals[j],p_cs[j]));
                }
                // GxB_print(p,5);
                //t = (p== [max_i p_1(i)])
                double max_p;
                GRB_TRY(GrB_Vector_reduce_FP64(&max_p,NULL,maxmon,p,NULL));
                // printf("max_p:%f\n",max_p);
                GRB_TRY(GrB_Vector_select_FP64(t,NULL,NULL,GrB_VALUEEQ_FP64,p,max_p,NULL));
                // GxB_print(t,5);
                GRB_TRY(GrB_Vector_nvals(&nvals_t,t));
                // LAGraph_Free((void**)p_cs,msg);
                free(p_cs);
                free(p_vals);
            }
            // GxB_print(t,5);

            //S(i:)=t
            GRB_TRY (LAGraph_Malloc ((void **) &coor, nvals_t, sizeof (GrB_Index), msg));//declare statically
            GRB_TRY (LAGraph_Malloc ((void **) &vals, nvals_t, sizeof (bool), msg)) ;
            GRB_TRY(GrB_Vector_extractTuples_BOOL(coor,vals,&nvals_t,t));
            GRB_TRY (GxB_Matrix_unpack_CSR (S, &Sp, &Sj, (void ** )&Sx,
                &Sp_size, &Sj_size, &Sx_size, NULL, &S_jumbled, NULL)) ;
            Sj[i] = coor[0];
            Sx[i] = true;
            GRB_TRY (GxB_Matrix_pack_CSR (S, &Sp, &Sj, (void**)&Sx,
                Sp_size, Sj_size, Sx_size, NULL, S_jumbled, NULL));
            free(coor);
            free(vals);
            // GxB_print(S,5);
            GRB_TRY(GrB_Vector_eWiseMult_BinaryOp(srxt,NULL,NULL,timesf64,sr,t,NULL));
            GRB_TRY(GrB_Vector_nvals(&nvals_srxt,srxt));
            // GxB_print(srxt,5);
           if(nvals_srxt==0){
                // GRB_TRY(GrB_Vector_nvals(&vertices_changed,k));
                changed  = true;
            }
            // vertices_changed -=1;
            // vc = vertices_changed;
        }
        iter++;
    }
    // GxB_print(S,5);
    double Q;
    double gamma = 1;
    GRB_TRY(LAGr_Modularity2(&Q,gamma,A,S,msg));
    printf("Iterations: %d\n", iter);
    printf("Q:%.15g\n",Q);
    // (*S_result) = S ;
    S = NULL;
    LG_FREE_ALL;
    return 0;
}
