//------------------------------------------------------------------------------
// LAGraph_Louvain.c: Runs the Louvain Algorithm on a given graph (Under construction)
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
#define LG_FREE_ALL                   \
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


double rd() {
    uint64_t r53 = ((uint64_t)(rand()) << 21) ^ (rand() >> 2);
    return (double)r53 / 9007199254740991.0; // 2^53 - 1
};
int LAGraph_Louvain(
    //output
    GrB_Matrix S,
    // input
    LAGraph_Graph G,
    char* msg
)
{
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];
    GrB_set (GrB_GLOBAL, false, GxB_BURBLE);
    //assignment of monoids, bops, and semis   
    GrB_Monoid plusmon = GrB_PLUS_MONOID_FP64;
    GrB_Monoid maxmon = GrB_MAX_MONOID_FP64;

    GrB_BinaryOp plusf64 = GrB_PLUS_FP64;
    GrB_BinaryOp timesf64 = GrB_TIMES_FP64;
    GrB_BinaryOp minusf64 = GrB_MINUS_FP64;


    GrB_Semiring stdmxm = GrB_PLUS_TIMES_SEMIRING_FP64;
    GrB_Semiring anypB = GxB_ANY_PAIR_FP64 ;

    bool *Sx; //try S as double not bool
    GrB_Index *Sp, *Sj, Sp_size, Sj_size, Sx_size ;
    bool S_jumbled, S_iso;
    GrB_Vector t_q, sr, q, q1, t, p,v;
    GrB_Vector srxt;
    GrB_Vector k;
    GrB_Vector x;
    GrB_Index nvals_srxt,nvals_t;
    GrB_Index *coor = NULL;
    bool * vals;
    GrB_Index *p_cs=NULL;
    double * p_vals;
    GxB_Container S_container = NULL;

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
    GRB_TRY(GrB_Matrix_new(&S, GrB_BOOL, n,n));
    GRB_TRY(GrB_Vector_new(&x,GrB_BOOL,n));
    GRB_TRY(GrB_assign (x, NULL, NULL, 1, GrB_ALL, n, NULL)) ;
    // GxB_print(i,5);
    GRB_TRY(GrB_Matrix_diag(&S,x,0));
    GRB_TRY(GrB_set(S,false,GxB_ISO));
    GrB_set (S, GxB_SPARSE, GxB_SPARSITY_CONTROL);
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
    GRB_TRY(GxB_Container_new(&S_container));

    // int64_t vc = vertices_changed;
    bool changed = true;
    int max_iter = 20;
    int iter =0;
    while(changed && iter < max_iter){
        changed = false;
        double k_i;
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
            GRB_TRY(GxB_unload_Matrix_into_Container(S,S_container,NULL));
            GRB_TRY(GrB_Vector_setElement_BOOL(S_container->x,false,i));
            GRB_TRY(GxB_load_Matrix_from_Container(S,S_container,NULL));


            double alpha = -k_i/m;

            //q <- k O(n)
            GrB_free(&q);
            GRB_TRY(GrB_Vector_dup(&q,k));//change to make faster
            

            //q<k> *= alpha O(n)
            GRB_TRY(GrB_Vector_apply_BinaryOp2nd_FP64(q,k,NULL,timesf64,q,alpha,NULL));

            //q += v
            GRB_TRY(GrB_eWiseAdd(q,NULL,NULL,plusf64,q,v,NULL));
            // GxB_print(q,5);

            //q_1<t_q> = q +.x S O(n)
            GRB_TRY(GrB_Vector_clear(q1));
            // GxB_print(S,5);
            GRB_TRY(GrB_vxm(q1,t_q,NULL,stdmxm,q,S,GrB_DESC_S));
            // GxB_print(q1,5);
            //t = (q1 == [max_i q_1(i)])
            double max_q1=0;
            GRB_TRY(GrB_Vector_reduce_FP64(&max_q1,NULL,GrB_MAX_MONOID_FP64,q1,NULL));
            // printf("%ld",max_q1);
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
                free(p_cs);
                free(p_vals);
            }
            // GxB_print(t,5);

            //S(i:)=t
            GRB_TRY (LAGraph_Malloc ((void **) &coor, nvals_t, sizeof (GrB_Index), msg));
            GRB_TRY (LAGraph_Malloc ((void **) &vals, nvals_t, sizeof (bool), msg)) ;
            GRB_TRY(GrB_Vector_extractTuples_BOOL(coor,vals,&nvals_t,t));
            GRB_TRY(GxB_unload_Matrix_into_Container(S,S_container,NULL));
            GRB_TRY(GrB_Vector_setElement(S_container->i,coor[0],i));
            GRB_TRY(GrB_Vector_setElement_BOOL(S_container->x,true,i));
            GRB_TRY(GxB_load_Matrix_from_Container(S,S_container,NULL));
            // GRB_TRY (GxB_Matrix_unpack_CSR (S, &Sp, &Sj, (void ** )&Sx,
            //     &Sp_size, &Sj_size, &Sx_size, NULL, &S_jumbled, NULL)) ;
            // Sj[i] = coor[0];
            // Sx[i] = true;
            // GRB_TRY (GxB_Matrix_pack_CSR (S, &Sp, &Sj, (void**)&Sx,
            //     Sp_size, Sj_size, Sx_size, NULL, S_jumbled, NULL));
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
    // LG_FREE_ALL;
    LAGraph_Random_Finalize(msg);
    return 0;
}
