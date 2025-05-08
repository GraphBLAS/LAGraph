#include "LG_internal.h"
#include <LAGraphX.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#undef LG_FREE_ALL
#define LG_FREE_ALL\
{   \
    GrB_free(&score) ; \
    GrB_free(&scoreA) ; \
    GrB_free(&neighbor_max) ; \
    GrB_free(&new_members) ; \
    GrB_free(&new_membersA) ; \
    GrB_free(&new_neighbors) ; \
    GrB_free(&candidates) ; \
    GrB_free(&empty) ; \
    GrB_free(&Seed) ; \
    GrB_free(&degree) ; \
}
#define DEBUG 0
#define dbg(x) if (DEBUG) GxB_print(x,5)
typedef GrB_Matrix mat;
typedef GrB_Vector vec ;
typedef GrB_Scalar sca; 
typedef GrB_Index idx;
int LAGraph_IsolateSets(
    //output
    GrB_Vector *isolate_set,
    //input
    LAGraph_Graph G,
    uint64_t seed,
    char* msg
){
    LG_CLEAR_MSG ;

    char MATRIX_TYPE[LAGRAPH_MSG_LEN];
    GrB_set (GrB_GLOBAL, DEBUG, GxB_BURBLE);
    GrB_Vector iset = NULL ;            // independent set (output vector)
    GrB_Vector score = NULL;
    GrB_Vector scoreA = NULL ;           // random score for each node
    // GrB_Vector scoreAA = NULL;
    GrB_Vector neighbor_max = NULL ;    // value of max neighbor score
    GrB_Vector new_members = NULL ;     // set of new members to add to iset
    GrB_Vector new_membersA = NULL;
    GrB_Vector new_neighbors = NULL ;   // new neighbors to new iset members
    GrB_Vector candidates = NULL ;      // candidate nodes
    GrB_Vector empty = NULL ;           // an empty vector
    GrB_Vector Seed = NULL ;            // random number seed vector
    GrB_Vector degree = NULL ;          // (float) G->out_degree
    GrB_Matrix A ;                      // G->A, the adjacency matrix
    GrB_Index n ;                       // # of nodes

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT(isolate_set != NULL, GrB_NULL_POINTER);
    A = G->A;
    dbg(A);

    GRB_TRY (GrB_Matrix_nrows(&n,A));
    GRB_TRY (GrB_Vector_new(&iset,GrB_BOOL,n));
    GRB_TRY (GrB_Vector_new (&neighbor_max, GrB_FP32, n)) ;
    GRB_TRY (GrB_Vector_new (&degree, GrB_FP32, n)) ;
    GRB_TRY (GrB_Vector_new (&new_members, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new (&new_neighbors, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new(&new_membersA,GrB_BOOL,n));
    GRB_TRY (GrB_Vector_new (&candidates, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new (&empty, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new (&Seed, GrB_UINT64, n)) ;
    GRB_TRY (GrB_Vector_new(&score,GrB_FP32, n));
    GRB_TRY (GrB_Vector_new (&scoreA, GrB_FP32, n)) ;

    //rand
    seed = (uint64_t)time(NULL);
    printf("%ld",seed);
    GRB_TRY (GrB_assign (Seed, NULL, NULL, 1, GrB_ALL, n, NULL));
    GRB_TRY (LAGraph_Random_Seed (Seed, seed, msg)) ;
    dbg(Seed);


    GrB_Index ncandidates ;
    GRB_TRY (GrB_assign (candidates, NULL, NULL, (bool) true, GrB_ALL,n, NULL)) ;
    GRB_TRY (GrB_Vector_nvals (&ncandidates, candidates)) ;

    GRB_TRY(GrB_assign(degree,NULL,NULL,G->out_degree,GrB_ALL,n,NULL));
    dbg(degree);
    GRB_TRY (GrB_assign (score, NULL, NULL, Seed, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_eWiseMult (score, NULL, NULL, GrB_DIV_FP32, score, degree,NULL)) ;
    dbg(score);

    dbg(candidates);
    GRB_TRY (GrB_vxm (scoreA, candidates, NULL,
        GrB_MAX_FIRST_SEMIRING_FP32, score, A, GrB_DESC_RS)) ;
    dbg(scoreA);
    GRB_TRY (GrB_vxm (neighbor_max, candidates, NULL,GrB_MAX_FIRST_SEMIRING_FP32, scoreA, A, GrB_DESC_RS)) ;
    dbg(neighbor_max);
    dbg(score);

    GRB_TRY (GrB_eWiseAdd (new_members, NULL, NULL, GrB_GE_FP32,
        score, neighbor_max, NULL)) ;
    dbg(new_members);
    GRB_TRY (GrB_select (new_members, NULL, NULL, GrB_VALUEEQ_BOOL,
        new_members, (bool) true, NULL)) ;
    dbg(new_members);
    GRB_TRY (GrB_assign (iset, new_members, NULL,true,GrB_ALL,n,NULL)) ;
    // GRB_TRY (GrB_assign (candidates, new_members, NULL, empty,
    //     GrB_ALL, n, GrB_DESC_S)) ;
    // GxB_print(candidates,5);
    // GrB_Index n_new_members ;
    // GRB_TRY (GrB_Vector_nvals (&n_new_members, new_members)) ;

    // GRB_TRY (GrB_vxm (new_membersA, candidates, NULL,
    //     LAGraph_any_one_bool, new_members, A, GrB_DESC_RS)) ; 
    // GRB_TRY (GrB_vxm (new_neighbors, candidates, NULL,
    //     LAGraph_any_one_bool, new_membersA, A, GrB_DESC_RS)) ;
    // GxB_print(new_neighbors,5);
    (*isolate_set) = iset;
    iset = NULL;
    LG_FREE_ALL;
    return 0;
}