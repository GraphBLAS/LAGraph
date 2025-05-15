#include "../../src/benchmark/LAGraph_demo.h"
#include "LG_internal.h"
#include "LAGraphX.h"

#define DEFAULT_SIZE 2000
#define DEFAULT_DENSITY 0.7
#define DEFAULT_SEED 42

#define LG_FREE_ALL            \
{                              \
    GrB_Matrix_free (&A) ;     \
    GrB_Matrix_free (&D) ;     \
    GrB_Matrix_free (&Res) ;   \
    GrB_Scalar_free (&s) ;     \
}

int main(int argc, char **argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix D = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Matrix Res = NULL ;
    GrB_Scalar s = NULL ;

    bool burble = false ; 
    demo_init (burble) ;

    GrB_Index n = (argc > 2 ? atoi (argv [2]) : DEFAULT_SIZE) ;
    double density = (argc > 3 ? atof (argv [3]) : DEFAULT_DENSITY) ;
    uint64_t seed = (argc > 4 ? atoll (argv [4]) : DEFAULT_SEED) ;

    int ntrials = 10 ;

    for (int i = 0 ; i < ntrials ; i++) {
        LG_TRY (LAGraph_Random_Matrix (&D, GrB_FP64, n, n, density, seed, msg)) ;
        LG_TRY (LAGraph_Random_Matrix (&A, GrB_FP64, n, n, density, seed + 32, msg)) ;
        GRB_TRY (LG_SET_FORMAT_HINT (D, LG_SPARSE)) ;
        GRB_TRY (LG_SET_FORMAT_HINT (A, LG_SPARSE)) ;

        GRB_TRY (GrB_Scalar_new (&s, GrB_FP64)) ;
        GRB_TRY (GrB_Scalar_setElement_UINT64 (s, 2.0)) ;

        for (GrB_Index i = 0 ; i < n ; i++) {
            double val ;
            if (GrB_Matrix_extractElement (&val, D, i, i)) {
                GrB_Matrix_setElement (D, s, i, i) ;
            }
        }

        GRB_TRY (GrB_Matrix_wait (D, GrB_MATERIALIZE)) ;

        GrB_Index D_nvals, A_nvals ;
        GRB_TRY (GrB_Matrix_nvals (&A_nvals, A)) ;
        GRB_TRY (GrB_Matrix_nvals (&D_nvals, D)) ;
        printf ("nvals: A: %" PRIu64 ", D: %" PRIu64 "\n", A_nvals, D_nvals) ;

        GRB_TRY (GrB_Matrix_new (&Res, GrB_FP64, n, n)) ;
        GrB_Scalar_free (&s) ;
        GRB_TRY (GrB_Scalar_new (&s, GrB_UINT64)) ;
        GRB_TRY (GrB_Scalar_setElement_UINT64 (s, (uint64_t) 0)) ;
        GRB_TRY (GrB_Matrix_select_Scalar (Res, NULL, NULL, GrB_DIAG, D, s, NULL)) ;

        GrB_Matrix_free (&D) ;
        GRB_TRY (GrB_Matrix_new (&D, GrB_FP64, n, n)) ;        

        GRB_TRY (GrB_mxm (D, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_UINT64, A, Res, GrB_DESC_T0T1)) ;
        
        GrB_Matrix_free (&A) ;
        GrB_Matrix_free (&D) ;
        GrB_Matrix_free (&Res) ;
        GrB_Scalar_free (&s) ;
    }
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
