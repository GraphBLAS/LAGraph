#include "../../src/benchmark/LAGraph_demo.h"
#include "LG_internal.h"
#include "LAGraphX.h"

#define DEFAULT_SIZE 1500
#define DEFAULT_DENSITY 0.7
#define DEFAULT_SEED 42

#define BINOP 1
#define IJ 0
#define BITMAP 0

#define LG_FREE_ALL            \
{                              \
    GrB_Matrix_free (&A) ;     \
    GrB_Matrix_free (&Res) ;   \
    GrB_Scalar_free (&s) ;     \
}

int main(int argc, char **argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Matrix Res = NULL ;
    GrB_Scalar s = NULL ;

    bool burble = false ; 
    demo_init (burble) ;
    
    GrB_Index n = (argc > 2 ? atoi (argv [2]) : DEFAULT_SIZE) ;
    double density = (argc > 3 ? atof (argv [3]) : DEFAULT_DENSITY) ;
    uint64_t seed = (argc > 4 ? atoll (argv [4]) : DEFAULT_SEED) ;

    int ntrials = 1 ;

    for (int i = 0 ; i < ntrials ; i++) {

        LG_TRY (LAGraph_Random_Matrix (&A, GrB_UINT64, n, n, density, seed, msg)) ;
        GRB_TRY (LG_SET_FORMAT_HINT (A, BITMAP ? LG_BITMAP : LG_SPARSE)) ;

        GrB_Index A_nvals ;
        GRB_TRY (GrB_Matrix_nvals (&A_nvals, A)) ;
        printf ("nvals: A: %" PRIu64 "\n", A_nvals) ;

        GRB_TRY (GrB_Scalar_new (&s, GrB_UINT64)) ;
        
        GRB_TRY (GrB_Matrix_new (&Res, GrB_BOOL, n, n)) ;
    
        GRB_TRY (GrB_Scalar_setElement_UINT64 (s, (uint64_t) 1)) ;
        
        if (BINOP) {
            GRB_TRY (GrB_Matrix_apply_BinaryOp1st_Scalar (Res, NULL, NULL, GrB_PLUS_UINT64, s, A, NULL)) ;
        } else {
            GRB_TRY (GrB_Matrix_apply_IndexOp_Scalar (Res, NULL, NULL, (IJ ? GrB_DIAG : GrB_ROWLE), A, s, NULL)) ;
        }
        GrB_Matrix_free (&A) ;
        GrB_Matrix_free (&Res) ;
        GrB_Scalar_free (&s) ;
    }
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
