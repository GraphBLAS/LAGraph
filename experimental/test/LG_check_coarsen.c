#include "LG_internal.h"
#include "LG_test.h"
#include "LG_test.h"
#include "LG_Xtest.h"

// #undef LG_FREE_ALL
// #undef LG_FREE_WORK

int LG_check_coarsen
(
    // outputs:
    GrB_Matrix *coarsened,    // coarsened adjacency
    // inputs:
    GrB_Matrix A,          // input adjacency (for the purposes of testing, is FP64)
    GrB_Vector parent,     // parent mapping (must be compressed, i.e. p[p[i]] = p[i] for all i)
    int preserve_mapping,  // whether to preserve the original namespace of nodes
    char *msg
)
{
#if 0
    Steps:
    - Create output matrix
        - If preserve_mapping, same as input dim
        - Else, count p[i] != i,  
    extractTuples on A to get all edges
    for each edge:
        - if they have the same parent, skip
        - otherwise, 
#endif
    GrB_Matrix result = NULL ;

    GrB_Index n ;
    GrB_Index n_new = n ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;

    if (preserve_mapping) {
        GRB_TRY (GrB_Matrix_new (&result, GrB_FP64, n, n)) ;
    } else {
        for (GrB_Index i = 0 ; i < n ; i++) {
            uint64_t par ;
            GRB_TRY (GrB_Vector_extractElement (&par, parent, i)) ;
            if (par != i) {
                n_new-- ;
            }
        }
        GRB_TRY (GrB_Matrix_new (&result, GrB_FP64, n_new, n_new)) ;
    }

}