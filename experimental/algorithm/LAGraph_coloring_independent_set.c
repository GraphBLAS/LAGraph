#include "LG_internal.h" // contains all internal grb operations
#include "LAGraphX.h"    // algorithm added to LAGraphX.h

#undef  LG_FREE_WORK
#define LG_FREE_WORK                \
    GrB_free (&local_color) ;       \
    GrB_free (&weight) ;            \
    GrB_free (&in_curr_subset) ;    \
    GrB_free (&max_weights) ;

#undef  LG_FREE_ALL
#define LG_FREE_ALL                 \
    LG_FREE_WORK ;

int LAGraph_coloring_independent_set
(
    // output
    GrB_Vector *color,
    int *num_colors,

    // input
    LAGraph_Graph G,
    char *msg
)
{
#if LAGRAPH_SUITESPARSE

    bool verbose = false;
    GrB_Vector local_color = NULL;
    GrB_Vector weight = NULL;
    GrB_Vector in_curr_subset = NULL;
    GrB_Vector max_weights = NULL;

    GrB_Index n;
    GRB_TRY (GrB_Matrix_nrows(&n, G->A)) ;

    GrB_Type Int = (n < UINT32_MAX) ? GrB_UINT32 : GrB_UINT64 ;

    /* initialize local copy of color to SPARSE vector */
    GRB_TRY(GrB_Vector_new(&local_color, Int, n));

    // lg_set_format_hint -> bitmap
    GRB_TRY (LG_SET_FORMAT_HINT (local_color, LG_BITMAP)) ;

    /* weights initialized randomly
    *  seed of 20 was chosen arbitrarily */   
    GRB_TRY(GrB_Vector_new(&weight, GrB_UINT64, n));
    GRB_TRY(GrB_assign (weight, NULL, NULL, 0, GrB_ALL, n, NULL));

    // LG_TRY(LAGraph_Random_Seed(weight, 2, msg));
    LG_TRY (LAGraph_Random_Seed(weight, 20, msg)) ;

    GRB_TRY(GrB_Vector_new(&in_curr_subset, GrB_BOOL, n));

    GRB_TRY(GrB_Vector_new(&max_weights, GrB_UINT64, n));
    double tlast = LAGraph_WallClockTime ( ) ;
    double tnow = 0 ;
    LG_SET_BURBLE(true) ;

    /* algorithm start */
    int64_t curr_color;
    for (curr_color = 1; curr_color < n+1; curr_color++) {
        /* mxv - find maximum of all neighboring weights */

        // FUTURE WORK: try using a set of sparse candidate nodes, not yet colored
        GRB_TRY(GrB_mxv(max_weights, local_color, GrB_NULL,
            GrB_MAX_SECOND_SEMIRING_UINT64, G->A, weight, GrB_DESC_RSC));

        /* eWiseAdd - 1 if current weight > max neighboring weight */
        GRB_TRY(GrB_eWiseMult(in_curr_subset, GrB_NULL, GrB_NULL, GrB_GT_UINT64, weight, max_weights, GrB_NULL));

        /* select - select all entries in in_curr_subset that are true, and delete falses */
        GRB_TRY(GrB_select(in_curr_subset, GrB_NULL, GrB_NULL, GrB_VALUEEQ_BOOL, in_curr_subset, true, GrB_NULL));
        tnow = LAGraph_WallClockTime ( ) ;
        double tthis = tnow - tlast ;
        tlast = tnow ;

        GrB_Index nvals_local_color;
        GRB_TRY(GrB_Vector_nvals(&nvals_local_color, local_color));

        /* check if in_curr_subset is empty then break */
        GrB_Index nvals_in_curr_subset;
        GRB_TRY(GrB_Vector_nvals(&nvals_in_curr_subset, in_curr_subset));
        if (nvals_in_curr_subset == 0) { 
            GrB_Index nvals_local_color;
            GRB_TRY(GrB_Vector_nvals(&nvals_local_color, local_color));
            if (nvals_local_color < n) {
                LG_ASSERT_MSG (false, LAGRAPH_CONVERGENCE_FAILURE,
                    "LAGraph_coloring_independent_set: in_curr_subset is empty, but nvals (local_color) < n") ;
            }
            break;
        }

        /* assign - write current color to C vector according to in_curr_subset mask */
        GRB_TRY(GrB_assign(local_color, in_curr_subset, GrB_NULL, curr_color, GrB_ALL, n, GrB_DESC_S));
        
        /* assign - write 0 to weight according to in_curr_subset mask */
        GRB_TRY(GrB_assign(weight, in_curr_subset, GrB_NULL, 0, GrB_ALL, n, GrB_DESC_S));
    }
    
    LG_SET_BURBLE(false) ;

    (*num_colors) = curr_color - 1;
    (*color) = local_color;
    local_color = NULL ;
    LG_FREE_ALL ;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
