#include "LG_internal.h" // contains all internal grb operations
#include "LAGraphX.h"    // algorithm added to LAGraphX.h

#undef  LG_FREE_WORK
#define LG_FREE_WORK                \
{                                   \
    GrB_free (&local_color) ;       \
    GrB_free (&curr_MIS) ;          \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
}

int LAGraph_coloring_MIS
(
    // output
    GrB_Vector *color,
    int *num_colors,

    // input
    LAGraph_Graph G,
    char *msg
)
{
    bool verbose = false;
    GrB_Vector local_color = NULL;
    GrB_Vector curr_MIS = NULL;
    GrB_Index n;
    GRB_TRY (GrB_Matrix_nrows(&n, G->A)) ;

    GrB_Type Int = (n < UINT32_MAX) ? GrB_UINT32 : GrB_UINT64 ;

    /* initialize local copy of color to SPARSE vector */
    GRB_TRY(GrB_Vector_new(&local_color, Int, n));

    LAGRAPH_TRY(LAGraph_DeleteSelfEdges(G, msg)) ;
    LAGRAPH_TRY(LAGraph_Cached_OutDegree(G, msg)) ;    

    /* algorithm start */
    GrB_Index colored_nodes;
    int64_t curr_color;
    for (curr_color = 1; curr_color < n+1; curr_color++) {

        /* compute MIS for current color */
        LAGRAPH_TRY(LAGraph_MaximalIndependentSet(&curr_MIS, G, 20, local_color, msg));

        /* assign - write current color to C vector according to in_curr_subset mask */
        GRB_TRY(GrB_assign(local_color, curr_MIS, GrB_NULL, curr_color, GrB_ALL, n, GrB_DESC_S));

        GrB_Vector_nvals(&colored_nodes, local_color);
        if (colored_nodes == n) {
            break;
        }

    }
    
    (*num_colors) = curr_color;
    (*color) = local_color;
    local_color = NULL ;
    LG_FREE_ALL ;
    return (GrB_SUCCESS) ;
}
