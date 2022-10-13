#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_Coarsen_Matching
(
    // outputs:
    GrB_Matrix *coarsened, // coarsened adjacency
    // inputs:
    LAGraph_Graph G,
    int matching_type,
    char *msg
)
{
    LG_CLEAR_MSG ;

    return (GrB_SUCCESS) ;
}