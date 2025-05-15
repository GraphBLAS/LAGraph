#define LAGRAPH_COLORING_RETURN_VALUES
#define LAGRAPH_COLORING_INVALID_COLORING   (-5501)

#include "LG_internal.h"
#include "LG_test.h"
#include "LG_Xtest.h"

#undef  LG_FREE_WORK
#define LG_FREE_WORK                            \
{                                               \
    LAGraph_Free ((void **) &Ap, NULL) ;        \
    LAGraph_Free ((void **) &Ai, NULL) ;        \
    LAGraph_Free ((void **) &Ax, NULL) ;        \
}

int LG_check_coloring
(
    LAGraph_Graph G,
    GrB_Vector C,
    char *msg
)
{
#if LAGRAPH_SUITESPARSE
    // ------------------------------------------------
    // check if coloring is valid
    // ------------------------------------------------

    /* extract graph in CSC form
    *  CSC form:
    *  Ap: start and end indices of Ai that represent a column of A
    *  Ai: values that represent the row indices of A where there is a value
    *       - in our case, these are the neighbors' IDs
    *  Ax: the values of that edge, stored in the same order as Ai
    *       - in our case, all values are 1
    *  note: CSC is same as CSR for undirected graphs
    *        maybe use the one that's more efficient
    *        ( prevent converting from one to other )
    *
    *  convert Ap_size from bytes to indices
    *   - make sure to loop only up to Ap_size - 1
    *     when checking [Ap_index] to [Ap_index + 1]
    * 
    *  traverse through unpacked matrix and
    *  check current node's color against its neighbors
    *   - Ap_index: current node
    *   - Ai_index: a neighbor
    */
   
    GrB_Index *Ap = NULL;
    GrB_Index *Ai = NULL;
    void *Ax = NULL;
    GrB_Index Ap_size, Ai_size, Ax_size;
    GRB_TRY(GxB_Matrix_unpack_CSC(G->A, &Ap, &Ai, &Ax, &Ap_size, &Ai_size, &Ax_size, NULL, NULL, NULL));
    
    Ap_size = Ap_size / sizeof(GrB_Index);

    GrB_Index Ap_index;
    GrB_Index Ai_index;
    GrB_Index Ai_index_start, Ai_index_end;

    for (GrB_Index i = 0; i < Ap_size - 1; i++) {
        int color;
        if (GrB_Vector_extractElement(&color, C, i) != GrB_SUCCESS) {
            printf("error: node %" PRIu64 " has no assigned color!\n", i);
            return -1;
        }
    }

    int current_color, neighbor_color;
    for (Ap_index = 0; Ap_index < Ap_size - 1; Ap_index++) {
        
        Ai_index_start = Ap[Ap_index];
        Ai_index_end = Ap[Ap_index + 1];

        GRB_TRY(GrB_Vector_extractElement(&current_color, C, Ap_index));

        for (Ai_index = Ai_index_start; Ai_index < Ai_index_end; Ai_index++) {

            // skip self-edges
            if (Ai[Ai_index] != Ap_index) {
                GRB_TRY(GrB_Vector_extractElement(&neighbor_color, C, Ai[Ai_index]));
                LG_ASSERT_MSG(neighbor_color != current_color, LAGRAPH_COLORING_INVALID_COLORING, "found 2 connected nodes with the same color");
            }
        }
    }

    LG_FREE_WORK;
    return (GrB_SUCCESS);
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
