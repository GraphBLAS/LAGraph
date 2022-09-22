#include "LG_internal.h"
#include "LAGraphX.h"

#define F_UNARY(f)  ((void (*)(void *, const void *)) f)
#define F_BINARY(f) ((void (*)(void *, const void *, const void *)) f)

void one_binop_func (uint64_t *out, const uint64_t *in1, const uint64_t *in2) { (*out) = 1 ; } 
// void sub_two_unop_func (uint64_t *out, const uint64_t *in) { (*out) = ((*in) - 2); }

/*
Maximal matching algorithm interpreted as a max independent set on a line graph
Heavily influenced by MaximalIndependentSet.c
*/
int LAGraph_MaximalMatching
(
    // outputs:
    GrB_Vector *matching,
    // inputs:
    GrB_Matrix E,       // incidence
    int matching_type,  // 0 (random), 1 (heavy weight matching), 2 (light weight matching)
    uint64_t seed,      // random number seed
    char *msg
)
{
    LG_CLEAR_MSG ;

    GrB_Matrix E_t = NULL ;                 // E transpose. Maybe it's better to use 'A' descriptor instead of storing this explicitly?
    GrB_Vector score = NULL ;               // score for each edge. Computed according to matching_type
    GrB_Vector candidates = NULL ;          // set of candidate edges
    GrB_Vector Seed = NULL ;                // random number seed vector
    GrB_Vector node_degree = NULL ;         // intermediate result for computing edge degree; degree of the node
    GrB_Vector degree = NULL ;              // edge degree; number of incident edges
    GrB_Vector max_node_neighbor = NULL ;   // intermediate result for computing max edge neighbor; max edge touching a node
    GrB_Vector max_neighbor = NULL ;        // max neighbor of an edge (including itself)
    GrB_Vector new_members = NULL ;         // new edges to include in the matching
    GrB_Vector discard = NULL ;             // edges to discard from consideration
    // GrB_Vector ones = NULL ;             // all ones
    GrB_Vector result = NULL;               // resulting matching    

    GrB_Semiring PLUS_SECOND_SEMIRING_UINT64 = NULL ;
    GrB_Semiring PLUS_ONE_SEMIRING_UINT64 = NULL ;

    // GrB_UnaryOp SUB_TWO_UINT64 = NULL ;
    GrB_BinaryOp ONE_UINT64 = NULL ;

    GrB_Index num_edges ;
    GrB_Index num_nodes ;
    
    GRB_TRY (GrB_Matrix_nrows (&num_nodes, E)) ;
    GRB_TRY (GrB_Matrix_ncols (&num_edges, E)) ;
    GRB_TRY (GrB_Matrix_new (&E_t, GrB_UINT64, num_nodes, num_edges)) ;
    GRB_TRY (GrB_transpose (E_t, NULL, NULL, E, NULL)) ;

    GRB_TRY (GrB_Vector_new (&candidates, GrB_BOOL, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&Seed, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&score, GrB_FP32, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&node_degree, GrB_UINT64, num_nodes)) ;
    GRB_TRY (GrB_Vector_new (&degree, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&max_node_neighbor, GrB_FP32, num_nodes)) ;
    GRB_TRY (GrB_Vector_new (&max_neighbor, GrB_FP32, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&new_members, GrB_BOOL, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&discard, GrB_BOOL, num_edges)) ;
    // GRB_TRY (GrB_Vector_new (&ones, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&result, GrB_BOOL, num_edges)) ;

    // GRB_TRY (GrB_UnaryOp_new (&SUB_TWO_UINT64, F_UNARY(sub_two_unop_func), GrB_UINT64, GrB_UINT64)) ;
    GRB_TRY (GrB_BinaryOp_new (&ONE_UINT64, F_BINARY(one_binop_func), GrB_UINT64, GrB_UINT64, GrB_UINT64)) ;

    GRB_TRY (GrB_Semiring_new (&PLUS_SECOND_SEMIRING_UINT64, GrB_PLUS_MONOID_UINT64, GrB_SECOND_UINT64)) ;
    GRB_TRY (GrB_Semiring_new (&PLUS_ONE_SEMIRING_UINT64, GrB_PLUS_MONOID_UINT64, ONE_UINT64)) ;

    LG_TRY (LAGraph_Random_Seed (Seed, seed, msg)) ;

    // initially all edges are considered
    GRB_TRY (GrB_assign (candidates, NULL, NULL, true, GrB_ALL, 
        num_edges, NULL)) ;
    
    // GRB_TRY (GrB_assign (ones, NULL, NULL, 1, GrB_ALL, 
    //     num_edges, NULL)) ;
    
    GrB_Index ncandidates ;
    GrB_Index nfailures = 0 ; // counts how many iterations have failed due to invalid matchings

    GRB_TRY (GrB_Vector_nvals(&ncandidates, candidates)) ;

    LAGRAPH_TRY (LAGraph_Vector_Print (candidates, LAGraph_SHORT, stdout, msg)) ;

    while (ncandidates > 0) {
        // first just generate the scores again
        if (matching_type == 0) {
            // assign random scores to edges, weighed by degree

            // for each node, counts incident edges
            // we disallow edges that are already in the matching from contributing by using candidates vector
            // so, just use the Seed vector
            GRB_TRY (GrB_mxv (node_degree, NULL, NULL, PLUS_ONE_SEMIRING_UINT64, E, candidates, GrB_DESC_R)) ;

            // for each edge, sums incident edges for each node. Each edge has an excess of 2 degree, but it doesn't matter since
            // we care about relative degree
            GRB_TRY (GrB_mxv (degree, candidates, NULL, PLUS_SECOND_SEMIRING_UINT64, E_t, node_degree, GrB_DESC_RS)) ;

            // weighs the score based on degree
            GRB_TRY (GrB_eWiseMult (score, candidates, NULL, GrB_DIV_FP32, Seed, degree, GrB_DESC_S)) ;

        } else {
            // first get weights of edges in vector
            // assign random scores, but somehow want to weigh by edge weight
            // or maybe, a combination of edge weight + degree?
            if (matching_type == 1) {
                // heavy
            } else {
                // light
            }
        }

        // the actual edge selection is common regardless of matching type

        // intermediate result. Max score edge touching each node
        GRB_TRY (GrB_mxv (max_node_neighbor, NULL, NULL, GrB_MAX_SECOND_SEMIRING_FP32, E, score, GrB_DESC_R)) ;

        // Max edge touching each candidate edge, including itself
        GRB_TRY (GrB_mxv (max_neighbor, candidates, NULL, GrB_MAX_SECOND_SEMIRING_FP32, E_t, max_node_neighbor, GrB_DESC_RS)) ;

        GRB_TRY (GrB_eWiseAdd (new_members, NULL, NULL, GrB_GE_FP32, score, max_neighbor, GrB_DESC_R)) ;
    }

    (*matching) = result ;

    return (GrB_SUCCESS) ;
}
