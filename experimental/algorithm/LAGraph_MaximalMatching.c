#include "LG_internal.h"
#include "LAGraphX.h"

#define F_UNARY(f)  ((void (*)(void *, const void *)) f)
#define F_BINARY(f) ((void (*)(void *, const void *, const void *)) f)

// void one_binop_func (uint64_t *out, const uint64_t *in1, const uint64_t *in2) { (*out) = 1 ; } 
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

    const GrB_Index MAX_FAILURES = 50 ;

    GrB_Matrix E_t = NULL ;                     // E transpose. Maybe it's better to use 'A' descriptor instead of storing this explicitly?
    GrB_Vector score = NULL ;                   // score for each edge. Computed according to matching_type
    GrB_Vector candidates = NULL ;              // set of candidate edges
    GrB_Vector Seed = NULL ;                    // random number seed vector
    GrB_Vector node_degree = NULL ;             // intermediate result for computing edge degree; degree of the node. Only computed once
    GrB_Vector degree = NULL ;                  // edge degree; number of incident edges. Only computed once
    GrB_Vector max_node_neighbor = NULL ;       // intermediate result for computing max edge neighbor; max edge touching a node
    GrB_Vector max_neighbor = NULL ;            // max neighbor of an edge (including itself)
    GrB_Vector new_members = NULL ;             // new edges to include in the matching
    GrB_Vector new_neighbors = NULL ;           // union of new members and their neighbor edges
    GrB_Vector new_members_nodes = NULL ;       // nodes touching an edge in new_members
    GrB_Vector new_members_node_degree = NULL ; // node degree considering only new members
    // GrB_Vector discard = NULL ;              // edges to discard from consideration
    GrB_Vector result = NULL ;                  // resulting matching    

    GrB_Vector empty = NULL ;

    GrB_Index num_edges ;
    GrB_Index num_nodes ;
    
    GRB_TRY (GrB_Matrix_nrows (&num_nodes, E)) ;
    GRB_TRY (GrB_Matrix_ncols (&num_edges, E)) ;
    // TODO: match this type with E
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
    GRB_TRY (GrB_Vector_new (&new_neighbors, GrB_BOOL, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&new_members_nodes, GrB_BOOL, num_nodes)) ;
    GRB_TRY (GrB_Vector_new (&new_members_node_degree, GrB_UINT64, num_nodes)) ;
    // GRB_TRY (GrB_Vector_new (&discard, GrB_BOOL, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&result, GrB_BOOL, num_edges)) ;

 
    LG_TRY (LAGraph_Random_Seed (Seed, seed, msg)) ;

    // initially all edges are considered
    GRB_TRY (GrB_assign (candidates, NULL, NULL, true, GrB_ALL, 
        num_edges, NULL)) ;
     
    GrB_Index ncandidates ;
    GrB_Index nfailures = 0 ; // counts how many iterations have failed due to invalid matchings

    GRB_TRY (GrB_Vector_nvals(&ncandidates, candidates)) ;

    LAGRAPH_TRY (LAGraph_Vector_Print (candidates, LAGraph_SHORT, stdout, msg)) ;

    // TODO: fix semirings to match type of E. For now, let's assume E has integral edge weights.
    // for each node, counts incident edges
    // we disallow edges that are already in the matching from contributing by using candidates vector
    GRB_TRY (GrB_mxv (node_degree, NULL, NULL, LAGraph_plus_one_uint64, E, candidates, NULL)) ;

    // for each edge, sums incident edges for each node. Each edge has an excess of 2 degree, but it doesn't matter since
    // we care about relative degree
    GRB_TRY (GrB_mxv (degree, candidates, NULL, LAGraph_plus_second_uint64, E_t, node_degree, NULL)) ;

    while (ncandidates > 0) {
        // first just generate the scores again
        if (matching_type == 0) {

            // weighs the score based on degree
            // use RS to clear out old scores
            GRB_TRY (GrB_eWiseMult (score, candidates, NULL, GrB_DIV_FP32, Seed, degree, GrB_DESC_RS)) ;

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
        GRB_TRY (GrB_mxv (max_node_neighbor, NULL, NULL, GrB_MAX_SECOND_SEMIRING_FP32, E, score, NULL)) ;

        // Max edge touching each candidate edge, including itself
        GRB_TRY (GrB_mxv (max_neighbor, candidates, NULL, GrB_MAX_SECOND_SEMIRING_FP32, E_t, max_node_neighbor, GrB_DESC_RS)) ;

        GRB_TRY (GrB_eWiseAdd (new_members, candidates, NULL, GrB_GE_FP32, score, max_neighbor, GrB_DESC_R)) ;
        // sparsifies new_members
        GRB_TRY (GrB_select (new_members, NULL, NULL, GrB_VALUEEQ_BOOL, new_members, true, NULL)) ; 

        // check if any node has > 1 edge touching it. 
        GRB_TRY (GrB_mxv (new_members_node_degree, new_members, NULL, LAGraph_plus_second_uint64, E, new_members, GrB_DESC_RS)) ;
        GrB_Index max_degree ; 
        GRB_TRY (GrB_reduce (&max_degree, NULL, GrB_MAX_MONOID_UINT64, new_members_node_degree, NULL)) ;
        if (max_degree >= 2) {
            nfailures++ ;
            if (nfailures > MAX_FAILURES) {
                break ;
            }
            // regen seed
            LG_TRY (LAGraph_Random_Seed (Seed, seed, msg)) ;
            continue ;
        }
        // add new members to result and remove from candidates
        // also want to remove all adjacent edges in new_members from candidates
        GRB_TRY (GrB_assign (result, new_members, NULL, true, GrB_ALL, num_edges, GrB_DESC_S)) ; 
        // to include neighbor edges, need to compute new_neighbors
        GRB_TRY (GrB_mxv (new_members_nodes, NULL, NULL, LAGraph_any_one_bool, E, new_members, NULL)) ;
        GRB_TRY (GrB_mxv (new_neighbors, NULL, NULL, LAGraph_any_one_bool, E_t, new_members_nodes, NULL)) ;

        GRB_TRY (GrB_assign (candidates, new_neighbors, NULL, empty, GrB_ALL, num_edges, GrB_DESC_S)) ;

        GrB_Index last_ncandidates = ncandidates ;
        // do something about stalling?

        ncandidates = GrB_Vector_nvals(&ncandidates, candidates) ;
    }
    (*matching) = result ;

    return (GrB_SUCCESS) ;
}