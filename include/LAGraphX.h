//------------------------------------------------------------------------------
// LAGraphX.h: include file for LAGraph experimental code
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#ifndef LAGRAPHX_H
#define LAGRAPHX_H

#include <GraphBLAS.h>
#include <LAGraph.h>

//==============================================================================
// Experimental methods: in experimental/algorithm and experimental/utility
//==============================================================================

// Do not rely on these in production.  These methods are still under
// development, and is intended only for illustration not benchmarking.  Do not
// use for benchmarking, without asking the authors.

//------------------------------------------------------------------------------
// LAGRAPH_OK: call LAGraph or GraphBLAS and check the result
//------------------------------------------------------------------------------

// To use LAGRAPH_OK, the #include'ing file must declare a scalar GrB_Info
// info, and must define LAGraph_FREE_ALL as a macro that frees all workspace
// if an error occurs.  The method can be a GrB_Info scalar as well, so that
// LAGRAPH_OK(info) works.  The function that uses this macro must return
// GrB_Info, or int.

#define LAGRAPH_ERROR(message,info)                                         \
{                                                                           \
    fprintf (stderr, "LAGraph error: %s\n[%d]\nFile: %s Line: %d\n",        \
        message, info, __FILE__, __LINE__) ;                                \
    LAGraph_FREE_ALL ;                                                      \
    return (info) ;                                                         \
}

#define LAGRAPH_OK(method)                                                  \
{                                                                           \
    info = method ;                                                         \
    if (! (info == GrB_SUCCESS || info == GrB_NO_VALUE))                    \
    {                                                                       \
        LAGRAPH_ERROR ("", info) ;                                          \
    }                                                                       \
}

//****************************************************************************
// Utilities
//****************************************************************************


//****************************************************************************
// Random number generator
//****************************************************************************

int LAGraph_Random_Init (char *msg) ;
int LAGraph_Random_Finalize (char *msg) ;

#if defined ( COVERAGE )
// for testing only
LAGRAPH_PUBLIC bool random_hack ;
#endif

int LAGraph_Random_Seed     // construct a random seed vector
(
    // input/output
    GrB_Vector Seed,    // vector of random number seeds
    // input
    int64_t seed,       // scalar input seed
    char *msg
) ;

int LAGraph_Random_Next     // random int64 vector of seeds
(
    // input/output
    GrB_Vector Seed,
    char *msg
) ;

int LAGraph_Random_INT64    // random int64 vector
(
    // output
    GrB_Vector X,       // already allocated on input
    // input/output
    GrB_Vector Seed,
    char *msg
) ;

GrB_Info LAGraph_Random_FP64    // random double vector
(
    // output
    GrB_Vector X,       // already allocated on input
    // input/output
    GrB_Vector Seed,
    char *msg
) ;

GrB_Info LAGraph_Random_FP32    // random float vector
(
    // output
    GrB_Vector X,       // already allocated on input
    // input/output
    GrB_Vector Seed,
    char *msg
) ;

//****************************************************************************
// Algorithms
//****************************************************************************

//****************************************************************************
/**
 * Given a symmetric graph A with no-self edges, compute all k-trusses of A.
 *
 * @param[out]  Cset    size n, output k-truss subgraphs.
 * @param[out]  kmax    smallest k where k-truss is empty
 * @param[out]  ntris   Array of size n (on input), ntris [k] is num triangles in k-truss
 * @param[out]  nedges  Array of size n (on input), nedges [k] is num edges in k-truss
 * @param[out]  nstepss Array of size n (on input), nstepss [k] is num steps for k-truss
 * @param[in]   G       input graph, A, not modified.  Must be undirected
 *                      or directed with symmetric structure, no self edges.
 *
 * @todo Need change return value to int (not GrB_Info)
 *
 * @retval GrB_SUCCESS      if completed successfully (equal or not)
 * @retval GrB_NULL_POINTER if kmax, ntris, nedges, nsteps is NULL
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 */
int LAGraph_AllKTruss   // compute all k-trusses of a graph
(
    // outputs
    GrB_Matrix *Cset,   // size n, output k-truss subgraphs
    int64_t *kmax,      // smallest k where k-truss is empty
    int64_t *ntris,     // size max(n,4), ntris [k] is #triangles in k-truss
    int64_t *nedges,    // size max(n,4), nedges [k] is #edges in k-truss
    int64_t *nstepss,   // size max(n,4), nstepss [k] is #steps for k-truss
    // input
    LAGraph_Graph G,    // input graph
    char *msg
) ;

//****************************************************************************
/**
 * Given an undirected graph G with no-self edges, LAGraph_KTruss finds the
 * k-truss subgraph of G.
 *
 * @param[out]  C       k-truss subgraph, of type GrB_UINT32
 * @param[in]   G       input graph, not modified
 * @param[in]   k       the truss to find
 *
 * @retval GrB_SUCCESS      if completed successfully (equal or not)
 * @retval GrB_NULL_POINTER if C or C_type is NULL
 * @return Any GraphBLAS errors that may have been encountered
 */
int LAGraph_KTruss      // compute the k-truss of a graph
(
    // outputs:
    GrB_Matrix *C,      // output k-truss subgraph, C
    // inputs:
    LAGraph_Graph G,    // input graph
    uint32_t k,         // find the k-truss, where k >= 3
    char *msg
) ;

//****************************************************************************
// Connected components
//****************************************************************************

/**
 * Determine connected components in an undirected graph.
 *
 * @param[out] result    array of component identifiers for each vertex (allocated
 *                       by the algorithm, ownership returned to caller).
 * @param[in]  A         the graph (symmetric)
 * @param[in]  sanitize  If true, test to ensure A is symmetric
 *
 * @retval GrB_SUCCESS      if completed successfully
 * @retval GrB_NULL_POINTER if result is NULL
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 */
GrB_Info LAGraph_cc_lacc (
    GrB_Vector *result,
    GrB_Matrix A,
    bool sanitize
) ;

//****************************************************************************
// Bellman Ford variants
//****************************************************************************

/**
 * Bellman-Ford single source shortest paths, returning just the shortest path
 * lengths.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[in]   A            matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output or A is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 */
GrB_Info LAGraph_BF_basic
(
    GrB_Vector *pd_output,
    const GrB_Matrix A,
    const GrB_Index s
) ;

/**
 * Bellman-Ford single source shortest paths, returning just the shortest path
 * lengths.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[in]   A            matrix for the graph (optional-ish)
 * @param[in]   AT           transpose of A (optional-ish)
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output is NULL or both A and AT are NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 *
 */
GrB_Info LAGraph_BF_basic_pushpull
(
    GrB_Vector *pd_output,
    const GrB_Matrix A,
    const GrB_Matrix AT,
    const GrB_Index s
) ;

/**
 * Bellman-Ford single source shortest paths, returning just the shortest path
 * lengths.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[in]   AT           transposed adjacency matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output or AT is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 *
 */
GrB_Info LAGraph_BF_basic_mxv
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    const GrB_Matrix AT,        //transposed adjacency matrix for the graph
    const GrB_Index s           //given index of the source
);

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[out]  ppi_output   the pointer to the vector of parent (created internally)
 * @param[out]  ph_output    the pointer to the vector of hops (created internally)
 * @param[in]   A            adjacency matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output, ppi_output, ph_output, or A is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 *
 */
GrB_Info LAGraph_BF_full
(
    GrB_Vector *pd_output,
    GrB_Vector *ppi_output,
    GrB_Vector *ph_output,
    const GrB_Matrix A,
    const GrB_Index s
) ;

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[out]  ppi_output   the pointer to the vector of parent (created internally)
 * @param[out]  ph_output    the pointer to the vector of hops (created internally)
 * @param[in]   A            adjacency matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output, ppi_output, ph_output, or A is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 *
 */
GrB_Info LAGraph_BF_full1
(
    GrB_Vector *pd_output,
    GrB_Vector *ppi_output,
    GrB_Vector *ph_output,
    const GrB_Matrix A,
    const GrB_Index s
);

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[out]  ppi_output   the pointer to the vector of parent (created internally)
 * @param[out]  ph_output    the pointer to the vector of hops (created internally)
 * @param[in]   A            adjacency matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output, ppi_output, ph_output, or A is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 *
 */
GrB_Info LAGraph_BF_full1a
(
    GrB_Vector *pd_output,
    GrB_Vector *ppi_output,
    GrB_Vector *ph_output,
    const GrB_Matrix A,
    const GrB_Index s
);

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[out]  ppi_output   the pointer to the vector of parent (created internally)
 * @param[out]  ph_output    the pointer to the vector of hops (created internally)
 * @param[in]   A            adjacency matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output, ppi_output, ph_output, or A is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 *
 */
GrB_Info LAGraph_BF_full2
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    GrB_Vector *ppi_output,     //the pointer to the vector of parent
    GrB_Vector *ph_output,      //the pointer to the vector of hops
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Index s           //given index of the source
);

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[out]  ppi_output   the pointer to the vector of parent (created internally)
 * @param[out]  ph_output    the pointer to the vector of hops (created internally)
 * @param[in]   AT           transpose of the adjacency matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output, ppi_output, ph_output, or AT is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 *
 */
GrB_Info LAGraph_BF_full_mxv
(
    GrB_Vector *pd_output,
    GrB_Vector *ppi_output,
    GrB_Vector *ph_output,
    const GrB_Matrix AT,
    const GrB_Index s
);

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree (integer weights).
 *
 * @param[out]  pd       pointer to distance vector d, d(k) = shortest distance
 *                       between s and k if k is reachable from s
 * @param[out]  ppi      pointer to parent index vector pi, pi(k) = parent of
 *                       node k in the shortest path tree
 * @param[in]   s        index of the source
 * @param[in]   n        number of nodes
 * @param[in]   nz       number of edges
 * @param[in]   I        row index vector (size n)
 * @param[in]   J        column index vector (size nz)
 * @param[in]   W        weight vector (size nz), W(i) = weight of edge (I(i),J(i))
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd, ppi, I, J, or W is NULL
 * @retval GrB_INVALID_VALUE  if s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails.
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 *
 */
GrB_Info LAGraph_BF_pure_c
(
    int32_t **pd,

    int64_t **ppi,

    const int64_t s,
    const int64_t n,
    const int64_t nz,
    const int64_t *I,
    const int64_t *J,
    const int32_t *W
);

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree (double weights).
 *
 * @param[out]  pd       pointer to distance vector d, d(k) = shortest distance
 *                       between s and k if k is reachable from s
 * @param[out]  ppi      pointer to parent index vector pi, pi(k) = parent of
 *                       node k in the shortest path tree
 * @param[in]   s        index of the source
 * @param[in]   n        number of nodes
 * @param[in]   nz       number of edges
 * @param[in]   I        row index vector (size n)
 * @param[in]   J        column index vector (size nz)
 * @param[in]   W        weight vector (size nz), W(i) = weight of edge (I(i),J(i))
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd, ppi, I, J, or W is NULL
 * @retval GrB_INVALID_VALUE  if s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails.
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 *
 */
GrB_Info LAGraph_BF_pure_c_double
(
    double **pd,

    int64_t **ppi,

    const int64_t s,
    const int64_t n,
    const int64_t nz,
    const int64_t *I,
    const int64_t *J,
    const double  *W
);

//****************************************************************************
/**
 * Community detection using label propagation algorithm
 *
 * @param[out]  CDLP_handle  community vector
 * @param[out]  CDLP_type    type of element stored in the community vector
 * @param[in]   A            adjacency matrix for the graph
 * @param[in]   symmetric    denote whether the matrix is symmetric
 * @param[in]   sanitize     if true, verify that A is binary
 * @param[in]   itermax      max number of iterations (0 computes nothing)
 * @param[out]  t            array of two doubles allocated by caller:
 *                           [0]=sanitize time, [1]=cdlp time in seconds
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If t, CDLP_handle or CDLP_type is NULL
 * @retval GrB_INVALID_OBJECT If A is not stored in CSR format (FIXME)
 * @retval GrB_OUT_OF_MEMORY  if allocation fails.
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 */
GrB_Info LAGraph_cdlp
(
    GrB_Vector *CDLP_handle,
    GrB_Type *CDLP_type,
    const GrB_Matrix A,
    bool symmetric,
    bool sanitize,
    int itermax,
    double *t
);

//****************************************************************************
/**
 * Sparse deep neural network inference. Performs ReLU inference using input
 * feature vectors Y0.
 *
 * @param[out]  Yhandle      Y, created on output
 * @param[in]   W            W [0..nlayers-1], each nneurons-by-nneurons
 * @param[in]   Bias         Bias [0..nlayers-1], diagonal nneurons-by-nneurons
 * @param[in]   nlayers      number of layers
 * @param[in]   Y0           input features: nfeatures-by-nneurons
 *
 * @retval GrB_SUCCESS         if completed successfully
 * @retval GrB_PANIC           vanilla version has not been implemented yet (NOT_IMPLEMENTED?)
 * @retval GrB_NULL_POINTER    If Yhandle, W, Bias, or Y0 is NULL
 * @retval GrB_DOMAIN_MISMATCH if type of Y0 is not FP32 or FP64, or the types of
 *                             W or Bias arent the same as Y0
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 */
GrB_Info LAGraph_dnn
(
    // output
    GrB_Matrix *Yhandle,
    // input: not modified
    GrB_Matrix *W,
    GrB_Matrix *Bias,
    int nlayers,
    GrB_Matrix Y0
);

//****************************************************************************
/**
 * Compute all-pairs shortest paths using Floyd-Warshall method
 *
 * @param[in]   G       input graph, with edge weights
 * @param[out]  D       output graph, created on output
 * @param[out]  D_type  type of scalar stored in D (see source for explanation)
 *
 * @retval GrB_SUCCESS         if completed successfully
 * @retval GrB_PANIC           vanilla version has not been implemented yet (NOT_IMPLEMENTED?)
 * @retval GrB_NULL_POINTER    If D or D_type is NULL
 * @retval GrB_INVALID_VALUE   If G is not square
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 */
GrB_Info LAGraph_FW
(
    const GrB_Matrix G,
    GrB_Matrix *D,
    GrB_Type   *D_type
);

//****************************************************************************
/**
 * Compute the local clustering coefficient for all nodes in a graph.
 *
 * @param[out]  LCC_handle   output vector holding coefficients
 * @param[out]  LCC_type     type scalars stored in LCC
 * @param[in]   A            adjacency matrix for the graph
 * @param[in]   symmetric    denote whether the matrix is symmetric
 * @param[in]   sanitize     if true, verify that A is binary
 * @param[out]  t            array of two doubles
 *                           [0]=sanitize time, [1]=lcc time in seconds
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_PANIC           vanilla version has not been implemented yet (NOT_IMPLEMENTED?)
 * @retval GrB_NULL_POINTER   If LCC_handle or LCC_type is NULL
 * @retval GrB_INVALID_VALUE  If A is not stored in CSR format (FIXME, OBJECT elsewhere)
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 */
GrB_Info LAGraph_lcc
(
    GrB_Vector *LCC_handle,
    GrB_Type   *LCC_type,
    const GrB_Matrix A,
    bool symmetric,
    bool sanitize,
    double t [2]
);

//****************************************************************************
GrB_Info LAGraph_msf (
    GrB_Matrix *result,     // output: an unsymmetrical matrix, the spanning forest
    GrB_Matrix A,           // input matrix
    bool sanitize           // if true, ensure A is symmetric
) ;

//****************************************************************************
GrB_Info LAGraph_scc (
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A            // input matrix
) ;

//****************************************************************************
int LAGraph_VertexCentrality_Triangle       // vertex triangle-centrality
(
    // outputs:
    GrB_Vector *centrality,     // centrality(i): triangle centrality of i
    uint64_t *ntriangles,       // # of triangles in the graph
    // inputs:
    int method,                 // 0, 1, 2, or 3
    LAGraph_Graph G,            // input graph
    char *msg
) ;

//****************************************************************************
int LAGraph_MaximalIndependentSet       // maximal independent set
(
    // outputs:
    GrB_Vector *mis,            // mis(i) = true if i is in the set
    // inputs:
    LAGraph_Graph G,            // input graph
    int64_t seed,               // random number seed
    GrB_Vector ignore_node,     // if NULL, no nodes are ignored.  Otherwise
                                // ignore_node(i) = true if node i is to be
                                // ignored, and not treated as a candidate
                                // added to maximal independent set.
    char *msg
) ;

int LG_CC_FastSV5           // SuiteSparse:GraphBLAS method, with GxB extensions
(
    // output
    GrB_Vector *component,  // output: array of component identifiers
    // inputs
    LAGraph_Graph G,        // input graph, modified then restored
    char *msg
) ;

#endif
