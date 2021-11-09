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

// Do not rely on these in production.  These methods are still under development,
// and is intended only for illustration not benchmarking.  Do not use for
// benchmarking, without asking the authors.

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
// ascii header prepended to all *.grb files
#define LAGRAPH_BIN_HEADER 512

//****************************************************************************
/**
 * LAGraph_binread: read a matrix from a binary file.
 *
 * @todo document the format
 *
 * @param[out]  A       Matrix to read from the file.  Is allocated by this
 *                      method.
 * @param[out]  A_type  The type of the scalars stored in A.  Only built-in
 *                      types are supported.
 * @param[in]   f       A handle to an open file containing the binary data.
 *
 * @retval  0   If operation finishes successfully
 * @retval -1   For various errors:
 *                  A or f is NULL,
 *                  unsupported scalar type,
 *                  Suitesparse:GraphBLAS 5 is not being used
 *                  any fread error
 *                  out of memory
 *                  any GraphBLAS error
 */
int LAGraph_binread
(
    GrB_Matrix *A,
    GrB_Type   *A_type,
    FILE       *f
) ;

//****************************************************************************
/**
 * LAGraph_tsvread: read a matrix from a tsv file
 *
 * Each line in the file specifies a single entry: i, j, x.
 * The indices i and j are assumed to be one-based.  The dimensions of the
 * matrix must be provided by the caller.  This format is used for matrices at
 * http://graphchallenge.org.  The Matrix Market format is recommended instead;
 * it is more flexible and easier to use, since that format includes the matrix
 * type and size in the file itself.  See LAGraph_mmread and LAGraph_mmwrite.
 *
 * @param[out]  A       Matrix read from the file. It is allocated by this
 *                      method
 * @param[in]   f       A handle to an open file containing the tsv data
 * @param[in]   type    The type of the matrix to create (casting may occur?)
 * @param[in]   nrows   Number of rows to set in the matrix
 * @param[in]   ncols   Number of cols to set in the matrix
 *
 * @retval  0   If operation finishes successfully (GrB_SUCCESS)
 * @return  Various GrB error codes from different issues: null pointer, out
 *          of memory, etc.
 */
GrB_Info LAGraph_tsvread
(
    GrB_Matrix *A,
    FILE       *f,
    GrB_Type    type,
    GrB_Index   nrows,
    GrB_Index   ncols
) ;

//****************************************************************************
/**
 * LAGraph_grread: read a matrix from a binary format based on the Galois graph
 *                 reader format
 *
 * The file format consists of a header, with the following content:
 *      uint64_t version : either 1 or 2.  1: nodes are 2^32, 2: nodes are
 *          64 bit.  This value is returned to the caller, but is otherwise
 *          unused.
 *      uint64_t esize : the size of the edge weight, as sizeof (edgetype).
 *          For example, if the file contains edge weights of type int32_t,
 *          esize is sizeof (int32_t) == 4.  The caller must specify the
 *          corresponding GrB_Type, and its size must match esize.
 *      uint64_t n : the number of node in the graph.  The GrB_Matrix is
 *          n-by-n.  Rectangular matrices are not supported by this format.
 *      uint64_t e : the number of edges in the graph
 *
 * This header is followed by a matrix in CSR format:
 *      Gp : an array of size ((n+1) * sizeof (uint64_t)) bytes, but Gp [0] = 0
 *          does not appear in the file.  This section of the file is thus
 *          (n * sizeof (uint64_t)) bytes in length.
 *      Gj : an array of size (e * sizeof (int32_t)), containing the adjaceny
 *          lists.  Note that the indices are 32 bit, not 64 bit, and thus
 *          this format is limited to graphs with n < 2^32.
 *      Gx : an array of size (e * esize), containing the edge weights.
 *
 * @param[out]  G          Matrix read from the file. It is allocated by this
 *                         method.
 * @param[out]  G_version  The version of the file
 * @param[in]   filename   C-string containing the pathname of the file to open.
 * @param[in]   gtype      The type of the matrix to create. If the are no edge
 *                         weights, this is set to GrB_NULL, and a GrB_BOOL matrix
 *                         is created with all edge weights equal to 'true'.
 *
 * @retval  0   If operation finishes successfully (GrB_SUCCESS)
 * @return  Various GrB error codes from different issues: null pointer, out
 *          of memory, etc.
 */
 GrB_Info LAGraph_grread
(
    GrB_Matrix *G,
    uint64_t *G_version,
    const char *filename,
    GrB_Type gtype
);


//****************************************************************************
/**
 * LAGraph_dense_relabel: relabel sparse IDs to dense row/column indices
 *
 * Converts array of sparse IDs (ids) to row/column indices between 0...(nids-1).
 * The order of IDs is kept, therefore ids can be used for index -> ID
 * conversion: ids[index]=id.
 *
 * Gives back two binary matrices for conversion between ID- and index-based
 * vertices.
 *
 * @param[out]  Id2index_handle  A[id, index] = 1 (unfilled if NULL)
 * @param[out]  Index2id_handle  B[index, id] = 1 (unfilled if NULL)
 * @param[out]  id2index_handle  v[id] = index    (unfilled if NULL)
 * @param[in]   ids              array of unique identifiers (under LAGRAPH_INDEX_MAX)
 * @param[in]   nids             Number of identifiers in ids
 * @param[out]  id_dimension     Number of rows in Id2index and id2index outputs (unfilled if NULL)
 *
 * @retval  0   If operation finishes successfully (GrB_SUCCESS)
 * @return  Various GrB error codes from different issues: null pointer, out
 *          of memory, etc.
 */
GrB_Info LAGraph_dense_relabel
(
    GrB_Matrix *Id2index_handle,
    GrB_Matrix *Index2id_handle,
    GrB_Vector *id2index_handle,
    const GrB_Index *ids,
    GrB_Index nids,
    GrB_Index *id_dimension
) ;

//****************************************************************************
/**
 * Checks if two vectors are identically equal (same size, pattern,
 * and values) according to a user specified comparator op.
 *
 * @note For the standard API, there is no way to determine the type of a
 *       vector.  Checking for the same type requires the GxB_Vector_type
 *       function, which is an extension in SuiteSparse:GraphBLAS.
 * @note If either or both contain NaN's, result will be false
 * @todo Add GrB_Type input parameters to check for type with the standard API
 *
 * @param[out]   result   Set to true on return is vectors are "equal"
 * @param[in]    A        First vector to compare
 * @param[in]    B        Second vector to compare
 * @param[in]    userop   Binary operator to use for the comparisons
 * @param[out]   msg      If an error code is returned, this may hold an error msg.
 *
 * @retval 0     if completed successfully (equal or not)
 * @retval -1001 result or userop is NULL
 * @return Any GraphBLAS errors that may have been encountered
 */
GrB_Info LAGraph_Vector_IsEqual_op    // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Vector A,
    GrB_Vector B,
    GrB_BinaryOp userop,    // comparator to use (required)
    char *msg
);


//****************************************************************************
/**
 * Checks if two vectors are identically equal (same size, pattern,
 * and values) according to an equal operator of a type specified by
 * the user.
 *
 * @note For the standard API, there is no way to determine the type of a
 *       vector. Checking for the same type requires the GxB_Vector_type
 *       function, which is an extension in SuiteSparse:GraphBLAS.
 * @note If either or both contain NaN's, result will be false
 * @todo Add GrB_Type input parameters to check for type with the standard API
 *
 * @param[out]   result   Set to true on return is vectors are "equal"
 * @param[in]    A        First vector to compare
 * @param[in]    B        Second vector to compare
 * @param[in]    type     The type of the GrB_EQ_<type> operator to compare with
 * @param[out]   msg      If an error code is returned, this may hold an error msg.
 *
 * @retval 0     if completed successfully (equal or not)
 * @retval -1001 result or type is NULL
 * @retval -1002 type is not supported
 * @return Any GraphBLAS errors that may have been encountered
 */
GrB_Info LAGraph_Vector_IsEqual_type // return GrB_SUCCESS if successful
(
    bool         *result,          // true if A == B, false if A != B or error
    GrB_Vector    A,
    GrB_Vector    B,
    GrB_Type      type,            // use GrB_EQ_type operator to compare A and B
    char         *msg
) ;


//****************************************************************************
/**
 * Checks if two vectors are identically equal (same size, type (if accessible),
 * pattern, and values) according to an equal operator of a type determined
 * internally.
 *
 * @note For the standard API, there is no way to determine the type of a
 *       vector and GrB_EQ_FP64 is used. Checking for the same type requires the
 *       GxB_Vector_type function, which is an extension in SuiteSparse:GraphBLAS.
 * @note If either or both contain NaN's, result will be false
 * @todo Add GrB_Type input parameters to check for type with the standard API
 *
 * @param[out]   result   Set to true on return is vectors are "equal"
 * @param[in]    A        First vector to compare
 * @param[in]    B        Second vector to compare
 * @param[out]   msg      If an error code is returned, this may hold an error msg.
 *
 * @retval 0     if completed successfully (equal or not)
 * @retval -1001 A, result or type is NULL
 * @retval -1002 type is not supported
 * @return Any GraphBLAS errors that may have been encountered
 */
int LAGraph_Vector_IsEqual         // returns 0 if successful, < 0 if failure
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Vector A,
    GrB_Vector B,
    char *msg
);

//****************************************************************************
/**
 * Return the pattern of a matrix (spones(A) in MATLAB) as a boolean matrix.
 *
 * @note This method cannot handle matrices containing user defined types
 *
 * @param[out]   C        Boolean matrix with the pattern of A
 * @param[in]    A        Input matrx
 * @param[in]    C_type   Type to use for elements stored in C
 *
 * @retval GrB_SUCCESS    if completed successfully (equal or not)
 * @retval ??? LAGraph_OK errors?
 * @return Any GraphBLAS errors that may have been encountered
 */
GrB_Info LAGraph_pattern    // return GrB_SUCCESS if successful
(
    GrB_Matrix *C,          // a boolean matrix with the structure of A
    GrB_Matrix A,
    GrB_Type C_type         // return type for C
) ;

//****************************************************************************
/**
 * Remove all entries from the diagonal of the specified matrix
 *
 * @retval GrB_SUCCESS if successful
 * @return Any LAGraph_OK or GrB errors that may have occured
 */
GrB_Info LAGraph_prune_diag
(
    GrB_Matrix A
) ;

//****************************************************************************
// Random number generator
//****************************************************************************

uint64_t LAGraph_rand64 (uint64_t *seed);
double LAGraph_rand_double (uint64_t *seed);

int LAGraph_Random_Init (char *msg) ;
int LAGraph_Random_Finalize (char *msg) ;

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
// Unused and/or deprecated methods
//****************************************************************************

GrB_Info LAGraph_log
(
    char *caller,           // calling function
    char *message1,         // message to include (may be NULL)
    char *message2,         // message to include (may be NULL)
    int nthreads,           // # of threads used
    double t                // time taken by the test
) ;

GrB_Info LAGraph_1_to_n     // create an integer vector v = 1:n
(
    GrB_Vector *v_handle,   // vector to create
    GrB_Index n             // size of vector to create
);

GrB_Info LAGraph_ispattern  // return GrB_SUCCESS if successful
(
    bool *result,           // true if A is all one, false otherwise
    GrB_Matrix A,
    GrB_UnaryOp userop      // for A with arbitrary user-defined type.
                            // Ignored if A and B are of built-in types
);

GrB_Info LAGraph_isall      // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp op         // GrB_EQ_<type>, for the type of A and B,
                            // to check for equality.  Or use any desired
                            // operator.  The operator should return GrB_BOOL.
);

GrB_Info LAGraph_Vector_isall      // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Vector A,
    GrB_Vector B,
    GrB_BinaryOp op         // GrB_EQ_<type>, for the type of A and B,
                            // to check for equality.  Or use any desired
                            // operator.  The operator should return GrB_BOOL.
);

GrB_Info LAGraph_Vector_to_dense
(
    GrB_Vector *vdense,     // output vector
    GrB_Vector v,           // input vector
    void *id                // pointer to value to fill vdense with
);

//****************************************************************************
// Algorithms
//****************************************************************************

//****************************************************************************
/**
 * Given a symmetric graph A with no-self edges, compute all k-trusses of A.
 *
 * @param[out]  Cset    size n, output k-truss subgraphs (optional, if NULL).
 *                      If not NULL, must contain an array of n GrB_Matrix
 *                      handles.
 * @param[in]   A       input adjacency matrix, A, not modified
 * @param[out]  kmax    smallest k where k-truss is empty
 * @param[out]  ntris   Array of size n (on input), ntris [k] is num triangles in k-truss
 * @param[out]  nedges  Array of size n (on input), nedges [k] is num edges in k-truss
 * @param[out]  nstepss Array of size n (on input), nstepss [k] is num steps for k-truss
 *
 * @todo Need a vanilla version based on GraphBLAS 2.0 spec.
 * @todo Need change return value to int (not GrB_Info)
 *
 * @retval GrB_SUCCESS      if completed successfully (equal or not)
 * @retval GrB_NULL_POINTER if kmax, ntris, nedges, nsteps is NULL
 * @retval GrB_NO_VALUE     vanilla version has not been implemented yet (NOT_IMPLEMENTED?)
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 */
GrB_Info LAGraph_allktruss
(
    GrB_Matrix *Cset,
    GrB_Matrix A,
    // output statistics
    int64_t *kmax,
    int64_t *ntris,
    int64_t *nedges,
    int64_t *nstepss
) ;

//****************************************************************************
/**
 * Given a symmetric graph A with no-self edges, ktruss_graphblas finds the
 * k-truss subgraph of A.
 *
 * @param[out]  C       k-truss subgraph.
 * @param[ou]   C_type  type of elements stored in C.
 * @param[in]   A       input adjacency matrix, A, not modified
 * @param[in]   k       the truss to find
 * @param[out]  nsteps  The number of steps taken (ignored if NULL)
 *
 * @retval GrB_SUCCESS      if completed successfully (equal or not)
 * @retval GrB_NULL_POINTER if C or C_type is NULL
 * @retval GrB_NO_VALUE     vanilla version has not been implemented yet (NOT_IMPLEMENTED?)
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 */
GrB_Info LAGraph_ktruss
(
    GrB_Matrix *C,
    GrB_Type   *C_type,
    const GrB_Matrix A,
    const uint32_t k,
    int32_t *nsteps
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
 * @retval GrB_PANIC        vanilla version has not been implemented yet (NOT_IMPLEMENTED?)
 * @return Any GraphBLAS errors that may have been encountered through LAGRAPH_OK.
 */
GrB_Info LAGraph_cc_boruvka (
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

#endif
