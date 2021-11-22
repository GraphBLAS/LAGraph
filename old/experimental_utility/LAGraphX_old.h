

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
// Random number generator
//****************************************************************************

uint64_t LAGraph_rand64 (uint64_t *seed);
double LAGraph_rand_double (uint64_t *seed);

