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
// info, and must define LAGRAPH_FREE_ALL as a macro that frees all workspace
// if an error occurs.  The method can be a GrB_Info scalar as well, so that
// LAGRAPH_OK(info) works.  The function that uses this macro must return
// GrB_Info, or int.

#define LAGRAPH_ERROR(message,info)                                         \
{                                                                           \
    fprintf (stderr, "LAGraph error: %s\n[%d]\nFile: %s Line: %d\n",        \
        message, info, __FILE__, __LINE__) ;                                \
    LAGRAPH_FREE_ALL ;                                                      \
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

GrB_Info LAGraph_log
(
    char *caller,           // calling function
    char *message1,         // message to include (may be NULL)
    char *message2,         // message to include (may be NULL)
    int nthreads,           // # of threads used
    double t                // time taken by the test
) ;

// LAGraph_BinRead: read a matrix from a binary file
int LAGraph_binread         // returns 0 if successful, -1 if failure
(
    GrB_Matrix *A,          // matrix to read from the file
    GrB_Type   *A_type,     // type of the scalar stored in A
    FILE       *f           // file to read it from, already open
) ;

GrB_Info LAGraph_tsvread        // returns GrB_SUCCESS if successful
(
    GrB_Matrix *Chandle,        // C, created on output
    FILE *f,                    // file to read from (already open)
    GrB_Type type,              // the type of C to create
    GrB_Index nrows,            // C is nrows-by-ncols
    GrB_Index ncols
) ;

GrB_Info LAGraph_dense_relabel   // relabel sparse IDs to dense row/column indices
(
    GrB_Matrix *Id2index_handle, // output matrix: A(id, index)=1 (unfilled if NULL)
    GrB_Matrix *Index2id_handle, // output matrix: B(index, id)=1 (unfilled if NULL)
    GrB_Vector *id2index_handle, // output vector: v(id)=index (unfilled if NULL)
    const GrB_Index *ids,        // array of unique identifiers (under LAGRAPH_INDEX_MAX)
    GrB_Index nids,              // number of identifiers
    GrB_Index *id_dimension      // number of rows in Id2index matrix, id2index vector (unfilled if NULL)
) ;

GrB_Info LAGraph_Vector_IsEqual_op    // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Vector A,
    GrB_Vector B,
    GrB_BinaryOp userop,    // comparator to use (required)
    char *msg
);

GrB_Info LAGraph_Vector_IsEqual_type // return GrB_SUCCESS if successful
(
    bool         *result,          // true if A == B, false if A != B or error
    GrB_Vector    A,
    GrB_Vector    B,
    GrB_Type      type,            // use GrB_EQ_type operator to compare A and B
    char         *msg
) ;

int LAGraph_Vector_IsEqual         // returns 0 if successful, < 0 if failure
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Vector A,
    GrB_Vector B,
    char *msg
);

//int LAGraph_TSVRread
//(
//    GrB_Matrix *C,             // C, created on output
//    FILE       *f,                    // file to read from (already open)
//    GrB_Type    type,              // the type of C to create
//    GrB_Index   nrows,            // C is nrows-by-ncols
//    GrB_Index   ncols,
//    char       *msg
//) ;

//GrB_Info LAGraph_IsPattern  // return GrB_SUCCESS if successful
//(
//    bool        *result,           // true if A is all one, false otherwise
//    GrB_Matrix   A,
//    GrB_UnaryOp  userop,      // for A with arbitrary user-defined type.
//    char        *msg
//) ;

//****************************************************************************
// Algorithms
//****************************************************************************

GrB_Info LAGraph_allktruss      // compute all k-trusses of a graph
(
    GrB_Matrix *Cset,           // size n, output k-truss subgraphs (optional)
    GrB_Matrix A,               // input adjacency matrix, A, not modified
    // output statistics
    int64_t *kmax,              // smallest k where k-truss is empty
    int64_t *ntris,             // size n, ntris [k] is #triangles in k-truss
    int64_t *nedges,            // size n, nedges [k] is #edges in k-truss
    int64_t *nstepss            // size n, nstepss [k] is #steps for k-truss
) ;

//****************************************************************************
// Connected components
//****************************************************************************

GrB_Info LAGraph_cc_lacc (
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A,           // input matrix
    bool sanitize           // if true, ensure A is symmetric
) ;

GrB_Info LAGraph_cc_boruvka (
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A,           // input matrix
    bool sanitize           // if true, ensure A is symmetric
) ;

//****************************************************************************
// Bellman Ford variants
//****************************************************************************

GrB_Info LAGraph_BF_basic
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Index s           //given index of the source
) ;

GrB_Info LAGraph_BF_basic_pushpull
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Matrix AT,        //transpose of A (optional)
    const GrB_Index s           //given index of the source
) ;

GrB_Info LAGraph_BF_basic_mxv
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    const GrB_Matrix AT,        //transposed adjacency matrix for the graph
    const GrB_Index s           //given index of the source
);

GrB_Info LAGraph_BF_full
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    GrB_Vector *ppi_output,     //the pointer to the vector of parent
    GrB_Vector *ph_output,      //the pointer to the vector of hops
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Index s           //given index of the source
) ;

GrB_Info LAGraph_BF_full1
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    GrB_Vector *ppi_output,     //the pointer to the vector of parent
    GrB_Vector *ph_output,      //the pointer to the vector of hops
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Index s           //given index of the source
);

GrB_Info LAGraph_BF_full1a
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    GrB_Vector *ppi_output,     //the pointer to the vector of parent
    GrB_Vector *ph_output,      //the pointer to the vector of hops
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Index s           //given index of the source
);

GrB_Info LAGraph_BF_full2
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    GrB_Vector *ppi_output,     //the pointer to the vector of parent
    GrB_Vector *ph_output,      //the pointer to the vector of hops
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Index s           //given index of the source
);

GrB_Info LAGraph_BF_full_mxv
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    GrB_Vector *ppi_output,     //the pointer to the vector of parent
    GrB_Vector *ph_output,      //the pointer to the vector of hops
    const GrB_Matrix AT,        //transposed adjacency matrix for the graph
    const GrB_Index s           //given index of the source
);

GrB_Info LAGraph_BF_pure_c
(
    int32_t **pd,     // pointer to distance vector d, d(k) = shorstest distance
                     // between s and k if k is reachable from s
    int64_t **ppi,   // pointer to parent index vector pi, pi(k) = parent of
                     // node k in the shortest path tree
    const int64_t s, // given source node index
    const int64_t n, // number of nodes
    const int64_t nz,// number of edges
    const int64_t *I,// row index vector
    const int64_t *J,// column index vector
    const int32_t *W // weight vector, W(i) = weight of edge (I(i),J(i))
);

GrB_Info LAGraph_BF_pure_c_double
(
    double **pd,     // pointer to distance vector d, d(k) = shorstest distance
                     // between s and k if k is reachable from s
    int64_t **ppi,   // pointer to parent index vector pi, pi(k) = parent of
                     // node k in the shortest path tree
    const int64_t s, // given source node index
    const int64_t n, // number of nodes
    const int64_t nz,// number of edges
    const int64_t *I,// row index vector
    const int64_t *J,// column index vector
    const double  *W // weight vector, W(i) = weight of edge (I(i),J(i))
);

//****************************************************************************

GrB_Info LAGraph_cdlp
(
    GrB_Vector *CDLP_handle, // output vector
    GrB_Type *CDLP_type,     // scalar type of output vector
    const GrB_Matrix A,      // input matrix
    bool symmetric,          // denote whether the matrix is symmetric
    bool sanitize,           // if true, ensure A is binary
    int itermax,             // max number of iterations,
    double *t                // t [0] = sanitize time, t [1] = cdlp time,
                             // in seconds
);

//****************************************************************************

GrB_Info LAGraph_dnn    // returns GrB_SUCCESS if successful
(
    // output
    GrB_Matrix *Yhandle,    // Y, created on output
    // input: not modified
    GrB_Matrix *W,      // W [0..nlayers-1], each nneurons-by-nneurons
    GrB_Matrix *Bias,   // Bias [0..nlayers-1], diagonal nneurons-by-nneurons
    int nlayers,        // # of layers
    GrB_Matrix Y0       // input features: nfeatures-by-nneurons
);

//****************************************************************************

GrB_Info LAGraph_FW
(
    const GrB_Matrix G,     // input graph, with edge weights
    GrB_Matrix *D,          // output graph, created on output
    GrB_Type   *D_type
);

//****************************************************************************

GrB_Info LAGraph_ktruss         // compute the k-truss of a graph
(
    GrB_Matrix *C,              // output: k-truss subgraph, C
    GrB_Type   *C_type,
    const GrB_Matrix A,         // input adjacency matrix, A, not modified
    const uint32_t k,           // find the k-truss, where k >= 3
    int32_t *nsteps             // # of steps taken (ignored if NULL)
) ;

#endif
