//------------------------------------------------------------------------------
// LAGraph.h: user-visible include file for LAGraph
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// LAGraph is a package of graph algorithms based on GraphBLAS.  GraphBLAS
// defines a set of sparse matrix operations on an extended algebra of
// semirings, using an almost unlimited variety of operators and types.  When
// applied to sparse adjacency matrices, these algebraic operations are
// equivalent to computations on graphs.  GraphBLAS provides a powerful and
// expressive framework creating graph algorithms based on the elegant
// mathematics of sparse matrix operations on a semiring.

// However, GraphBLAS itself does not have graph algorithms.  The purpose of
// LAGraph is to provide a robust, easy-to-use high-performance library of
// graph algorithms that rely on GraphBLAS.

//------------------------------------------------------------------------------

#ifndef LAGRAPH_H
#define LAGRAPH_H

//==============================================================================
// LAGraph version
//==============================================================================

#define LAGRAPH_VERSION_MAJOR 0
#define LAGRAPH_VERSION_MINOR 9
#define LAGRAPH_VERSION_UPDATE 9
#define LAGRAPH_DATE "Feb 6, 2022"

//==============================================================================
// include files
//==============================================================================

#include <time.h>
#include <ctype.h>
#include <limits.h>
#include <GraphBLAS.h>
#if defined ( _OPENMP )
    #include <omp.h>
#endif

//==============================================================================
// GraphBLAS platform specifics
//==============================================================================

// GraphBLAS C API specification, OpenMP, and vanilla vs
// SuiteSparse:GraphBLAS GxB extensions.

#if ( GRB_VERSION < 2 )
    #error "The GraphBLAS library must support the v2.0 C API Specification"
#endif

#if ( _MSC_VER && !__INTEL_COMPILER )
    #ifdef LG_LIBRARY
        // compiling LAGraph itself, exporting symbols to user apps
        #define LAGRAPH_PUBLIC extern __declspec ( dllexport )
    #else
        // compiling the user application, importing symbols from LAGraph
        #define LAGRAPH_PUBLIC extern __declspec ( dllimport )
    #endif
#else
    // for other compilers
    #define LAGRAPH_PUBLIC extern
#endif

#if defined ( __cplusplus )
    // C++ does not have the restrict keyword
    #define LAGRAPH_RESTRICT
#elif ( _MSC_VER && !__INTEL_COMPILER )
    // Microsoft Visual Studio uses __restrict instead of restrict for C
    #define LAGRAPH_RESTRICT __restrict
#else
    // use the restrict keyword for ANSI C99 compilers
    #define LAGRAPH_RESTRICT restrict
#endif

// vanilla vs SuiteSparse:
#if !defined ( LG_VANILLA )
    // by default, set LG_VANILLA to false
    #define LG_VANILLA 0
#endif

#if ( !LG_VANILLA ) && defined ( GxB_SUITESPARSE_GRAPHBLAS )
    // use SuiteSparse, and its GxB* extensions
    #define LG_SUITESPARSE 1
#else
    // use any GraphBLAS library (possibly SuiteSparse) but with no GxB*
    #define LG_SUITESPARSE 0
#endif

// maximum length of the name of a GrB type, including the null-terminator
#if LG_SUITESPARSE
#define LAGRAPH_MAX_NAME_LEN GxB_MAX_NAME_LEN
#else
#define LAGRAPH_MAX_NAME_LEN 128
#endif

//==============================================================================
// LAGraph error handling
//==============================================================================

// Each LAGraph function returns an int.  Negative values indicate an error,
// zero denotes success, and positive values denote success but with some kind
// of algorithm-specific note or warning.  In addition, all LAGraph functions
// have a final parameter that is a pointer to a user-allocated string in which
// an algorithm-specific error message can be returned.  If NULL, no error
// message is returned.  This is not itself an error condition, it just
// indicates that the caller does not need the message returned.  If the
// message string is provided but no error occurs, an empty string is returned.

// LAGRAPH_MSG_LEN: The maximum required length of a message string
#define LAGRAPH_MSG_LEN 256

// For example, the following call computes the breadth-first-search of an
// LAGraph_Graph G, starting at a given source node.  It returns a status of
// zero if it succeeds and a negative value on failure.

/*
    GrB_Vector level, parent ;
    char msg [LAGRAPH_MSG_LEN] ;
    int status = LAGraph_BreadthFirstSearch (&level, &parent, G, src, msg) ;
    if (status < 0)
    {
        printf ("error: %s\n", msg) ;
        // take corrective action ...
    }
*/

//------------------------------------------------------------------------------
// LAGraph error codes:
//------------------------------------------------------------------------------

// LAGraph methods return an integer to denote their status:  zero if they are
// successful (which is the value of GrB_SUCCESS), negative on error, or
// positive for an informational value (such as GrB_NO_VALUE).  Integers in the
// range -999 to 999 are reserved for GraphBLAS GrB_Info return values.

//  successful results:
//  GrB_SUCCESS = 0             // all is well
//  GrB_NO_VALUE = 1            // A(i,j) requested but not there
//  GxB_EXHAUSTED = 2           // iterator is exhausted (SuiteSparse only)

//  errors:
//  GrB_UNINITIALIZED_OBJECT = -1   // object has not been initialized
//  GrB_NULL_POINTER = -2           // input pointer is NULL
//  GrB_INVALID_VALUE = -3          // generic error; some value is bad
//  GrB_INVALID_INDEX = -4          // row or column index is out of bounds
//  GrB_DOMAIN_MISMATCH = -5        // object domains are not compatible
//  GrB_DIMENSION_MISMATCH = -6     // matrix dimensions do not match
//  GrB_OUTPUT_NOT_EMPTY = -7       // output matrix already has values
//  GrB_NOT_IMPLEMENTED = -8        // method not implemented
//  GrB_PANIC = -101                // unknown error
//  GrB_OUT_OF_MEMORY = -102        // out of memory
//  GrB_INSUFFICIENT_SPACE = -103,  // output array not large enough
//  GrB_INVALID_OBJECT = -104       // object is corrupted
//  GrB_INDEX_OUT_OF_BOUNDS = -105  // row or col index out of bounds
//  GrB_EMPTY_OBJECT = -106         // an object does not contain a value

// LAGraph uses the GrB_* error codes in these cases:
//  GrB_INVALID_INDEX: if a source node id is out of range
//  GrB_INVALID_VALUE: if an enum to select a method or option is out of range
//  GrB_NOT_IMPLEMENTED: if a type is not supported, or when SuiteSparse
//      GraphBLAS is required.

// Many LAGraph methods share common error cases, described below.  These
// return values are in the range -1000 to -1999.  Return values of -2000 or
// greater may be used by specific LAGraph methods, which denote errors not in
// this list:

//      LAGRAPH_INVALID_GRAPH:  the input graph is invalid; the details are
//      given in the error msg string returned by the method.

//      LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED: the method requires an undirected
//      graph, or a directed graph with an adjancency matrix that is known to
//      have a symmetric structure.  LAGraph_Property_ASymmetricStructure can
//      be used to determine this property.

//      LAGRAPH_IO_ERROR:  a file input or output method failed, or an input
//      file has an incorrect format that cannot be parsed.

//      LAGRAPH_PROPERTY_MISSING:  some methods require one or more cached
//      properties to be computed before calling them (AT, rowdegree, or
//      coldegree.  Use LAGraph_Property_AT, LAGraph_Property_RowDegree,
//      and/or LAGraph_Property_ColDegree to compute them.

//      LAGRAPH_NO_SELF_EDGES_ALLOWED:  some methods requires that the graph
//      have no self edges, which correspond to the entries on the diagonal of
//      the adjacency matrix.  If the G->ndiag property is nonzero or unknown,
//      this error condition is returned.  Use LAGraph_Property_NDiag to
//      compute G->ndiag, or LAGraph_DeleteDiag to remove all diagonal entries
//      from G->A.

#define LAGRAPH_INVALID_GRAPH                   (-1000)
#define LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED    (-1001)
#define LAGRAPH_IO_ERROR                        (-1002)
#define LAGRAPH_PROPERTY_MISSING                (-1003)
#define LAGRAPH_NO_SELF_EDGES_ALLOWED           (-1004)

//------------------------------------------------------------------------------
// LAGraph_TRY: try an LAGraph method and check for errors
//------------------------------------------------------------------------------

// In a robust application, the return values from each call to LAGraph and
// GraphBLAS should be checked, and corrective action should be taken if an
// error occurs.  The LAGraph_TRY and GrB_TRY macros assist in this effort.

// LAGraph and GraphBLAS are written in C, and so they cannot rely on the
// try/catch mechanism of C++.  To accomplish a similar goal, each LAGraph file
// must #define its own file-specific macro called LAGraph_CATCH.  The typical
// usage of macro is to free any temporary matrices/vectors or workspace when
// an error occurs, and then "throw" the error by returning to the caller.  A
// user application may also #define LAGraph_CATCH and use these macros.

// A typical example of a user function that calls LAGraph might #define
// LAGraph_CATCH as follows.  Suppose workvector is a GrB_vector used for
// computations internal to the mybfs function, and W is a (double *) space
// created by malloc.

#if example_usage_only

    // an example user-defined LAGraph_CATCH macro
    #define LAGraph_CATCH(status)                                   \
    {                                                               \
        /* an LAGraph error has occurred */                         \
        printf ("LAGraph error: (%d): file: %s, line: %d\n%s\n",    \
            status, __FILE__, __LINE__, msg) ;                      \
        /* free any internal workspace and return the status */     \
        GrB_free (*parent) ;                                        \
        GrB_free (workvector) ;                                     \
        LAGraph_Free ((void **) &W) ;                               \
        return (status) ;                                           \
    }

    // an example user function that uses LAGraph_TRY / LAGraph_CATCH
    int mybfs (LAGraph_Graph G, GrB_Vector *parent, int64_t src)
    {
        GrB_Vector workvector = NULL ;
        double *W = NULL ;
        char msg [LAGRAPH_MSG_LEN] ;
        (*parent) = NULL ;
        LAGraph_TRY (LAGraph_BreadthFirstSearch (NULL, parent, G, src, msg)) ;
        // ...
        return (GrB_SUCCESS) ;
    }

#endif

#define LAGraph_TRY(LAGraph_method)             \
{                                               \
    int LAGraph_status = LAGraph_method ;       \
    if (LAGraph_status < GrB_SUCCESS)           \
    {                                           \
        LAGraph_CATCH (LAGraph_status) ;        \
    }                                           \
}

//------------------------------------------------------------------------------
// GrB_TRY: try a GraphBLAS method and check for errors
//------------------------------------------------------------------------------

// LAGraph provides a similar functionality for calling GraphBLAS methods.
// GraphBLAS returns info = 0 (GrB_SUCCESS) or 1 (GrB_NO_VALUE) on success, and
// a value < 0 on failure.  The user application must #define GrB_CATCH to use
// GrB_TRY.  Note that GraphBLAS_info is internal to this macro.  If the
// user application or LAGraph method wants a copy, a statement such as
// info = GraphBLAS_info ; where info is defined outside of this macro.

#define GrB_TRY(GrB_method)                                                  \
{                                                                            \
    GrB_Info GraphBLAS_info = GrB_method ;                                   \
    if (GraphBLAS_info < GrB_SUCCESS)                                        \
    {                                                                        \
        GrB_CATCH (GraphBLAS_info) ;                                         \
    }                                                                        \
}

//==============================================================================
// LAGraph memory management
//==============================================================================

// LAGraph provides wrappers for the malloc/calloc/realloc/free set of memory
// management functions, initialized by LAGraph_Init or LAGraph_Xinit.  By
// default, they are pointers to the ANSI C11 malloc/calloc/realloc/free
// functions.  Unlike all other LAGraph utility functions, these methods do not
// return an int, and do not have a final char *msg parameter.  Instead, they
// closely match the sytax of malloc/calloc/realloc/free.  LAGraph_Calloc and
// LAGraph_free have the same syntax as calloc and free.  LAGraph_Malloc has
// the syntax of calloc instead of malloc.  LAGraph_Realloc is very different
// from realloc, since the ANSI C11 realloc syntax is difficult to use safely.

// Only LAGraph_Malloc_function and LAGraph_Free_function are required.
// LAGraph_Calloc_function may be NULL, in which case LAGraph_Malloc and memset
// are used.  Likewise, LAGraph_Realloc_function may be NULL, in which case
// LAGraph_Malloc, memcpy, and LAGraph_Free are used.

LAGRAPH_PUBLIC void * (* LAGraph_Malloc_function  ) (size_t)         ;
LAGRAPH_PUBLIC void * (* LAGraph_Calloc_function  ) (size_t, size_t) ;
LAGRAPH_PUBLIC void * (* LAGraph_Realloc_function ) (void *, size_t) ;
LAGRAPH_PUBLIC void   (* LAGraph_Free_function    ) (void *)         ;

// LAGraph_Malloc:  allocate a block of memory (wrapper for malloc)
LAGRAPH_PUBLIC
void *LAGraph_Malloc        // returns pointer to allocated block of memory
(                           // or NULL if the allocation fails
    size_t nitems,          // number of items
    size_t size_of_item     // size of each item
) ;

// LAGraph_Calloc:  allocate a block of memory (wrapper for calloc)
LAGRAPH_PUBLIC
void *LAGraph_Calloc        // returns pointer to allocated block of memory
(                           // or NULL if the allocation fails
    size_t nitems,          // number of items
    size_t size_of_item     // size of each item
) ;

LAGRAPH_PUBLIC
void *LAGraph_Realloc       // returns pointer to reallocated block of memory,
(                           // or original block if reallocation fails.
    size_t nitems_new,      // new number of items in the object
    size_t nitems_old,      // old number of items in the object
    size_t size_of_item,    // size of each item
    // input/output
    void *p,                // old block to reallocate
    // output
    bool *ok                // true if successful, false otherwise
) ;

// LAGraph_Free:  free a block of memory (wrapper for free)
LAGRAPH_PUBLIC
void LAGraph_Free           // free a block of memory and set p to NULL
(
    void **p                // pointer to object to free, does nothing if NULL
) ;

//==============================================================================
// LAGraph data structures
//==============================================================================

// In addition to relying on underlying GraphBLAS objects (GrB_Matrix,
// GrB_Vector, GrB_Descriptor, ...), LAGraph adds the LAGraph_Graph.  This
// object contains a representation of a graph and its associated data.

// LAGRAPH_UNKNOWN is used for all scalars whose value is not known
#define LAGRAPH_UNKNOWN (-1)

//------------------------------------------------------------------------------
// LAGraph_Kind: the kind of a graph
//------------------------------------------------------------------------------

// Currently, only two types of graphs are supported: undirected graphs and
// directed graphs.  Edge weights are assumed to be present.  Unweighted graphs
// can be represented by setting all entries present in the sparsity structure
// to the same value, typically 1.  Additional types of graphs will be added in
// the future.

typedef enum
{
    LAGRAPH_ADJACENCY_UNDIRECTED = 0, // A is square and symmetric; both upper
                                      // and lower triangular parts are
                                      // present.  A(i,j) is the edge (i,j)

    LAGRAPH_ADJACENCY_DIRECTED = 1,   // A is square; A(i,j) is the edge (i,j)

    // possible future kinds of graphs:
    // LAGRAPH_ADJACENCY_UNDIRECTED_UNWEIGHTED
    // LAGRAPH_ADJACENCY_DIRECTED_UNWEIGHTED
    // LAGRAPH_ADJACENCY_UNDIRECTED_TRIL
    // LAGRAPH_ADJACENCY_UNDIRECTED_TRIU
    // LAGRAPH_BIPARTITE
    // LAGRAPH_BIPARTITE_DIRECTED
    // LAGRAPH_BIPARTITE_UNDIRECTED
    // LAGRAPH_INCIDENCE_*
    // LAGRAPH_MULTIGRAPH_*
    // LAGRAPH_HYPERGRAPH
    // LAGRAPH_HYPERGRAPH_DIRECTED
    // ...

    LAGRAPH_KIND_UNKNOWN = LAGRAPH_UNKNOWN      // the graph kind is unknown
}
LAGraph_Kind ;

//------------------------------------------------------------------------------
// LAGraph_BooleanProperty: true, false, or unknown
//------------------------------------------------------------------------------

typedef enum
{
    LAGRAPH_FALSE = 0,
    LAGRAPH_TRUE = 1,
    LAGRAPH_BOOLEAN_UNKNOWN = LAGRAPH_UNKNOWN
}
LAGraph_BooleanProperty ;

//------------------------------------------------------------------------------
// LAGraph_Graph: the primary graph data structure of LAGraph
//------------------------------------------------------------------------------

// The LAGraph_Graph object contains a GrB_Matrix A as its primary component.
// For graphs represented with adjacency matrices, A(i,j) denotes the edge
// (i,j).  Unlike GrB_* objects in GraphBLAS, the LAGraph_Graph data structure
// is not opaque.  User applications have full access to its contents.

// An LAGraph_Graph G contains two kinds of components:

// (1) primary components of the graph, which fully define the graph:
//      A           the adjacency matrix of the graph
//      kind        the kind of graph (undirected, directed, bipartite, ...)

// (2) cached properties of the graph, which can be recreated any time:
//      AT          AT = A'
//      rowdegree   rowdegree(i) = # of entries in A(i,:)
//      coldegree   coldegree(j) = # of entries in A(:,j)
//      A_structure_is_symmetric: true if the structure of A is symmetric
//      ndiag       the number of entries on the diagonal of A

struct LAGraph_Graph_struct
{

    //--------------------------------------------------------------------------
    // primary components of the graph
    //--------------------------------------------------------------------------

    GrB_Matrix   A;          // the adjacency matrix of the graph
    LAGraph_Kind kind;       // the kind of graph

    // possible future components:
    // multigraph ..
    // GrB_Matrix *Amult ; // array of size nmatrices
    // int nmatrices ;
    // GrB_Vector VertexWeights ;

    //--------------------------------------------------------------------------
    // cached properties of the graph
    //--------------------------------------------------------------------------

    // All of these components may be deleted or set to 'unknown' at any time.
    // For example, if AT is NULL, then the transpose of A has not been
    // computed.  A scalar property of type LAGraph_BooleanProperty would be
    // set to LAGRAPH_UNKNOWN to denote that its value is unknown.

    // If present, the properties must be valid and accurate.  If the graph
    // changes, these properties can either be recomputed or deleted to denote
    // the fact that they are unknown.  This choice is up to individual LAGraph
    // methods and utilities.

    // LAGraph methods can set non-scalar properties only if they are
    // constructing the graph.  They cannot modify them or create them if the
    // graph is declared as a read-only object in the parameter list of the
    // method.

    GrB_Matrix AT;          // AT = A', the transpose of A

    GrB_Vector rowdegree;   // a GrB_INT64 vector of length m, if A is m-by-n.
           // where rowdegree(i) is the number of entries in A(i,:).  If
           // rowdegree is sparse and the entry rowdegree(i) is not present,
           // then it is assumed to be zero.

    GrB_Vector coldegree ;  // a GrB_INT64 vector of length n, if A is m-by-n.
            // where coldegree(j) is the number of entries in A(:,j).  If
            // coldegree is sparse and the entry coldegree(j) is not present,
            // then it is assumed to be zero.  If A is known to have a
            // symmetric structure, the convention is that the degree is held in
            // rowdegree, and coldegree is left as NULL.

    LAGraph_BooleanProperty A_structure_is_symmetric ;    // For an undirected
            // graph, this property will always be implicitly true and can be
            // ignored.  The matrix A for a directed weighted graph will
            // typically be unsymmetric, but might have a symmetric structure.
            // In that case, this scalar property can be set to true.

    int64_t ndiag ; // # of entries on the diagonal of A, or LAGRAPH_UNKNOWN if
            // unknown.  For the adjacency matrix of a directed or undirected
            // graph, this is the number of self-edges in the graph.

    // possible future cached properties:
    // GrB_Vector rowsum, colsum ;
    // rowsum (i) = sum (A (i,:)), regardless of kind
    // colsum (j) = sum (A (:,j)), regardless of kind
    // LAGraph_BooleanProperty connected ;   // true if G is a connected graph
} ;

typedef struct LAGraph_Graph_struct *LAGraph_Graph ;

//==============================================================================
// LAGraph utilities
//==============================================================================

// LAGraph_Init: start GraphBLAS and LAGraph
LAGRAPH_PUBLIC
int LAGraph_Init (char *msg) ;

// LAGraph_Xinit: start GraphBLAS and LAGraph, and set malloc/etc functions
LAGRAPH_PUBLIC
int LAGraph_Xinit
(
    // pointers to memory management functions
    void * (* user_malloc_function  ) (size_t),
    void * (* user_calloc_function  ) (size_t, size_t),
    void * (* user_realloc_function ) (void *, size_t),
    void   (* user_free_function    ) (void *),
    char *msg
) ;

// LAGraph user-defined semirings, created by LAGraph_Init or LAGraph_Xinit:
LAGRAPH_PUBLIC GrB_Semiring

    // LAGraph_plus_first_T: using the GrB_PLUS_MONOID_T monoid and the
    // corresponding GrB_FIRST_T multiplicative operator.
    LAGraph_plus_first_int8   ,
    LAGraph_plus_first_int16  ,
    LAGraph_plus_first_int32  ,
    LAGraph_plus_first_int64  ,
    LAGraph_plus_first_uint8  ,
    LAGraph_plus_first_uint16 ,
    LAGraph_plus_first_uint32 ,
    LAGraph_plus_first_uint64 ,
    LAGraph_plus_first_fp32   ,
    LAGraph_plus_first_fp64   ,

    // LAGraph_plus_second_T: using the GrB_PLUS_MONOID_T monoid and the
    // corresponding GrB_SECOND_T multiplicative operator.
    LAGraph_plus_second_int8   ,
    LAGraph_plus_second_int16  ,
    LAGraph_plus_second_int32  ,
    LAGraph_plus_second_int64  ,
    LAGraph_plus_second_uint8  ,
    LAGraph_plus_second_uint16 ,
    LAGraph_plus_second_uint32 ,
    LAGraph_plus_second_uint64 ,
    LAGraph_plus_second_fp32   ,
    LAGraph_plus_second_fp64   ,

    // LAGraph_plus_one_T: using the GrB_PLUS_MONOID_T monoid and the
    // corresponding GrB_ONEB_T multiplicative operator.
    LAGraph_plus_one_int8   ,
    LAGraph_plus_one_int16  ,
    LAGraph_plus_one_int32  ,
    LAGraph_plus_one_int64  ,
    LAGraph_plus_one_uint8  ,
    LAGraph_plus_one_uint16 ,
    LAGraph_plus_one_uint32 ,
    LAGraph_plus_one_uint64 ,
    LAGraph_plus_one_fp32   ,
    LAGraph_plus_one_fp64   ,

    // LAGraph_structural_T: using the GrB_MIN_MONOID_T for non-boolean types
    // or GrB_LOR_MONOID_BOOL for boolean, and the GrB_ONEB_T multiplicative
    // op.  These semirings are very useful for unweighted graphs, or for
    // algorithms that operate only on the sparsity structure of unweighted
    // graphs.
    LAGraph_structural_bool   ,
    LAGraph_structural_int8   ,
    LAGraph_structural_int16  ,
    LAGraph_structural_int32  ,
    LAGraph_structural_int64  ,
    LAGraph_structural_uint8  ,
    LAGraph_structural_uint16 ,
    LAGraph_structural_uint32 ,
    LAGraph_structural_uint64 ,
    LAGraph_structural_fp32   ,
    LAGraph_structural_fp64   ;

// LAGraph_Finalize: finish LAGraph
LAGRAPH_PUBLIC
int LAGraph_Finalize (char *msg) ;

// LAGraph_MIN/MAX: suitable for integers, and non-NaN floating point
#define LAGraph_MIN(x,y) (((x) < (y)) ? (x) : (y))
#define LAGraph_MAX(x,y) (((x) > (y)) ? (x) : (y))

// LAGraph_New: create a new graph
LAGRAPH_PUBLIC
int LAGraph_New
(
    LAGraph_Graph *G,   // the graph to create, NULL if failure
    GrB_Matrix    *A,   // the adjacency matrix of the graph, may be NULL.
                        // A is moved into G as G->A, and A itself is set
                        // to NULL to denote that is now a part of G.
                        // That is, { G->A = A ; A = NULL ; } is performed.
                        // When G is deleted, G->A is freed.  If A is NULL,
                        // the graph is invalid until G->A is set.
    LAGraph_Kind kind,  // the kind of graph. This may be LAGRAPH_UNKNOWN,
                        // which must then be revised later before the
                        // graph is used.
    char *msg
) ;

// LAGraph_Delete: free a graph and all its contents
LAGRAPH_PUBLIC
int LAGraph_Delete
(
    LAGraph_Graph *G,   // the graph to delete; G set to NULL on output.
                        // All internal GrB_Matrix and GrB_Vector objects are
                        // freed, including G->A.  To keep G->A while deleting
                        // the graph G, use:
                        // { A = G->A ; G->A = NULL ; LAGraph_Delete (&G, msg);}
    char *msg
) ;

// LAGraph_DeleteProperties: free any internal cached properties of a graph
LAGRAPH_PUBLIC
int LAGraph_DeleteProperties
(
    LAGraph_Graph G,    // G stays valid, only cached properties are freed
    char *msg
) ;

// LAGraph_Property_AT: construct G->AT for a graph
LAGRAPH_PUBLIC
int LAGraph_Property_AT
(
    LAGraph_Graph G,    // graph for which to compute G->AT
    char *msg
) ;

// LAGraph_Property_ASymmetricStructure: determine G->A_structure_is_symmetric
LAGRAPH_PUBLIC
int LAGraph_Property_ASymmetricStructure
(
    LAGraph_Graph G,    // graph to determine the symmetry of structure of A
    char *msg
) ;

// LAGraph_Property_RowDegree: determine G->rowdegree
LAGRAPH_PUBLIC
int LAGraph_Property_RowDegree
(
    LAGraph_Graph G,    // graph to determine G->rowdegree
    char *msg
) ;

// LAGraph_Property_ColDegree: determine G->coldegree
LAGRAPH_PUBLIC
int LAGraph_Property_ColDegree
(
    LAGraph_Graph G,    // graph to determine G->coldegree
    char *msg
) ;

// LAGraph_Property_NDiag: determine G->ndiag
LAGRAPH_PUBLIC
int LAGraph_Property_NDiag
(
    LAGraph_Graph G,    // graph to compute G->ndiag
    char *msg
) ;

//  LAGraph_DeleteDiag: remove all diagonal entries fromG->A
LAGRAPH_PUBLIC
int LAGraph_DeleteDiag
(
    LAGraph_Graph G,    // diagonal entries removed, most properties cleared
    char *msg
) ;

// LAGraph_CheckGraph: determine if a graph is valid
LAGRAPH_PUBLIC
int LAGraph_CheckGraph
(
    LAGraph_Graph G,    // graph to check
    // TODO: int level, // 0:quick, O(1), 1:a bit more, 2: still more, 3:
                        // exhaustive!
    char *msg
) ;

// LAGraph_GetNumThreads: determine # of OpenMP threads to use
LAGRAPH_PUBLIC
int LAGraph_GetNumThreads
(
    int *nthreads,      // # of threads to use
    char *msg
) ;

// LAGraph_SetNumThreads: set # of OpenMP threads to use
LAGRAPH_PUBLIC
int LAGraph_SetNumThreads
(
    int nthreads,       // # of threads to use
    char *msg
) ;

// LAGraph_Tic: start the timer
LAGRAPH_PUBLIC
int LAGraph_Tic
(
    double tic [2],     // tic [0]: seconds, tic [1]: nanoseconds
    char *msg
) ;

// LAGraph_Toc: return time since last call to LAGraph_Tic
LAGRAPH_PUBLIC
int LAGraph_Toc
(
    double *t,              // time since last call to LAGraph_Tic
    const double tic [2],   // tic from last call to LAGraph_Tic
    char *msg
) ;


/****************************************************************************
 *
 * LAGraph_MMRead:  read a matrix from a Matrix Market file.
 *
 * The file format used here is compatible with all variations of the Matrix
 * Market "coordinate" and "array" format (http://www.nist.gov/MatrixMarket),
 * for sparse and dense matrices repsectively.  The format is fully described
 * in LAGraph/Doc/MatrixMarket.pdf, and summarized here (with extensions for
 * LAGraph).
 *
 * The first line of the file starts with %%MatrixMarket, with the following
 * format:
 *
 *      %%MatrixMarket matrix <fmt> <type> <storage>
 *
 *      <fmt> is one of: coordinate or array.  The former is a sparse matrix in
 *      triplet form.  The latter is a dense matrix in column-major form.
 *      Both formats are returned as a GrB_Matrix.
 *
 *      <type> is one of: real, complex, pattern, or integer.  The real,
 *      integer, and pattern formats are returned as GrB_FP64, GrB_INT64, and
 *      GrB_BOOL, respectively, but these types are modified by the %%GraphBLAS
 *      structured comment described below.  Complex matrices are currently not
 *      supported.
 *
 *      <storage> is one of: general, Hermitian, symmetric, or skew-symmetric.
 *      The Matrix Market format is case-insensitive, so "hermitian" and
 *      "Hermitian" are treated the same).
 *
 *      Not all combinations are permitted.  Only the following are meaningful:
 *
 *      (1) (coordinate or array) x (real, integer, or complex)
 *          x (general, symmetric, or skew-symmetric)
 *
 *      (2) (coordinate or array) x (complex) x (Hermitian)
 *
 *      (3) (coodinate) x (pattern) x (general or symmetric)
 *
 * The second line is an optional extension to the Matrix Market format:
 *
 *      %%GraphBLAS type <entrytype>
 *
 *      <entrytype> is one of the 11 built-in types (bool, int8_t, int16_t,
 *      int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t, float, or
 *      double.
 *
 *      If this second line is included, it overrides the default GraphBLAS
 *      types for the Matrix Market <type> on line one of the file: real,
 *      pattern, and integer.  The Matrix Market complex <type> is not yet
 *      supported.
 *
 * Any other lines starting with "%" are treated as comments, and are ignored.
 * Comments may be interspersed throughout the file.  Blank lines are ignored.
 * The Matrix Market header is optional in this routine (it is not optional in
 * the Matrix Market format).  If not present, the <fmt> defaults to
 * coordinate, <type> defaults to real, and <storage> defaults to general.  The
 * remaining lines are space delimited, and free format (one or more spaces can
 * appear, and each field has arbitrary width).
 *
 * The Matrix Market file <fmt> can be coordinate or array:
 *
 *      coordinate:  for this format, the first non-comment line must appear,
 *      and it must contain three integers:
 *
 *              nrows ncols nvals
 *
 *          For example, a 5-by-12 matrix with 42 entries would have:
 *
 *              5 12 42
 *
 *          Each of the remaining lines defines one entry.  The order is
 *          arbitrary.  If the Matrix Market <type> is real or integer, each
 *          line contains three numbers: row index, column index, and value.
 *          For example, if A(3,4) is equal to 5.77, a line:
 *
 *              3 4 5.77
 *
 *          would appear in the file.  The indices in the Matrix Market are
 *          1-based, so this entry becomes A(2,3) in the GrB_Matrix returned to
 *          the caller.  If the <type> is pattern, then only the row and column
 *          index appears.  If <type> is complex, four values appear.  If
 *          A(8,4) has a real part of 6.2 and an imaginary part of -9.3, then
 *          the line is:
 *
 *              8 4 6.2 -9.3
 *
 *          and since the file is 1-based but a GraphBLAS matrix is always
 *          0-based, one is subtracted from the row and column indices in the
 *          file, so this entry becomes A(7,3).
 *
 *      array: for this format, the first non-comment line must appear, and it
 *      must contain just two integers:
 *
 *              nrows ncols
 *
 *          A 5-by-12 matrix would thus have the line
 *
 *              5 12
 *
 *          Each of the remaining lines defines one entry, in column major
 *          order.  If the <type> is real or integer, this is the value of the
 *          entry.  An entry if <type> of complex consists of two values, the
 *          real and imaginary part (not yet supported).  The <type> cannot be
 *          pattern in this case.
 *
 *      For both coordinate and array formats, real and complex values may use
 *      the terms INF, +INF, -INF, and NAN to represent floating-point infinity
 *      and NaN values, in either upper or lower case.
 *
 * The <storage> token is general, symmetric, skew-symmetric, or Hermitian:
 *
 *      general: the matrix has no symmetry properties (or at least none that
 *      were exploited when the file was created).
 *
 *      symmetric:  A(i,j) == A(j,i).  Only entries on or below the diagonal
 *      appear in the file.  Each off-diagonal entry in the file creates two
 *      entries in the GrB_Matrix that is returned.
 *
 *      skew-symmetric:  A(i,j) == -A(i,j).  There are no entries on the
 *      diagonal.  Only entries below the diagonal appear in the file.  Each
 *      off-diagonal entry in the file creates two entries in the GrB_Matrix
 *      that is returned.
 *
 *      Hermitian: square complex matrix with A(i,j) = conj (A(j,i)).  All
 *      entries on the diagonal are real.  Each off-diagonal entry in the file
 *      creates two entries in the GrB_Matrix that is returned.
 *
 * According to the Matrix Market format, entries are always listed in
 * column-major order.  This rule is follwed by LAGraph_MMWrite.  However,
 * LAGraph_MMRead can read the entries in any order.
 *
 * @param[out]    A     handle of the matrix to create
 * @param[in]     f     handle to an open file to read from
 * @param[in,out] msg   any error messages
 *
 * @retval GrB_SUCCESS if successful,
 * @retval LAGRAPH_IO_ERROR if the file could not
 *         be read or contains a matrix with an invalid format
 * @retval GrB_NOT_IMPLEMENTED if the type is not supported
 */

// LAGraph_MMRead: read a matrix in MatrixMarket format
LAGRAPH_PUBLIC
int LAGraph_MMRead
(
    GrB_Matrix *A,  // handle of matrix to create
    FILE *f,        // file to read from, already open
    char *msg
);

// LAGraph_MMWrite: write a matrix in MatrixMarket format, auto select type
LAGRAPH_PUBLIC
int LAGraph_MMWrite
(
    GrB_Matrix A,       // matrix to write to the file
    FILE *f,            // file to write it to, must be already open
    FILE *fcomments,    // optional file with extra comments, may be NULL
    char *msg
) ;

// LAGraph_Structure: return the structure of a matrix
LAGRAPH_PUBLIC
int LAGraph_Structure
(
    GrB_Matrix *C,  // a boolean matrix with same structure of A, with C(i,j)
                    // set to true if A(i,j) appears in the sparsity structure
                    // of A.
    GrB_Matrix A,
    char *msg
) ;

// LAGraph_Vector_Structure: return the structure of a vector (TODO)
int LAGraph_Vector_Structure
(
    GrB_Vector *w,  // a boolean vector with same structure of u, with w(i)
                    // set to true if u(i) appears in the sparsity structure
                    // of u.
    GrB_Vector u,
    char *msg
) ;

// LAGraph_NameOfType: return the name of a type
LAGRAPH_PUBLIC
int LAGraph_NameOfType
(
    char *name,     // name of the type: user provided array of size at
                    // least LAGRAPH_MAX_NAME_LEN.
    GrB_Type type,  // GraphBLAS type
    char *msg
) ;

// LAGraph_TypeFromName: return a GrB_Type from its name
LAGRAPH_PUBLIC
int LAGraph_TypeFromName
(
    GrB_Type *type, // GraphBLAS type
    char *name,     // name of the type
    char *msg
) ;

// LAGraph_SizeOfType: return sizeof(...) of a GraphBLAS GrB_Type
LAGRAPH_PUBLIC
int LAGraph_SizeOfType
(
    size_t *size,   // size of the type
    GrB_Type type,  // GraphBLAS type
    char *msg
) ;

// LAGraph_MatrixTypeName: return the name of the GrB_Type of a GrB_Matrix
LAGRAPH_PUBLIC
int LAGraph_MatrixTypeName
(
    char *name,     // name of the type of the matrix A (user-provided array
                    // of size at least LAGRAPH_MAX_NAME_LEN).
    GrB_Matrix A,   // matrix to query
    char *msg
) ;

// LAGraph_VectorTypeName: return the name of the GrB_Type of a GrB_Vector
LAGRAPH_PUBLIC
int LAGraph_VectorTypeName
(
    char *name,     // name of the type of the vector v (user-provided array
                    // of size at least LAGRAPH_MAX_NAME_LEN).
    GrB_Vector v,   // vector to query
    char *msg
) ;

// LAGraph_ScalarTypeName: return the name of the GrB_Type of a GrB_Scalar
LAGRAPH_PUBLIC
int LAGraph_ScalarTypeName
(
    char *name,     // name of the type of the scalar s (user-provided array
                    // of size at least LAGRAPH_MAX_NAME_LEN).
    GrB_Scalar s,   // scalar to query
    char *msg
) ;

// LAGraph_KindName: return the name of a kind
LAGRAPH_PUBLIC
int LAGraph_KindName
(
    char *name,     // name of the kind (user provided array of size at least
                    // LAGRAPH_MAX_NAME_LEN)
    LAGraph_Kind kind,  // graph kind
    char *msg
) ;

// LAGraph_SortByDegree: sort a graph by its row or column degree
LAGRAPH_PUBLIC
int LAGraph_SortByDegree
(
    // output
    int64_t **P_handle,     // P is returned as a permutation vector of size n
    // input
    LAGraph_Graph G,        // graph of n nodes
    bool byrow,             // if true, sort G->rowdegree, else G->coldegree
    bool ascending,         // sort in ascending or descending order
    char *msg
) ;

// LAGraph_SampleDegree: sample the degree median and mean
LAGRAPH_PUBLIC
int LAGraph_SampleDegree
(
    double *sample_mean,    // sampled mean degree
    double *sample_median,  // sampled median degree
    // input
    LAGraph_Graph G,        // graph of n nodes
    bool byrow,             // if true, sample G->rowdegree, else G->coldegree
    int64_t nsamples,       // number of samples
    uint64_t seed,          // random number seed
    char *msg
) ;

// LAGraph_DisplayGraph: print the contents of a graph
LAGRAPH_PUBLIC
int LAGraph_DisplayGraph
(
    LAGraph_Graph G,        // graph to display
    // TODO: use an enum for pr
    int pr,                 // 0: nothing, 1: terse, 2: summary, 3: all,
                            // 4: same as 2 but with %0.15g for doubles
                            // 5: same as 3 but with %0.15g for doubles
    FILE *f,                // file to write to, must already be open
    char *msg
) ;

// LAGraph_IsEqual: compare for exact equality, auto selection of type
LAGRAPH_PUBLIC
int LAGraph_IsEqual     // TODO rename LAGraph_Matrix_IsEqual
(
    bool *result,       // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    char *msg
) ;

// TODO add LAGraph_Matrix_IsEqual_op

//****************************************************************************
/**
 * Checks if two vectors are identically equal (same size, pattern,
 * and values) according to a user specified comparator op.
 *
 * @note If either or both contain NaN's, result will be false
 *
 * @param[out]   result   Set to true on return is vectors are "equal"
 * @param[in]    A        First vector to compare
 * @param[in]    B        Second vector to compare
 * @param[in]    userop   Binary operator to use for the comparisons
 * @param[out]   msg      If an error code is returned, this may hold an error msg.
 *
 * @retval GrB_SUCCESS          if completed successfully (equal or not)
 * @retval GrB_NULL_POINTER     result or userop is NULL
 * @return Any GraphBLAS errors that may have been encountered
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_Vector_IsEqual_op
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Vector A,
    GrB_Vector B,
    GrB_BinaryOp userop,    // comparator to use (required)
    char *msg
);

//****************************************************************************
/**
 * Checks if two vectors are identically equal (same size, type (if accessible),
 * pattern, and values) according to an equal operator of a type determined
 * internally.
 *
 * @note If either or both contain NaN's, result will be false
 *
 * @param[out]   result   Set to true on return is vectors are "equal"
 * @param[in]    A        First vector to compare
 * @param[in]    B        Second vector to compare
 * @param[out]   msg      If an error code is returned, this may hold an error msg.
 *
 * @retval GrB_SUCCESS          if completed successfully (equal or not)
 * @retval GrB_NULL_POINTER     A, result or type is NULL
 * @retval GrB_NOT_IMPLEMENTED  type is not supported
 * @return Any GraphBLAS errors that may have been encountered
 */
LAGRAPH_PUBLIC
int LAGraph_Vector_IsEqual
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Vector A,
    GrB_Vector B,
    char *msg
);

// LAGraph_Matrix_print: pretty-print a matrix
LAGRAPH_PUBLIC
int LAGraph_Matrix_print
(
    GrB_Matrix A,       // matrix to pretty-print to the file
    // TODO: use an enum for pr
    int pr,             // print level: -1 nothing, 0: one line, 1: terse,
                        //      2: summary, 3: all,
                        //      4: as 2 but with %0.15g for float/double
                        //      5: as 3 but with %0.15g for float/double
    FILE *f,            // file to write it to, must be already open; use
                        // stdout or stderr to print to those locations.
    char *msg
) ;

// LAGraph_Vector_print: pretty-print a matrix
LAGRAPH_PUBLIC
int LAGraph_Vector_print
(
    GrB_Vector v,       // vector to pretty-print to the file
    // TODO: use an enum for pr
    int pr,             // print level: -1 nothing, 0: one line, 1: terse,
                        //      2: summary, 3: all,
                        //      4: as 2 but with %0.15g for float/double
                        //      5: as 3 but with %0.15g for float/double
    FILE *f,            // file to write it to, must be already open; use
                        // stdout or stderr to print to those locations.
    char *msg
) ;

//------------------------------------------------------------------------------
// parallel sorting methods
//------------------------------------------------------------------------------

// LAGraph_Sort1: sort array A of size n
LAGRAPH_PUBLIC
int LAGraph_Sort1
(
    int64_t *A_0,       // size n array
    const int64_t n,
    int nthreads,       // # of threads to use
    char *msg
) ;

// LAGraph_Sort2: sort array A of size 2-by-n, using 2 keys (A [0:1][])
LAGRAPH_PUBLIC
int LAGraph_Sort2
(
    int64_t *A_0,       // size n array
    int64_t *A_1,       // size n array
    const int64_t n,
    int nthreads,       // # of threads to use
    char *msg
) ;

// LAGraph_Sort3: sort array A of size 3-by-n, using 3 keys (A [0:2][])
LAGRAPH_PUBLIC
int LAGraph_Sort3
(
    int64_t *A_0,       // size n array
    int64_t *A_1,       // size n array
    int64_t *A_2,       // size n array
    const int64_t n,
    int nthreads,       // # of threads to use
    char *msg
) ;

//==============================================================================
// LAGraph Basic algorithms
//==============================================================================

// Basic algorithm are meant to be easy to use.  They may encompass many
// underlying Advanced algorithms, each with various parameters that may be
// controlled.  For the Basic API, these parameters are determined
// automatically.  Graph properties may be determined, and as a result, the
// graph G is both an input and an output of these methods, since they may be
// modified.

/****************************************************************************
 *
 * Perform breadth-first traversal, computing parent vertex ID's
 * and/or level encountered.
 *
 * @param[out]    level      If non-NULL on input, on successful return, it
 *                           contains the levels of each vertex reached. The
 *                           src vertex is assigned level 0. If a vertex i is not
 *                           reached, parent(i) is not present.
 *                           The level vector is not computed if NULL.
 * @param[out]    parent     If non-NULL on input, on successful return, it
 *                           contains parent vertex IDs for each vertex reached.
 *                           The src vertex will have itself as its parent. If a
 *                           vertex i is not reached, parent(i) is not present.
 *                           The parent vector is not computed if NULL.
 * @param[in]     G          The graph, directed or undirected.
 * @param[in]     src        The index of the src vertex (0-based)
 * @param[in]     pushpull   if true, use push/pull; otherwise, use pushonly.
 *                           Push/pull requires G->AT, G->rowdegree,
 *                           and library support.
 *                           TODO: consider removing this option or reverse logic
 * @param[out]    msg        Error message if a failure code is returned.
 *
 * @TODO pick return values that do not conflict with GraphBLAS errors.
 *
 * @retval GrB_SUCCESS      successful
 * @retval LAGRAPH_INVALID_GRAPH Graph is invalid (LAGraph_CheckGraph failed)
 */

LAGRAPH_PUBLIC
int LAGraph_BreadthFirstSearch
(
    GrB_Vector    *level,
    GrB_Vector    *parent,
    LAGraph_Graph  G,
    GrB_Index      src,
    bool           pushpull,
    char          *msg
);

//****************************************************************************
// TODO the following is a draft:
typedef enum
{
    LAGRAPH_CENTRALITY_BETWEENNESS = 0,     // node or edge centrality
    LAGRAPH_CENTRALITY_PAGERANKGAP = 1,     // GAP-style PageRank
    LAGRAPH_CENTRALITY_PAGERANK = 2,        // PageRank (handle dangling nodes)
    // ...
}
LAGraph_Centrality_Kind ;

LAGRAPH_PUBLIC
int LAGraph_VertexCentrality    // TODO: write this
(
    // output:
    GrB_Vector *centrality,     // centrality(i): centrality metric of node i
    // inputs:
    LAGraph_Graph G,            // input graph
    LAGraph_Centrality_Kind kind,    // kind of centrality to compute
//  int accuracy,               // TODO?: 0:quick, 1:better, ... max:exact
    char *msg
) ;

/****************************************************************************
 *
 * Count the triangles in a graph.
 *
 * @param[out]    ntriangles On successful return, contains the number of tris.
 * @param[in,out] G          The graph, symmetric, no self loops.
 * @param[out]    msg        Error message if a failure code is returned.
 *
 * @TODO pick return values that do not conflict with GraphBLAS errors.
 *
 * @retval GrB_SUCCESS          successful
 * @retval LAGRAPH_INVALID_GRAPH Graph is invalid (LAGraph_CheckGraph failed)
 * @retval GrB_NULL_POINTER     ntriangles is NULL
 * @retval FIXME:RETVAL         G->ndiag (self loops) is nonzero
 * @retval FIXME:RETVAL         graph is not symmetric
 * @retval FIXME:RETVAL         G->rowdegree was not precalculated (for modes 3-6)
 */
LAGRAPH_PUBLIC
int LAGraph_TriangleCount   // returns 0 if successful, < 0 if failure
(
    uint64_t      *ntriangles,   // # of triangles
    // input:
    LAGraph_Graph  G,
    char          *msg
) ;

// TODO: this is a Basic method, since G is input/output.
LAGRAPH_PUBLIC
int LAGraph_ConnectedComponents
(
    // output
    GrB_Vector *component,  // component(i)=s if node i is in the component
                            // whose representative node is s
    // inputs
    LAGraph_Graph G,        // input graph
    char *msg
) ;

// TODO: add AIsAllPositive or related as a G->property.
// TODO: Should a Basic method pick delta automatically?
LAGRAPH_PUBLIC
int LAGraph_SingleSourceShortestPath
(
    // output:
    GrB_Vector *path_length,    // path_length (i) is the length of the shortest
                                // path from the source vertex to vertex i
    // inputs:
    LAGraph_Graph G,
    GrB_Index source,           // source vertex
    int32_t delta,              // delta value for delta stepping
                                // TODO: use GxB_Scalar for delta?
    // TODO: make this an enum, and add to LAGraph_Graph properties, and then
    // remove it from the inputs to this function
    //      case 0: A can have negative, zero, or positive entries
    //      case 1: A can have zero or positive entries
    //      case 2: A only has positive entries
    // TODO: add AIsAllPositive to G->A_is_something...
    bool AIsAllPositive,       // A boolean indicating whether the entries of
                               // matrix A are all positive
    char *msg
) ;

//==============================================================================
// LAGraph Advanced algorithms
//==============================================================================

// The Advanced methods require the caller to select the algorithm and choose
// any parameter settings.  G is not modified, and so it is an input-only
// parameter to these methods.  If an Advanced algorithm requires a graph
// property to be computed, it must be computed prior to calling the Advanced
// method.

LAGRAPH_PUBLIC
int LAGraph_VertexCentrality_Betweenness
(
    // output:
    GrB_Vector *centrality,     // centrality(i): betweeness centrality of i
    // inputs:
    LAGraph_Graph G,            // input graph
    const GrB_Index *sources,   // source vertices to compute shortest paths
    int32_t ns,                 // number of source vertices
    char *msg
) ;

LAGRAPH_PUBLIC
int LAGraph_VertexCentrality_PageRankGAP
(
    // outputs:
    GrB_Vector *centrality, // centrality(i): GAP-style pagerank of node i
    // inputs:
    LAGraph_Graph G,        // input graph
    float damping,          // damping factor (typically 0.85)
    float tol,              // stopping tolerance (typically 1e-4) ;
    int itermax,            // maximum number of iterations (typically 100)
    int *iters,             // output: number of iterations taken
    char *msg
) ;

/****************************************************************************
 *
 * Count the triangles in a graph. Advanced API
 *
 * @param[out]    ntriangles On successful return, contains the number of tris.
 * @param[in]     G          The graph, symmetric, no self loops, and for some methods
 *                           (3-6), must have the row degree property calculated
 * @param[in]     method     specifies which algorithm to use
 *                             0:  use the default method
 *                             1:  Burkhardt:  ntri = sum (sum ((A^2) .* A)) / 6
 *                             2:  Cohen:      ntri = sum (sum ((L * U) .* A)) / 2
 *                             3:  Sandia:     ntri = sum (sum ((L * L) .* L))
 *                             4:  Sandia2:    ntri = sum (sum ((U * U) .* U))
 *                             5:  SandiaDot:  ntri = sum (sum ((L * U') .* L)).
 *                             6:  SandiaDot2: ntri = sum (sum ((U * L') .* U)).
 * @param[in,out] presort    controls the presort of the graph. If set
 *                           to 2 on input, presort will be set to sort type used
 *                           on output:
 *                             0: no sort
 *                             1: sort by degree, ascending order
 *                            -1: sort by degree, descending order
 *                             2: auto selection:
 *                                Method   Presort return value
 *                                1        0
 *                                2        0
 *                                3        1
 *                                4       -1
 *                                5        1
 *                                6       -1
 * @param[out]    msg        Error message if a failure code is returned.
 *
 * @retval GrB_SUCCESS          successful
 * @retval FIXME:RETVAL         invalid method value
 * @retval LAGRAPH_INVALID_GRAPH Graph is invalid (LAGraph_CheckGraph failed)
 * @retval GrB_NULL_POINTER     ntriangles is NULL
 * @retval FIXME:RETVAL         G->ndiag (self loops) is nonzero
 * @retval FIXME:RETVAL         graph is not "known" to be symmetric
 * @retval FIXME:RETVAL         G->rowdegree was not precalculated (for modes 3-6)
 */

typedef enum
{
    LAGraph_TriangleCount_Default = 0,      // use default method
    LAGraph_TriangleCount_Burkhardt = 1,    // sum (sum ((A^2) .* A)) / 6
    LAGraph_TriangleCount_Cohen = 2,        // sum (sum ((L * U) .* A)) / 2
    LAGraph_TriangleCount_Sandia = 3,       // sum (sum ((L * L) .* L))
    LAGraph_TriangleCount_Sandia2 = 4,      // sum (sum ((U * U) .* U))
    LAGraph_TriangleCount_SandiaDot = 5,    // sum (sum ((L * U') .* L))
    LAGraph_TriangleCount_SandiaDot2 = 6,   // sum (sum ((U * L') .* U))
}
LAGraph_TriangleCount_Method ;

typedef enum
{
    LAGraph_TriangleCount_NoSort = 0,       // no sort
    LAGraph_TriangleCount_Ascending = 1,    // sort by degree, ascending order
    LAGraph_TriangleCount_Descending = -1,  // sort by degree, descending order
    LAGraph_TriangleCount_AutoSort = 2,     // auto selection: no sort if rule
    // is not triggered.  Otherwise: sort in ascending order for methods 3 and
    // 5, descending ordering for methods 4 and 6.  On output, presort is
    // modified to reflect the sorting method used (0, -1, or 1).  If presort
    // is NULL on input, no sort is performed.
}
LAGraph_TriangleCount_Presort ;

LAGRAPH_PUBLIC
int LAGraph_TriangleCount_Methods
(
    uint64_t       *ntriangles,
    LAGraph_Graph   G,
    LAGraph_TriangleCount_Method    method,
    LAGraph_TriangleCount_Presort *presort,
    char           *msg
) ;

#endif
