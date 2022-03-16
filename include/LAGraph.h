//------------------------------------------------------------------------------
// LAGraph.h: user-visible include file for LAGraph
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
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

// See also the LAGraph_Version utility method, which returns these values.
// These definitions must match the same definitions in LAGraph/CMakeLists.txt.
// FIXME: use config to create include/LAGraph.h from LAGraph/CMakeLists.txt
#define LAGRAPH_DATE "Mar 16, 2022"
#define LAGRAPH_VERSION_MAJOR 0
#define LAGRAPH_VERSION_MINOR 9
#define LAGRAPH_VERSION_UPDATE 14

//==============================================================================
// include files and helper macros
//==============================================================================

#include <GraphBLAS.h>
#if defined ( _OPENMP )
    #include <omp.h>
#endif

// LAGRAPH_MIN/MAX: suitable for integers, and non-NaN floating point
#define LAGRAPH_MIN(x,y) (((x) < (y)) ? (x) : (y))
#define LAGRAPH_MAX(x,y) (((x) > (y)) ? (x) : (y))

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
#if !defined ( LAGRAPH_VANILLA )
    // by default, set LAGRAPH_VANILLA to false
    #define LAGRAPH_VANILLA 0
#endif

#if ( !LAGRAPH_VANILLA ) && defined ( GxB_SUITESPARSE_GRAPHBLAS )
    // use SuiteSparse, and its GxB* extensions
    #define LAGRAPH_SUITESPARSE 1
#else
    // use any GraphBLAS library (possibly SuiteSparse) but with no GxB*
    #define LAGRAPH_SUITESPARSE 0
#endif

// maximum length of the name of a GrB type, including the null-terminator
#if LAGRAPH_SUITESPARSE
#define LAGRAPH_MAX_NAME_LEN GxB_MAX_NAME_LEN
#else
#define LAGRAPH_MAX_NAME_LEN 128
#endif

//==============================================================================
// LAGraph error handling
//==============================================================================

// All LAGraph functions return an int to indicate an error status (described
// below), where zero (GrB_SUCCESS) denotes success, negative values indicate
// an error, and positive values denote success but with some kind of
// algorithm-specific note or warning.

// In addition, all LAGraph functions also have a final parameter that is a
// pointer to a user-allocated string in which an algorithm-specific error
// message can be returned.  If NULL, no error message is returned.  This is
// not itself an error condition, it just indicates that the caller does not
// need the message returned.  If the message string is provided but no error
// occurs, an empty string is returned.

// LAGRAPH_MSG_LEN: The maximum required length of a message string
#define LAGRAPH_MSG_LEN 256

// For example, the following call computes the breadth-first-search of an
// LAGraph_Graph G, starting at a given source node.  It returns a status of
// zero if it succeeds and a negative value on failure.

/*
    GrB_Vector level, parent ;
    char msg [LAGRAPH_MSG_LEN] ;
    int status = LAGr_BreadthFirstSearch (&level, &parent, G, src, msg) ;
    if (status < 0)
    {
        printf ("error: %s\n", msg) ;
        // take corrective action ...
    }
*/

//------------------------------------------------------------------------------
// LAGraph error status:
//------------------------------------------------------------------------------

// All LAGraph methods return an int to denote their status:  zero if they are
// successful (which is the value of GrB_SUCCESS), negative on error, or
// positive for an informational value (such as GrB_NO_VALUE).  Integers in the
// range -999 to 999 are reserved for GraphBLAS GrB_Info return values:

//  successful results:
//  GrB_SUCCESS = 0             // all is well
//  GrB_NO_VALUE = 1            // A(i,j) requested but not there

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

// LAGraph returns any errors it receives from GraphBLAS, and also uses the
//  GrB_* error codes in these cases:
//  GrB_INVALID_INDEX: if a source node id is out of range
//  GrB_INVALID_VALUE: if an enum to select a method or option is out of range
//  GrB_NOT_IMPLEMENTED: if a type is not supported, or when SuiteSparse
//      GraphBLAS is required.

// Many LAGraph methods share common error cases, described below.  These
// return values are in the range -1000 to -1999.  Return values of -2000 or
// greater may be used by specific LAGraph methods, which denote errors not in
// this list:

//  LAGRAPH_INVALID_GRAPH:  the input graph is invalid; the details are given
//      in the error msg string returned by the method.

//  LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED: the method requires an undirected
//      graph, or a directed graph with an adjacency matrix that is known to
//      have a symmetric structure.  LAGraph_Property_SymmetricStructure can
//      be used to determine this property.

//  LAGRAPH_IO_ERROR:  a file input or output method failed, or an input file
//      has an incorrect format that cannot be parsed.

//  LAGRAPH_PROPERTY_MISSING:  some methods require one or more cached
//      properties to be computed before calling them (AT, rowdegree, or
//      coldegree.  Use LAGraph_Property_AT, LAGraph_Property_RowDegree,
//      and/or LAGraph_Property_ColDegree to compute them.

//  LAGRAPH_NO_SELF_EDGES_ALLOWED:  some methods requires that the graph have
//      no self edges, which correspond to the entries on the diagonal of the
//      adjacency matrix.  If the G->ndiag property is nonzero or unknown, this
//      error condition is returned.  Use LAGraph_Property_NDiag to compute
//      G->ndiag, or LAGraph_DeleteDiag to remove all diagonal entries from
//      G->A.

//  LAGRAPH_CONVERGENCE_FAILURE:  an iterative process failed to converge to
//      a good solution.

#define LAGRAPH_INVALID_GRAPH                   (-1000)
#define LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED    (-1001)
#define LAGRAPH_IO_ERROR                        (-1002)
#define LAGRAPH_PROPERTY_MISSING                (-1003)
#define LAGRAPH_NO_SELF_EDGES_ALLOWED           (-1004)
#define LAGRAPH_CONVERGENCE_FAILURE             (-1005)

// @retval GrB_SUCCESS if successful
// @retval a negative GrB_Info value on error (in range -999 to -1)
// @retval a positive GrB_Info value if successful but with extra information
//         (in range 1 to 999)

// @retval -1999 to -1000: a common LAGraph-specific error, given in the list
//          above
// @retval 1000 to 1999: if successful, with extra LAGraph-specific information

// @retval <= -2000 an LAGraph error specific to a particular LAGraph method
//         >= 2000

// @param[in,out] msg   any error messages

//------------------------------------------------------------------------------
// LAGRAPH_TRY: try an LAGraph method and check for errors
//------------------------------------------------------------------------------

// In a robust application, the return values from each call to LAGraph and
// GraphBLAS should be checked, and corrective action should be taken if an
// error occurs.  The LAGRAPH_TRY and GRB_TRY macros assist in this effort.

// LAGraph and GraphBLAS are written in C, and so they cannot rely on the
// try/catch mechanism of C++.  To accomplish a similar goal, each LAGraph file
// must #define its own file-specific macro called LAGRAPH_CATCH.  The typical
// usage of macro is to free any temporary matrices/vectors or workspace when
// an error occurs, and then "throw" the error by returning to the caller.  A
// user application may also #define LAGRAPH_CATCH and use these macros.

// A typical example of a user function that calls LAGraph might #define
// LAGRAPH_CATCH as follows.  Suppose workvector is a GrB_vector used for
// computations internal to the mybfs function, and W is a (double *) space
// created by malloc.

#if example_usage_only

    // an example user-defined LAGRAPH_CATCH macro
    #define LAGRAPH_CATCH(status)                                   \
    {                                                               \
        /* an LAGraph error has occurred */                         \
        printf ("LAGraph error: (%d): file: %s, line: %d\n%s\n",    \
            status, __FILE__, __LINE__, msg) ;                      \
        /* free any internal workspace and return the status */     \
        GrB_free (*parent) ;                                        \
        GrB_free (workvector) ;                                     \
        LAGraph_Free ((void **) &W, NULL) ;                         \
        return (status) ;                                           \
    }

    // an example user function that uses LAGRAPH_TRY / LAGRAPH_CATCH
    int mybfs (LAGraph_Graph G, GrB_Vector *parent, int64_t src)
    {
        GrB_Vector workvector = NULL ;
        double *W = NULL ;
        char msg [LAGRAPH_MSG_LEN] ;
        (*parent) = NULL ;
        LAGRAPH_TRY (LAGr_BreadthFirstSearch (NULL, parent, G, src, true,
            msg)) ;
        // ...
        return (GrB_SUCCESS) ;
    }

#endif

#define LAGRAPH_TRY(LAGraph_method)             \
{                                               \
    int LG_status = LAGraph_method ;            \
    if (LG_status < GrB_SUCCESS)                \
    {                                           \
        LAGRAPH_CATCH (LG_status) ;             \
    }                                           \
}

//------------------------------------------------------------------------------
// GRB_TRY: try a GraphBLAS method and check for errors
//------------------------------------------------------------------------------

// LAGraph provides a similar functionality for calling GraphBLAS methods.
// GraphBLAS returns info = 0 (GrB_SUCCESS) or 1 (GrB_NO_VALUE) on success, and
// a value < 0 on failure.  The user application must #define GRB_CATCH to use
// GRB_TRY.  Note that GraphBLAS_info is internal to this macro.  If the
// user application or LAGraph method wants a copy, a statement such as
// info = GraphBLAS_info ; where info is defined outside of this macro.

// GraphBLAS and LAGraph both use the convention that negative values are
// errors, and the LAGraph_status is a superset of the GrB_Info enum.  As a
// result, the user can define LAGRAPH_CATCH and GRB_TRY as the same operation.
// The main difference between the two would be the error message string.  For
// LAGraph, the string is the last parameter, and LAGRAPH_CATCH can optionally
// print it out.  For GraphBLAS, the GrB_error mechanism can return a string.

#define GRB_TRY(GrB_method)                     \
{                                               \
    GrB_Info LG_GrB_Info = GrB_method ;         \
    if (LG_GrB_Info < GrB_SUCCESS)              \
    {                                           \
        GRB_CATCH (LG_GrB_Info) ;               \
    }                                           \
}

//==============================================================================
// LAGraph memory management
//==============================================================================

// LAGraph provides wrappers for the malloc/calloc/realloc/free set of memory
// management functions, initialized by LAGraph_Init or LAGr_Init.  By default,
// the following are pointers to the ANSI C11 malloc/calloc/realloc/free
// functions.

LAGRAPH_PUBLIC void * (* LAGraph_Malloc_function  ) (size_t)         ;
LAGRAPH_PUBLIC void * (* LAGraph_Calloc_function  ) (size_t, size_t) ;
LAGRAPH_PUBLIC void * (* LAGraph_Realloc_function ) (void *, size_t) ;
LAGRAPH_PUBLIC void   (* LAGraph_Free_function    ) (void *)         ;

//------------------------------------------------------------------------------
// LAGraph_Malloc:  allocate a block of memory (wrapper for malloc)
//------------------------------------------------------------------------------

LAGRAPH_PUBLIC
int LAGraph_Malloc
(
    // output:
    void **p,               // pointer to allocated block of memory
    // input:
    size_t nitems,          // number of items
    size_t size_of_item,    // size of each item
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Calloc:  allocate a block of memory (wrapper for calloc)
//------------------------------------------------------------------------------

LAGRAPH_PUBLIC
int LAGraph_Calloc
(
    // output:
    void **p,               // pointer to allocated block of memory
    // input:
    size_t nitems,          // number of items
    size_t size_of_item,    // size of each item
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Realloc: reallocate a block of memory (wrapper for realloc)
//------------------------------------------------------------------------------

LAGRAPH_PUBLIC
int LAGraph_Realloc
(
    // input/output:
    void **p,               // old block to reallocate
    // input:
    size_t nitems_new,      // new number of items in the object
    size_t nitems_old,      // old number of items in the object
    size_t size_of_item,    // size of each item
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Free:  free a block of memory (wrapper for free)
//------------------------------------------------------------------------------

// If LAGraph_Malloc (&p, ...) is the pointer to the allocated block of memory,
// LAGraph_Free (&p, ...) is the method to free it.  The parameter is passed as
// &p so that p can be set to NULL on return, to guard against double-free.
// LAGraph_Free does nothing if &p or p are NULL on input.

LAGRAPH_PUBLIC
int LAGraph_Free            // free a block of memory and set p to NULL
(
    // input/output:
    void **p,               // pointer to object to free, does nothing if NULL
    char *msg
) ;

//==============================================================================
// LAGraph data structures
//==============================================================================

// In addition to relying on underlying GraphBLAS objects (GrB_Matrix,
// GrB_Vector, GrB_Descriptor, ...), LAGraph adds the LAGraph_Graph.  This
// object contains a representation of a graph and its associated data.  Unlike
// the GrB_* objects, the LAGraph_Graph is not opaque.

// LAGRAPH_UNKNOWN is used for all scalars whose value is not known
#define LAGRAPH_UNKNOWN (-1)

// FIXME: start here on Mar 16, 2022

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
    LAGraph_ADJACENCY_UNDIRECTED = 0, // A is square and symmetric; both upper
                                      // and lower triangular parts are
                                      // present.  A(i,j) is the edge (i,j)

    LAGraph_ADJACENCY_DIRECTED = 1,   // A is square; A(i,j) is the edge (i,j)

    // possible future kinds of graphs:
    // LAGraph_ADJACENCY_UNDIRECTED_UNWEIGHTED
    // LAGraph_ADJACENCY_DIRECTED_UNWEIGHTED
    // LAGraph_ADJACENCY_UNDIRECTED_TRIL
    // LAGraph_ADJACENCY_UNDIRECTED_TRIU
    // LAGraph_BIPARTITE
    // LAGraph_BIPARTITE_DIRECTED
    // LAGraph_BIPARTITE_UNDIRECTED
    // LAGraph_INCIDENCE_*
    // LAGraph_MULTIGRAPH_*
    // LAGraph_HYPERGRAPH
    // LAGraph_HYPERGRAPH_DIRECTED
    // ...

    LAGraph_KIND_UNKNOWN = LAGRAPH_UNKNOWN      // the graph kind is unknown
}
LAGraph_Kind ;

//------------------------------------------------------------------------------
// LAGraph_BooleanProperty: true, false, or unknown
//------------------------------------------------------------------------------

typedef enum
{
    LAGraph_FALSE = 0,
    LAGraph_TRUE = 1,
    LAGraph_BOOLEAN_UNKNOWN = LAGRAPH_UNKNOWN
}
LAGraph_BooleanProperty ;

//------------------------------------------------------------------------------
// LAGraph_BoundKind: exact, bound, approximate, or unknown
//------------------------------------------------------------------------------

// LAGraph_BoundKind describes the status of a graph property or other metric.
// If the metric is computed in floating-point arithmetic, it may have been
// computed with roundoff error, but it may still be declared as "exact" if the
// roundoff error is expected to be small, or if the metric was computed as
// carefully as possible (to within reasonable roundoff error).  The
// "approximate" state is used when the metric is a rough estimate, not because
// of roundoff error but because of other algorithmic approximations.  The
// decision of when to tag a metric as "exact" or "approximate" is up to the
// particular algorithm, which each algorithm must document.

// The "bound" state indicates that the metric is an upper or lower bound,
// depending on the particular metric.  If computed in floating-point
// arithmetic, an "upper bound" metric may be actually slightly lower than the
// actual upper bound, because of floating-point roundoff.

typedef enum
{
    LAGraph_EXACT = 0,      // the metric is exact (possibly ignoring roundoff)
    LAGraph_BOUND = 1,      // the metric is a bound (upper or lower, depending
                            // on the particular metric)
    LAGraph_APPROX = 2,     // the metric is a rough approximation
    LAGraph_BOUND_UNKNOWN = LAGRAPH_UNKNOWN,
}
LAGraph_BoundKind ;

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
//      structure_is_symmetric: true if the structure of A is symmetric
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

    GrB_Matrix AT ;         // AT = A', the transpose of A, with the same type.

    GrB_Vector rowdegree ;  // a GrB_INT64 vector of length m, if A is m-by-n.
           // where rowdegree(i) is the number of entries in A(i,:).  If
           // rowdegree is sparse and the entry rowdegree(i) is not present,
           // then it is assumed to be zero.

    GrB_Vector coldegree ;  // a GrB_INT64 vector of length n, if A is m-by-n.
            // where coldegree(j) is the number of entries in A(:,j).  If
            // coldegree is sparse and the entry coldegree(j) is not present,
            // then it is assumed to be zero.  If A is known to have a
            // symmetric structure, the convention is that the degree is held in
            // rowdegree, and coldegree is left as NULL.

    // If G is held as an incidence matrix, then G->A might be rectangular,
    // in the future, but the graph G may have a symmetric structure anyway.
    LAGraph_BooleanProperty structure_is_symmetric ;    // For an undirected
            // graph, this property will always be implicitly true and can be
            // ignored.  The matrix A for a directed weighted graph will
            // typically be unsymmetric, but might have a symmetric structure.
            // In that case, this scalar property can be set to true.
            // By default, this property is set to LAGRAPH_UNKNOWN.

    int64_t ndiag ; // # of entries on the diagonal of A, or LAGRAPH_UNKNOWN if
            // unknown.  For the adjacency matrix of a directed or undirected
            // graph, this is the number of self-edges in the graph.

    GrB_Scalar emin ;   // minimum edge weight: exact, lower bound, or estimate
    LAGraph_BoundKind emin_kind ;
            // EXACT: emin is exactly equal to the smallest entry, min(G->A)
            // BOUND: emin <= min(G->A)
            // APPROX: emin is a rough estimate of min(G->A)
            // UNKNOWN: emin is unknown

    GrB_Scalar emax ;   // maximum edge weight: exact, upper bound, or estimate
    LAGraph_BoundKind emax_kind ;
            // EXACT: emax is exactly equal to the largest entry, max(G->A)
            // BOUND: emax >= max(G->A)
            // APPROX: emax is a rough estimate of max(G->A)
            // UNKNOWN: emax is unknown

    // possible future cached properties:

    // Some algorithms may want to know if the graph has any edge weights
    // exactly equal to zero.  In some cases, this can be inferred from the
    // emin/emax bounds, or it can be indicated via the following property:
    // LAGraph_BooleanProperty nonzero ;  // If true, then all entries in
            // G->A are known to be nonzero.  If false, G->A may contain
            // entries in its structure that are identically equal to zero.  If
            // unknown, then G->A may or may not have entries equal to zero.
    // other edge weight metrics: median, standard deviation....  Might be
    // useful for computing Delta for a Basic SSSP.
    // GrB_Vector rowsum, colsum ;
    // rowsum(i) = sum(A(i,:)), regardless of kind
    // colsum(j) = sum(A(:,j)), regardless of kind
    // LAGraph_BooleanProperty connected ;   // true if G is a connected graph
} ;

typedef struct LAGraph_Graph_struct *LAGraph_Graph ;

//==============================================================================
// LAGraph utilities
//==============================================================================

//------------------------------------------------------------------------------
// LAGraph_Init: start GraphBLAS and LAGraph
//------------------------------------------------------------------------------

// This method must be called before calling any other GrB* or LAGraph* method.
// It initializes GraphBLAS with GrB_init and then performs LAGraph-specific
// initializations.  In particular, the LAGraph semirings listed below are
// created.

LAGRAPH_PUBLIC
int LAGraph_Init (char *msg) ;

// LAGraph semirings, created by LAGraph_Init or LAGr_Init:
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

//------------------------------------------------------------------------------
// LAGraph_Version: determine the version of LAGraph
//------------------------------------------------------------------------------

// The version number and date can also be obtained via compile-time constants
// from LAGraph.h.  However, it is possible to compile a user application that
// #includes one version of LAGraph.h and then links with another version of
// the LAGraph library later on, so the version number and date may differ from
// the compile-time constants.

// The LAGraph_Version method allows the library itself to be
// queried, after it is linked in with the user application.

// The version_number array is set to LAGRAPH_VERSION_MAJOR,
// LAGRAPH_VERSION_MINOR, and LAGRAPH_VERSION_UPDATE, in that order.
// The LAGRAPH_DATE string is copied into the user-provided version_date
// string, and is null-terminated.

LAGRAPH_PUBLIC
int LAGraph_Version
(
    // output:
    int version_number [3],     // user-provided array of size 3
    char version_date [LAGRAPH_MSG_LEN],    // user-provided array
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Finalize: finish LAGraph
//------------------------------------------------------------------------------

// LAGraph_Finalize must be called as the last LAGraph method.  It calls
// GrB_finalize and frees any LAGraph objects created by LAGraph_Init or
// LAGr_Init.  After calling this method, no LAGraph or GraphBLAS methods
// may be used.

LAGRAPH_PUBLIC
int LAGraph_Finalize (char *msg) ;

//------------------------------------------------------------------------------
// LAGraph_New: create a new graph
//------------------------------------------------------------------------------

// LAGraph_New creates a new graph G.  The properties G->AT, G->rowdegree, and
// G->coldegree are set to NULL, and scalar properties are set to
// LAGRAPH_UNKNOWN.

LAGRAPH_PUBLIC
int LAGraph_New
(
    // output:
    LAGraph_Graph *G,   // the graph to create, NULL if failure
    // input/output:
    GrB_Matrix    *A,   // the adjacency matrix of the graph, may be NULL.
                        // A is moved into G as G->A, and A itself is set
                        // to NULL to denote that is now a part of G.
                        // That is, { G->A = A ; A = NULL ; } is performed.
                        // When G is deleted, G->A is freed.  If A is NULL,
                        // the graph is invalid until G->A is set.
    // input:
    LAGraph_Kind kind,  // the kind of graph. This may be LAGRAPH_UNKNOWN,
                        // which must then be revised later before the
                        // graph can be used.
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Delete: free a graph and all its contents
//------------------------------------------------------------------------------

// LAGraph_Delete frees a graph G, including its adjacency matrix G->A and the
// cached properties G->AT, G->rowdegree, and G->coldegree.

LAGRAPH_PUBLIC
int LAGraph_Delete
(
    // input/output:
    LAGraph_Graph *G,   // the graph to delete; G set to NULL on output.
                        // All internal GrB_Matrix and GrB_Vector objects are
                        // freed, including G->A.  To keep G->A while deleting
                        // the graph G, use:
                        // { A = G->A ; G->A = NULL ; LAGraph_Delete (&G, msg);}
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_DeleteProperties: free any internal cached properties of a graph
//------------------------------------------------------------------------------

// LAGraph_DeleteProperties frees all cached properies of a graph G.  The graph
// is still valid.  This method should be used if G->A changes, since such
// changes will normally invalidate G->AT, G->rowdgree, and/or G->coldegree.

LAGRAPH_PUBLIC
int LAGraph_DeleteProperties
(
    // input/output:
    LAGraph_Graph G,    // G stays valid, only cached properties are freed
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Property_AT: construct G->AT for a graph
//------------------------------------------------------------------------------

// LAGraph_Property_AT constructs G->AT, the transpose of G->A.  This matrix is
// required by some of the algorithms.  Basic algorithms may construct G->AT if
// they require it.  The matrix G->AT is then available for subsequent use.
// If G->A changes, G->AT should be freed and recomputed.  If G->AT already
// exists, it is left unchanged.  As a result, if G->A changes, G->AT should
// be explictly freed.

LAGRAPH_PUBLIC
int LAGraph_Property_AT
(
    // input/output:
    LAGraph_Graph G,    // graph for which to compute G->AT
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Property_SymmetricStructure: determine G->structure_is_symmetric
//------------------------------------------------------------------------------

// LAGraph_Property_SymmetricStructure determines if the sparsity structure
// of G->A is symmetric (ignoring its values).  If G->kind denotes that the
// graph is undirected, this property is implicitly true (and not checked).
// Otherwise, this method determines if the structure of G->A for a directed
// graph G has a symmetric sparsity structure.  No work is performend if the
// property is already known.

LAGRAPH_PUBLIC
int LAGraph_Property_SymmetricStructure
(
    // input/output:
    LAGraph_Graph G,    // graph to determine the symmetry of structure of A
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Property_RowDegree: determine G->rowdegree
//------------------------------------------------------------------------------

// LAGraph_Property_RowDegree computes G->rowdegree.  No work is performed if
// it already exists in G.

LAGRAPH_PUBLIC
int LAGraph_Property_RowDegree
(
    // input/output:
    LAGraph_Graph G,    // graph to determine G->rowdegree
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Property_ColDegree: determine G->coldegree
//------------------------------------------------------------------------------

// LAGraph_Property_ColDegree computes G->coldegree.  No work is performed if
// it already exists in G.  If G is undirected, G->coldegree is never computed.
// Instead, G->rowdegree is used instead.  No work is performed it is already
// exists in G.

LAGRAPH_PUBLIC
int LAGraph_Property_ColDegree
(
    // input/output:
    LAGraph_Graph G,    // graph to determine G->coldegree
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Property_NDiag: determine G->ndiag
//------------------------------------------------------------------------------

// LAGraph_Property_NDiag computes G->ndiag, the number of diagonal entries
// that appear in the G->A matrix.  For an undirected or directed graph with an
// adjacency matrix G->A, these are the number of self-edges in G.  No work is
// performed it is already computed.

LAGRAPH_PUBLIC
int LAGraph_Property_NDiag
(
    // input/output:
    LAGraph_Graph G,    // graph to compute G->ndiag
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Property_EMin: determine G->emin
//------------------------------------------------------------------------------

// LAGraph_Property_EMin computes G->emin = min (G->A).

LAGRAPH_PUBLIC
int LAGraph_Property_EMin
(
    // input/output:
    LAGraph_Graph G,    // graph to determine G->emin
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Property_EMax: determine G->emax
//------------------------------------------------------------------------------

// LAGraph_Property_EMax computes G->emax = max (G->A).

LAGRAPH_PUBLIC
int LAGraph_Property_EMax
(
    // input/output:
    LAGraph_Graph G,    // graph to determine G->emax
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_DeleteDiag: remove all diagonal entries from G->A
//------------------------------------------------------------------------------

// LAGraph_DeleteDiag removes any diagonal entries from G->A.  Most properties
// are cleared or set to LAGRAPH_UNKNOWN.  G->ndiag is set to zero, and
// G->structure_is_symmetric is left unchanged.

LAGRAPH_PUBLIC
int LAGraph_DeleteDiag
(
    // input/output:
    LAGraph_Graph G,    // diagonal entries removed, most properties cleared
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_CheckGraph: determine if a graph is valid
//------------------------------------------------------------------------------

// LAGraph_CheckGraph determines if a graph is valid.  Only basic checks are
// performed on the cached properties, taking O(1) time.

LAGRAPH_PUBLIC
int LAGraph_CheckGraph
(
    // input/output:
    LAGraph_Graph G,    // graph to check
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_GetNumThreads: determine # of OpenMP threads to use
//------------------------------------------------------------------------------

// LAGraph_GetNumThreads determines the current number of OpenMP threads that
// can be used.  This is provided by SuiteSparse:GraphBLAS via a GxB extension,
// or by omp_get_max_threads() otherwise.  If OpenMP is not in use, then 1 is
// returned.

LAGRAPH_PUBLIC
int LAGraph_GetNumThreads
(
    // output:
    int *nthreads,      // # of threads to use
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_SetNumThreads: set # of OpenMP threads to use
//------------------------------------------------------------------------------

// LAGraph_SetNumThreads sets the current number of OpenMP threads that
// can be used.  This is provided by SuiteSparse:GraphBLAS via a GxB extension,
// or by omp_set_max_threads() otherwise.  If OpenMP is not in use, then this
// function silently returns with no error message.

LAGRAPH_PUBLIC
int LAGraph_SetNumThreads
(
    // input:
    int nthreads,       // # of threads to use
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Tic: start the timer
//------------------------------------------------------------------------------

// LAGraph_Tic starts a wallclock timer.  Normally, this is simply a wrapper
// for omp_get_wtime, if OpenMP is in use.  Otherwise, an OS-specific timing
// function is called.

LAGRAPH_PUBLIC
int LAGraph_Tic
(
    double tic [2],     // tic [0]: seconds, tic [1]: nanoseconds
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Toc: return time since last call to LAGraph_Tic
//------------------------------------------------------------------------------

// LAGraph_Toc returns the time since the last call to LAGraph_Tic with the
// same tic input array (which is not modified).

LAGRAPH_PUBLIC
int LAGraph_Toc
(
    // output:
    double *t,              // time since last call to LAGraph_Tic
    // input:
    const double tic [2],   // tic from last call to LAGraph_Tic
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_MMRead: read a matrix in MatrixMarket format
//------------------------------------------------------------------------------

/* The file format used here is compatible with all variations of the Matrix
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

LAGRAPH_PUBLIC
int LAGraph_MMRead
(
    // output:
    GrB_Matrix *A,  // handle of matrix to create
    // input:
    FILE *f,        // file to read from, already open
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_MMWrite: write a matrix in MatrixMarket format
//------------------------------------------------------------------------------

// LAGraph_MMWrite writes a matrix in MatrixMarket format.  Refer to
// LAGraph_MMRead for a description of the output file format.  The
// MatrixMarket header line always appears, followed by the second line
// containing the GraphBLAS type:
//      %%GraphBLAS type <entrytype>

LAGRAPH_PUBLIC
int LAGraph_MMWrite
(
    // input:
    GrB_Matrix A,       // matrix to write to the file
    FILE *f,            // file to write it to, must be already open
    FILE *fcomments,    // optional file with extra comments, may be NULL
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Matrix_Structure: return the structure of a matrix
//------------------------------------------------------------------------------

// LAGraph_Matrix_Structure returns the sparsity structure of a matrix A as a
// boolean (GrB_BOOL) matrix C.  If A(i,j) appears in the sparsity structure of
// A, then C(i,j) is set to true.  The sparsity structure of A and C are
// identical.

LAGRAPH_PUBLIC
int LAGraph_Matrix_Structure
(
    // output:
    GrB_Matrix *C,  // a boolean matrix with same structure of A, with C(i,j)
                    // set to true if A(i,j) appears in the sparsity structure
                    // of A.
    // input:
    GrB_Matrix A,
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Vector_Structure: return the structure of a vector
//------------------------------------------------------------------------------

// LAGraph_Vector_Structure return the sparsity structure of a vector u as a
// boolean (GrB_BOOL) vector w.  If u(i) appears in the sparsity structure of
// u, then w(i) is set to true.  The sparsity structure of u and w are
// identical.

LAGRAPH_PUBLIC
int LAGraph_Vector_Structure
(
    // output:
    GrB_Vector *w,  // a boolean vector with same structure of u, with w(i)
                    // set to true if u(i) appears in the sparsity structure
                    // of u.
    // input:
    GrB_Vector u,
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_NameOfType: return the name of a type
//------------------------------------------------------------------------------

// LAGraph_NameOfType returns the name of a GraphBLAS type as a string.  The
// names for the 11 built-in types (GrB_BOOL, GrB_INT8, etc) correspond to the
// names of the corresponding C types (bool, int8_t, etc).

LAGRAPH_PUBLIC
int LAGraph_NameOfType
(
    // output:
    char *name,     // name of the type: user provided array of size at
                    // least LAGRAPH_MAX_NAME_LEN.
    // input:
    GrB_Type type,  // GraphBLAS type
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_TypeFromName: return a GrB_Type from its name
//------------------------------------------------------------------------------

// LAGraph_TypeFromName returns the GrB_Type corresponding to its name.  That
// is, given the string "bool", this method returns GrB_BOOL.

LAGRAPH_PUBLIC
int LAGraph_TypeFromName
(
    // output:
    GrB_Type *type, // GraphBLAS type
    // input:
    char *name,     // name of the type: a null-terminated string
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_SizeOfType: return sizeof(...) of a GraphBLAS GrB_Type
//------------------------------------------------------------------------------

// LAGraph_SizeOfType returns sizeof(...) of a GraphBLAS GrB_Type.  For
// example, if given the GrB_Type of GrB_FP64, the value sizeof(double) is
// returned.

LAGRAPH_PUBLIC
int LAGraph_SizeOfType
(
    // output:
    size_t *size,   // size of the type
    // input:
    GrB_Type type,  // GraphBLAS type
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Matrix_TypeName: return the name of the GrB_Type of a GrB_Matrix
//------------------------------------------------------------------------------

// LAGraph_Matrix_TypeName returns the name of the GrB_Type of a GrB_Matrix.
// Currently, this method requires SuiteSparse:GraphBLAS.

LAGRAPH_PUBLIC
int LAGraph_Matrix_TypeName
(
    // output:
    char *name,     // name of the type of the matrix A (user-provided array
                    // of size at least LAGRAPH_MAX_NAME_LEN).
    // input:
    GrB_Matrix A,   // matrix to query
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Vector_TypeName: return the name of the GrB_Type of a GrB_Vector
//------------------------------------------------------------------------------

// LAGraph_Vector_TypeName returns the name of the GrB_Type of a GrB_Vector.
// Currently, this method requires SuiteSparse:GraphBLAS.

LAGRAPH_PUBLIC
int LAGraph_Vector_TypeName
(
    // output:
    char *name,     // name of the type of the vector v (user-provided array
                    // of size at least LAGRAPH_MAX_NAME_LEN).
    // input:
    GrB_Vector v,   // vector to query
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Scalar_TypeName: return the name of the GrB_Type of a GrB_Scalar
//------------------------------------------------------------------------------

// LAGraph_Scalar_TypeName returns the name of the GrB_Type of a GrB_Scalar.
// Currently, this method requires SuiteSparse:GraphBLAS.

LAGRAPH_PUBLIC
int LAGraph_Scalar_TypeName
(
    // output:
    char *name,     // name of the type of the scalar s (user-provided array
                    // of size at least LAGRAPH_MAX_NAME_LEN).
    // input:
    GrB_Scalar s,   // scalar to query
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_KindName: return the name of a kind
//------------------------------------------------------------------------------

// LAGraph_KindName: return the name of a graph kind.  For example, if given
// LAGaphH_ADJACENCY_UNDIRECTED, the string "undirected" is returned.

LAGRAPH_PUBLIC
int LAGraph_KindName
(
    // output:
    char *name,     // name of the kind (user provided array of size at least
                    // LAGRAPH_MAX_NAME_LEN)
    // input:
    LAGraph_Kind kind,  // graph kind
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_SortByDegree: sort a graph by its row or column degree
//------------------------------------------------------------------------------

// LAGraph_SortByDegree sorts the nodes of a graph by their row or column
// degrees.  The graph G->A itself is not changed.  Refer to LAGr_TriangleCount
// for an example of how to permute G->A after calling this function.

LAGRAPH_PUBLIC
int LAGraph_SortByDegree
(
    // output:
    int64_t **P_handle,     // P is returned as a permutation vector of size n
    // input:
    const LAGraph_Graph G,  // graph of n nodes
    bool byrow,             // if true, sort G->rowdegree, else G->coldegree
    bool ascending,         // sort in ascending or descending order
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_SampleDegree: sample the degree median and mean
//------------------------------------------------------------------------------

// LAGraph_SampleDegree computes an estimate of the median and mean of the row
// or column degree, by randomly sampling the G->rowdegree or G->coldegree
// vector.

LAGRAPH_PUBLIC
int LAGraph_SampleDegree
(
    // output:
    double *sample_mean,    // sampled mean degree
    double *sample_median,  // sampled median degree
    // input:
    const LAGraph_Graph G,  // graph of n nodes
    bool byrow,             // if true, sample G->rowdegree, else G->coldegree
    int64_t nsamples,       // number of samples
    uint64_t seed,          // random number seed
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_DisplayGraph: print the contents of a graph
//------------------------------------------------------------------------------

// LAGraph_DisplayGraph prints the contents of a graph to a file in a human-
// readable form.  This method is not meant for saving a graph to a file;
// see LAGraph_MMWrite for that method.

typedef enum
{
    LAGraph_SILENT = 0,     // nothing is printed
    LAGraph_SUMMARY = 1,    // print a terse summary
    LAGraph_SHORT = 2,      // short description, about 30 entries of a matrix
    LAGraph_COMPLETE = 3,   // print the entire contents of the object
    LAGraph_SHORT_VERBOSE = 4,    // LAGraph_SHORT but with "%.15g" for doubles
    LAGraph_COMPLETE_VERBOSE = 5  // LAGraph_COMPLETE, but "%.15g" for doubles
}
LAGraph_Print_Level ;

LAGRAPH_PUBLIC
int LAGraph_DisplayGraph
(
    // input:
    const LAGraph_Graph G,  // graph to display
    LAGraph_Print_Level pr, // print level (0 to 5)
    FILE *f,                // file to write to, must already be open
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Matrix_IsEqual: compare for exact equality
//------------------------------------------------------------------------------

// LAGraph_Matrix_IsEqual compares two matrices for exact equality.  If the two
// matrices must have different data types, the result is always false (no
// typecasting is performed).  Only the 11 built-in GrB* types are supported.

LAGRAPH_PUBLIC
int LAGraph_Matrix_IsEqual
(
    // output:
    bool *result,       // true if A == B, false if A != B or error
    // input:
    const GrB_Matrix A,
    const GrB_Matrix B,
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Matrix_IsEqual_op: check if two matrices are equal with given op
//------------------------------------------------------------------------------

// LAGraph_Matrix_IsEqual_op compares two matrices using the given binary
// operator.  The op may be built-in or user-defined.

LAGRAPH_PUBLIC
int LAGraph_Matrix_IsEqual_op
(
    // output:
    bool *result,           // true if A == B, false if A != B or error
    // input:
    const GrB_Matrix A,
    const GrB_Matrix B,
    const GrB_BinaryOp op,        // comparator to use
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Vector_IsEqual: check if two vectors are equal
//------------------------------------------------------------------------------

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
    // output:
    bool *result,           // true if A == B, false if A != B or error
    // input:
    const GrB_Vector A,
    const GrB_Vector B,
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Vector_IsEqual_op: check if two vectors are equal with given op
//------------------------------------------------------------------------------

/**
 * Checks if two vectors are identically equal (same size, pattern,
 * and values) according to a user specified comparator op.
 *
 *
 * @param[out]   result   Set to true on return is vectors are "equal"
 * @param[in]    A        First vector to compare
 * @param[in]    B        Second vector to compare
 * @param[in]    op       Binary operator to use for the comparisons
 * @param[out]   msg      If an error code is returned, this may hold an error msg.
 *
 * @retval GrB_SUCCESS          if completed successfully (equal or not)
 * @retval GrB_NULL_POINTER     result or op is NULL
 * @return Any GraphBLAS errors that may have been encountered
 */
LAGRAPH_PUBLIC
int LAGraph_Vector_IsEqual_op
(
    // output:
    bool *result,           // true if A == B, false if A != B or error
    // input:
    const GrB_Vector A,
    const GrB_Vector B,
    const GrB_BinaryOp op,        // comparator to use
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Matrix_Print: pretty-print a matrix
//------------------------------------------------------------------------------

// LAGraph_Matrix_Print displays a matrix in a human-readable form.  This
// method is not meant for saving a GrB_Matrix to a file; see LAGraph_MMWrite
// for that method.

LAGRAPH_PUBLIC
int LAGraph_Matrix_Print
(
    // input:
    const GrB_Matrix A,     // matrix to pretty-print to the file
    LAGraph_Print_Level pr, // print level (0 to 5)
    FILE *f,            // file to write it to, must be already open; use
                        // stdout or stderr to print to those locations.
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Vector_Print: pretty-print a matrix
//------------------------------------------------------------------------------

// LAGraph_Vector_Print displays a vector in a human-readable form.  This
// method is not meant for saving a GrB_Vector to a file.  To perform that
// operation, copy the GrB_Vector into an n-by-1 GrB_Matrix and use
// LAGraph_MMWrite.

LAGRAPH_PUBLIC
int LAGraph_Vector_Print
(
    // input:
    const GrB_Vector v,     // vector to pretty-print to the file
    LAGraph_Print_Level pr, // print level (0 to 5)
    FILE *f,            // file to write it to, must be already open; use
                        // stdout or stderr to print to those locations.
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Sort1: sort array of size n
//------------------------------------------------------------------------------

// LAGraph_Sort1 sorts an int64_t array of size n in ascending order.

LAGRAPH_PUBLIC
int LAGraph_Sort1
(
    // input/output:
    int64_t *A_0,       // size n array
    // input:
    const int64_t n,
    int nthreads,       // # of threads to use
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Sort2: sort two arrays of size n
//------------------------------------------------------------------------------

// LAGraph_Sort2 sorts two int64_t arrays A of size n in ascending order.
// The arrays are kept in the same order, where the pair (A_0 [k], A_1 [k]) is
// treated as a single pair.  The pairs are sorted by the first value A_0,
// with ties broken by A_1.

LAGRAPH_PUBLIC
int LAGraph_Sort2
(
    // input/output:
    int64_t *A_0,       // size n array
    int64_t *A_1,       // size n array
    // input:
    const int64_t n,
    int nthreads,       // # of threads to use
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGraph_Sort2: sort three arrays of size n
//------------------------------------------------------------------------------

// LAGraph_Sort3 sorts three int64_t arrays A of size n in ascending order.
// The arrays are kept in the same order, where the triplet (A_0 [k], A_1 [k],
// A_2 [k]) is treated as a single triplet.  The triplets are sorted by the
// first value A_0, with ties broken by A_1, and then by A_2 if the values of
// A_0 and A_1 are identical.

LAGRAPH_PUBLIC
int LAGraph_Sort3
(
    // input/output:
    int64_t *A_0,       // size n array
    int64_t *A_1,       // size n array
    int64_t *A_2,       // size n array
    // input:
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

// LAGraph Basic algorithms are named with the LAGraph_* prefix.

//------------------------------------------------------------------------------
// LAGraph_TriangleCount
//------------------------------------------------------------------------------

// This is a Basic algorithm (G->ndiag, G->rowdegree, G->structure_is_symmetric
// are computed, if not present).

/*
 * Count the triangles in a graph.
 *
 * @param[out]    ntriangles On successful return, contains the number of tris.
 * @param[in,out] G          The graph, symmetric, no self loops.
 * @param[out]    msg        Error message if a failure code is returned.
 */
LAGRAPH_PUBLIC
int LAGraph_TriangleCount
(
    // output:
    uint64_t      *ntriangles,   // # of triangles
    // input/output:
    LAGraph_Graph  G,
    char          *msg
) ;

//==============================================================================
// LAGraph Advanced algorithms and utilities
//==============================================================================

// The Advanced algorithms require the caller to select the algorithm and choose
// any parameter settings.  G is not modified, and so it is an input-only
// parameter to these methods.  If an Advanced algorithm requires a graph
// property to be computed, it must be computed prior to calling the Advanced
// method.

// Advanced algorithms are named with the LAGr_* prefix, to distinguish them
// from Basic algorithms.

//------------------------------------------------------------------------------
// LAGr_Init: start GraphBLAS and LAGraph, and set malloc/etc functions
//------------------------------------------------------------------------------

// LAGr_Init is identical to LAGraph_Init, except that it allows the user
// application to provide four memory management functions, replacing the
// standard malloc, calloc, realloc, and free.  The functions pointed to by
// user_malloc_function, user_calloc_function, user_realloc_function, and
// user_free_function have the same signature as the ANSI C malloc, calloc,
// realloc, and free functions, respectively.

// Only user_malloc_function and user_free_function are required.
// user_calloc_function may be NULL, in which case LAGraph_Calloc uses
// LAGraph_Malloc and memset.  Likewise, user_realloc_function may be NULL, in
// which case LAGraph_Realloc uses LAGraph_Malloc, memcpy, and LAGraph_Free.

LAGRAPH_PUBLIC
int LAGr_Init
(
    // input:
    void * (* user_malloc_function  ) (size_t),
    void * (* user_calloc_function  ) (size_t, size_t),
    void * (* user_realloc_function ) (void *, size_t),
    void   (* user_free_function    ) (void *),
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGr_BreadthFirstSearch: breadth-first search
//------------------------------------------------------------------------------

// This is an Advanced algorithm (G->AT and G->rowdgree are required).

/*
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
 * @param[out]    msg        Error message if a failure code is returned.
 *
 * @retval GrB_SUCCESS      successful
 * @retval LAGRAPH_INVALID_GRAPH Graph is invalid (LAGraph_CheckGraph failed)
 */

LAGRAPH_PUBLIC
int LAGr_BreadthFirstSearch
(
    // output:
    GrB_Vector    *level,
    GrB_Vector    *parent,
    // input:
    const LAGraph_Graph G,
    GrB_Index      src,
    char          *msg
) ;

//------------------------------------------------------------------------------
// LAGr_ConnectedComponents: connected components of an undirected graph
//------------------------------------------------------------------------------

// This is an Advanced algorithm (G->structure_is_symmetric must be known),

LAGRAPH_PUBLIC
int LAGr_ConnectedComponents
(
    // output:
    GrB_Vector *component,  // component(i)=s if node i is in the component
                            // whose representative node is s
    // input:
    const LAGraph_Graph G,  // input graph
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGr_SingleSourceShortestPath: single-source shortest path
//------------------------------------------------------------------------------

// This is an Advanced algorithm (G->emin is required).

// The graph G must have an adjacency matrix of type GrB_INT32, GrB_INT64,
// GrB_UINT32, GrB_UINT64, GrB_FP32, or GrB_FP64.  If G->A has any other type,
// GrB_NOT_IMPLEMENTED is returned.

// FUTURE: add a Basic algorithm that computes G->emin, G->emax, and then uses
// that information to compute an appropriate (estimated) Delta,

LAGRAPH_PUBLIC
int LAGr_SingleSourceShortestPath
(
    // output:
    GrB_Vector *path_length,    // path_length (i) is the length of the shortest
                                // path from the source vertex to vertex i
    // input:
    const LAGraph_Graph G,
    GrB_Index source,           // source vertex
    GrB_Scalar Delta,           // delta value for delta stepping
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGr_Betweenness: betweeness centrality metric
//------------------------------------------------------------------------------

// This is an Advanced algorithm (G->AT is required).

LAGRAPH_PUBLIC
int LAGr_Betweenness
(
    // output:
    GrB_Vector *centrality,     // centrality(i): betweeness centrality of i
    // input:
    const LAGraph_Graph G,      // input graph
    const GrB_Index *sources,   // source vertices to compute shortest paths
    int32_t ns,                 // number of source vertices
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGr_PageRank: pagerank
//------------------------------------------------------------------------------

// This is an Advanced algorithm (G->AT and G->rowdegree are required).

// LAGr_PageRank computes the standard pagerank of a
// directed graph G.  Sinks (nodes with no out-going edges) are handled.

LAGRAPH_PUBLIC
int LAGr_PageRank
(
    // output:
    GrB_Vector *centrality, // centrality(i): pagerank of node i
    int *iters,             // number of iterations taken
    // input:
    const LAGraph_Graph G,  // input graph
    float damping,          // damping factor (typically 0.85)
    float tol,              // stopping tolerance (typically 1e-4) ;
    int itermax,            // maximum number of iterations (typically 100)
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGr_PageRankGAP: GAP-style pagerank
//------------------------------------------------------------------------------

// This is an Advanced algorithm (G->AT and G->rowdegree are required).

// LAGr_PageRankGAP computes the GAP-style pagerank of a
// directed graph G.  Sinks (nodes with no out-going edges) are not handled.
// This method should be used for the GAP benchmark only, not for production.

LAGRAPH_PUBLIC
int LAGr_PageRankGAP
(
    // output:
    GrB_Vector *centrality, // centrality(i): GAP-style pagerank of node i
    int *iters,             // number of iterations taken
    // input:
    const LAGraph_Graph G,  // input graph
    float damping,          // damping factor (typically 0.85)
    float tol,              // stopping tolerance (typically 1e-4) ;
    int itermax,            // maximum number of iterations (typically 100)
    char *msg
) ;

//------------------------------------------------------------------------------
// LAGr_TriangleCount: triangle counting
//------------------------------------------------------------------------------

// This is an Advanced algorithm (G->ndiag, G->rowdegree,
// G->structure_is_symmetric are required).

/* Count the triangles in a graph. Advanced API
 *
 * @param[out]    ntriangles On successful return, contains the number of tris.
 * @param[in]     G          The graph, symmetric, no self loops, and for some methods
 *                           (3-6), must have the row degree property calculated
 * @param[in]     method     specifies which algorithm to use
 *                           0:  use the default method
 *                           1:  Burkhardt:  ntri = sum (sum ((A^2) .* A)) / 6
 *                           2:  Cohen:      ntri = sum (sum ((L * U) .* A)) / 2
 *                           3:  Sandia:     ntri = sum (sum ((L * L) .* L))
 *                           4:  Sandia2:    ntri = sum (sum ((U * U) .* U))
 *                           5:  SandiaDot:  ntri = sum (sum ((L * U') .* L)).
 *                           6:  SandiaDot2: ntri = sum (sum ((U * L') .* U)).
 * @param[in,out] presort    controls the presort of the graph. If set to 2 on
 *                           input, presort will be set to sort type used
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
 * @retval GrB_INVALID_VALUE    invalid method value
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
int LAGr_TriangleCount
(
    // output:
    uint64_t       *ntriangles,
    // input:
    const LAGraph_Graph G,
    LAGraph_TriangleCount_Method    method,
    LAGraph_TriangleCount_Presort *presort,
    char           *msg
) ;

#endif
