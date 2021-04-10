//------------------------------------------------------------------------------
// LAGraph.h: user-visible include file for LAGraph (TODO: rename LAGraph.h)
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
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
// include files
//==============================================================================

#include <time.h>
#include <ctype.h>

#include <GraphBLAS.h>
#include <LAGraph_platform.h>
#include <LAGraph_internal.h>

//#include <complex.h>
// "I" is defined by <complex.h>, but is used in LAGraph and GraphBLAS to
// denote a list of row indices; remove it here.
//#undef I

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
        LAGraph_Free ((void **) &W, W_size) ;                       \
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
        return (0) ;
    }

#endif

#define LAGraph_TRY(LAGraph_method)             \
{                                               \
    int LAGraph_status = LAGraph_method ;       \
    if (LAGraph_status < 0)                     \
    {                                           \
        LAGraph_CATCH (LAGraph_status) ;        \
    }                                           \
}

//------------------------------------------------------------------------------
// GrB_TRY: try a GraphBLAS method and check for errors
//------------------------------------------------------------------------------

// LAGraph provides a similar functionality for calling GraphBLAS methods.
// GraphBLAS returns info = 0 (GrB_SUCCESS) or 1 (GrB_NO_VALUE) on success, and
// a value > 1 on failure.  The user application must #define GrB_CATCH to use
// GrB_TRY.  Note that GraphBLAS_info is internal to this macro.  If the
// user application or LAGraph method wants a copy, a statement such as
// info = GraphBLAS_info ; where info is defined outside of this macro.

#define GrB_TRY(GrB_method)                                                  \
{                                                                            \
    GrB_Info GraphBLAS_info = GrB_method ;                                   \
    if (! (GraphBLAS_info == GrB_SUCCESS || GraphBLAS_info == GrB_NO_VALUE)) \
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

extern void * (* LAGraph_Malloc_function  ) (size_t)         ;
extern void * (* LAGraph_Calloc_function  ) (size_t, size_t) ;
extern void * (* LAGraph_Realloc_function ) (void *, size_t) ;
extern void   (* LAGraph_Free_function    ) (void *)         ;
extern bool LAGraph_Malloc_is_thread_safe ;

// LAGraph_Malloc:  allocate a block of memory (wrapper for malloc)
void *LAGraph_Malloc
(
    size_t nitems,          // number of items
    size_t size_of_item,    // size of each item
    // output:
    size_t *size_allocated  // # of bytes actually allocated
) ;

// LAGraph_Calloc:  allocate a block of memory (wrapper for calloc)
void *LAGraph_Calloc
(
    size_t nitems,          // number of items
    size_t size_of_item,    // size of each item
    // output:
    size_t *size_allocated  // # of bytes actually allocated
) ;

void *LAGraph_Realloc       // returns pointer to reallocated block of memory,
(                           // or original block if reallocation fails.
    size_t nitems_new,      // new number of items in the object
    size_t nitems_old,      // old number of items in the object
    size_t size_of_item,    // size of each item
    // input/output
    void *p,                // old block to reallocate
    size_t *size_allocated, // # of bytes actually allocated
    // output
    bool *ok                // true if successful, false otherwise
) ;

// LAGraph_Free:  free a block of memory (wrapper for free)
void LAGraph_Free
(
    void **p,               // pointer to object to free, does nothing if NULL
    size_t size_allocated   // # of bytes actually allocated
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
// directed graphs.  Both kinds of graphs can be weighted or unweighted (in the
// future).  Additional types of graphs will be added in the future.

typedef enum
{
    // A(i,j) is the edge (i,j)
    // undirected: A is square, symmetric (both tril and triu present)
    // directed: A is square, unsymmetric or might happen to symmetric
    LAGRAPH_ADJACENCY_UNDIRECTED = 0,
    LAGRAPH_ADJACENCY_DIRECTED = 1,

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

// The LAGraph_Graph object contains a GrB_Matrix A as its primary component,
// as the adjacency matrix of the graph.  Typically, A(i,j) denotes the edge
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
//      A_pattern_is_symmetric: true if the pattern of A is symmetric

struct LAGraph_Graph_struct
{

    //--------------------------------------------------------------------------
    // primary components of the graph
    //--------------------------------------------------------------------------

    GrB_Matrix   A;          // the adjacency matrix of the graph
    GrB_Type     A_type;     // the type of scalar stored in A
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
    // changes, these properties can either be recomputed or deleted to denoted
    // the fact that they are unknown.  This choice is up to individual LAGraph
    // methods and utilities (TODO: define a default rule, and give user
    // control over this decision).

    // LAGraph methods can set non-scalar properties only if they are
    // constructing the graph.  They cannot modify them or create them if the
    // graph is declared as a read-only object in the parameter list of the
    // method.

    // TODO: discuss this: "However, scalar properties can be set even if the
    // graph is a read-only parameter, but only if they are accessed with
    // OpenMP atomic read/write operations."  OK?

    GrB_Matrix AT;          // AT = A', the transpose of A
    GrB_Type   AT_type;     // The type of scalar stored in AT

    GrB_Vector rowdegree;   // a GrB_INT64 vector of length m, if A is m-by-n.
           // where rowdegree(i) is the number of entries in A(i,:).  If
           // rowdegree is sparse and the entry rowdegree(i) is not present,
           // then it is assumed to be zero.
    GrB_Type   rowdegree_type;   // the type of scalar stored in rowdegree

    GrB_Vector coldegree ;  // a GrB_INT64 vector of length n, if A is m-by-n.
            // where coldegree(j) is the number of entries in A(:,j).  If
            // coldegree is sparse and the entry coldegree(j) is not present,
            // then it is assumed to be zero.  If A is known to have a
            // symmetric pattern, the convention is that the degree is held in
            // rowdegree, and coldegree is left as NULL.
    GrB_Type   coldegree_type;   // the type of scalar stored in coldegree

    LAGraph_BooleanProperty A_pattern_is_symmetric ;    // For an undirected
            // graph, this property will always be implicitly true and can be
            // ignored.  The matrix A for a directed weighted graph will
            // typically by unsymmetric, but might have a symmetric pattern.
            // In that case, this scalar property can be set to true.

    int64_t ndiag ; // # of entries on the diagonal of A, or -1 if unknown.
            // For the adjacency matrix of a directed or undirected graph,
            // this is the # of self-edges in the graph.
            // TODO: discuss this.

    // possible future cached properties:
    // GrB_Vector rowsum, colsum ;
    // rowsum (i) = sum (A (i,:)), regardless of kind
    // colsum (j) = sum (A (:,j)), regardless of kind
    // LAGraph_BooleanProperty connected ;   // true if G is a connected graph

    //--------------------------------------------------------------------------
    // for internal use only -- do not modify this parameter
    //--------------------------------------------------------------------------

    size_t size ;
} ;

typedef struct LAGraph_Graph_struct *LAGraph_Graph ;

//==============================================================================
// LAGraph utilities
//==============================================================================

// LAGraph_Init: start GraphBLAS and LAGraph
int LAGraph_Init (char *msg) ;      // return 0 if success, -1 if failure

// LAGraph_Xinit: start GraphBLAS and LAGraph, and set malloc/etc functions
int LAGraph_Xinit           // returns 0 if successful, -1 if failure
(
    // pointers to memory management functions
    void * (* user_malloc_function  ) (size_t),
    void * (* user_calloc_function  ) (size_t, size_t),
    void * (* user_realloc_function ) (void *, size_t),
    void   (* user_free_function    ) (void *),
    bool user_malloc_is_thread_safe,
    char *msg
) ;

// LAGraph_Finalize: finish LAGraph
int LAGraph_Finalize (char *msg) ;  // returns 0 if successful, -1 if failure

// LAGraph_MIN/MAX: suitable for integers, and non-NaN floating point
#define LAGraph_MIN(x,y) (((x) < (y)) ? (x) : (y))
#define LAGraph_MAX(x,y) (((x) > (y)) ? (x) : (y))

// LAGraph_New: create a new graph
int LAGraph_New         // returns 0 if successful, -1 if failure
(
    LAGraph_Graph *G,       // the graph to create, NULL if failure
    GrB_Matrix    *A,       // the adjacency matrix of the graph, may be NULL
    GrB_Type       A_type,  // the type of scalar stored in A
    LAGraph_Kind   kind,    // the kind of graph, may be LAGRAPH_UNKNOWN
    char *msg
) ;

// LAGraph_Delete: free a graph and all its contents
int LAGraph_Delete      // returns 0 if successful, -1 if failure
(
    LAGraph_Graph *G,   // the graph to delete; G set to NULL on output
    char *msg
) ;

// LAGraph_DeleteProperties: free any internal cached properties of a graph
int LAGraph_DeleteProperties    // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // G stays valid, only cached properties are freed
    char *msg
) ;

// LAGraph_Property_AT: construct G->AT for a graph
int LAGraph_Property_AT     // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to compute G->AT for
    char *msg
) ;

// LAGraph_Property_ASymmetricPattern: determine G->A_pattern_is_symmetric
int LAGraph_Property_ASymmetricPattern  // 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to determine the symmetry of pattern of A
    char *msg
) ;

// LAGraph_Property_RowDegree: determine G->rowdegree
int LAGraph_Property_RowDegree  // 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to determine G->rowdegree
    char *msg
) ;

// LAGraph_Property_ColDegree: determine G->coldegree
int LAGraph_Property_ColDegree  // 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to determine G->coldegree
    char *msg
) ;

// LAGraph_CheckGraph: determine if a graph is valid
int LAGraph_CheckGraph      // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to check
    // TODO: int level,     // 0:quick, O(1), 1:a bit more, 2: still more, 3: exhaustive!
    char *msg
) ;

// LAGraph_GetNumThreads: determine # of OpenMP threads to use
int LAGraph_GetNumThreads   // returns 0 if successful, or -1 if failure
(
    int *nthreads,          // # of threads to use
    char *msg
) ;

// LAGraph_SetNumThreads: set # of OpenMP threads to use
int LAGraph_SetNumThreads   // returns 0 if successful, or -1 if failure
(
    int nthreads,           // # of threads to use
    char *msg
) ;

// LAGraph_Tic: start the timer
int LAGraph_Tic             // returns 0 if successful, -1 if failure
(
    double tic [2],         // tic [0]: seconds, tic [1]: nanoseconds
    char *msg
) ;

// LAGraph_Toc: return time since last call to LAGraph_Tic
int LAGraph_Toc             // returns 0 if successful, -1 if failure
(
    double *t,              // time since last call to LAGraph_Tic
    const double tic [2],   // tic from last call to LAGraph_Tic
    char *msg
) ;

// ascii header prepended to all *.grb files
#define LAGRAPH_BIN_HEADER 512

// LAGraph_BinRead: read a matrix from a binary file
int LAGraph_BinRead         // returns 0 if successful, -1 if failure
(
    GrB_Matrix *A,          // matrix to read from the file
    GrB_Type   *A_type,     // type of the scalar stored in A
    FILE *f,                // file to read it from, already open
    char *msg
) ;

// LAGraph_MMRead: read a matrix from a Matrix Market file
int LAGraph_MMRead          // returns 0 if successful, -1 if faillure
(
    GrB_Matrix *A,          // handle of matrix to create
    GrB_Type   *A_type,     // type of the scalar stored in A
    FILE *f,                // file to read from, already open
    char *msg
) ;

// TODO: compression

// LAGraph_Pattern: return the pattern of a matrix (spones(A) in MATLAB)
int LAGraph_Pattern     // return 0 if successful, -1 if failure
(
    GrB_Matrix *C,      // a boolean matrix with the pattern of A
    GrB_Matrix A,
    char *msg
) ;

// LAGraph_IsEqual: compare two matrices for exact equality
int LAGraph_IsEqual         // returns 0 if successful, -1 if failure
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp op,        // for A and B with arbitrary user-defined types.
                            // Ignored if A and B have built-in types.
    char *msg
) ;

// LAGraph_IsAll: compare two matrices (types can differ)
int LAGraph_IsAll           // returns 0 if successful, -1 if failure
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp op,        // GrB_EQ_<type>, for the type of A and B,
                            // to check for equality.  Or use any desired
                            // operator.  The operator should return GrB_BOOL.
    char *msg
) ;

// LAGraph_TypeName: return the name of a type
int LAGraph_TypeName        // returns 0 if successful, -1 if failure
(
    char **name,            // name of the type
    GrB_Type type,          // GraphBLAS type
    char *msg
) ;

// LAGraph_KindName: return the name of a kind
int LAGraph_KindName        // returns 0 if successful, -1 if failure
(
    char **name,            // name of the kind
    LAGraph_Kind kind,      // graph kind
    char *msg
) ;

// LAGraph_DisplayGraph: print the contents of a graph
int LAGraph_DisplayGraph    // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to display
    int pr,                 // 0: nothing, 1: terse, 2: summary, 3: all,
                            // 4: same as 2 but with %0.15g for doubles
                            // 5: same as 3 but with %0.15g for doubles
    char *msg
) ;

// LAGraph_SortByDegree: sort a graph by its row or column degree
int LAGraph_SortByDegree    // returns 0 if successful, -1 if failure
(
    // output
    int64_t **P_handle,     // P is returned as a permutation vector of size n
    size_t *P_size,         // size of P in bytes
    // input
    LAGraph_Graph G,        // graph of n nodes
    bool byrow,             // if true, sort G->rowdegree, else G->coldegree
    bool ascending,         // sort in ascending or descending order
    char *msg
) ;

// LAGraph_SampleDegree: sample the degree median and mean
int LAGraph_SampleDegree        // returns 0 if successful, -1 if failure
(
    double *sample_mean,        // sampled mean degree
    double *sample_median,      // sampled median degree
    // input
    LAGraph_Graph G,        // graph of n nodes
    bool byrow,             // if true, sample G->rowdegree, else G->coldegree
    int64_t nsamples,       // number of samples
    uint64_t seed,          // random number seed
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

// Basic method: computes and *saves* any properties.
// G is input/output.

int LAGraph_BreadthFirstSearch      // returns -1 on failure, 0 if successful
(
    // outputs:
    GrB_Vector *level,      // not computed if NULL
    GrB_Vector *parent,     // not computed if NULL
    // inputs:
    LAGraph_Graph G,        // graph to traverse
    GrB_Index src,          // source node
    bool pushpull,          // if true, use push/pull, otherwise pushonly
    char *msg
) ;

// the following is a draft:
typedef enum
{
    LAGRAPH_CENTRALITY_BETWEENNESS = 0,     // node or edge centrality
    LAGRAPH_CENTRALITY_PAGERANKGAP = 1,     // GAP-style PageRank
    LAGRAPH_CENTRALITY_PAGERANK = 2,        // PageRank (handle dangling nodes)
    // ...
}
LAGraph_Centrality_Kind ;

int LAGraph_VertexCentrality    // returns -1 on failure, 0 if successful
(
    // output:
    GrB_Vector *centrality,     // centrality(i): centrality metric of node i
    // inputs:
    LAGraph_Graph G,            // input graph
    LAGraph_Centrality_Kind kind,    // kind of centrality to compute
//  int accuracy,               // TODO?: 0:quick, 1:better, ... max:exact
    char *msg
) ;

// TODO: pick a default method (say 5, with presort = 2 (auto)).
// Compute G->ndiag, and G->rowdegree if needed.  Determine if G->A is
// symmetric, if not known.
int LAGraph_TriangleCount   // returns 0 if successful, < 0 if failure
(
    uint64_t *ntriangles,   // # of triangles
    // input:
    LAGraph_Graph G,
    char *msg
) ;

// TODO: this is a Basic method, since G is input/output.
int LAGraph_ConnectedComponents
(
    // output
    GrB_Vector *component,  // component(i)=k if node i is in the kth component
    // inputs
    LAGraph_Graph G,        // input graph, G->A can change (ptr, not contents)
    char *msg
) ;

// TODO: add AIsAllPositive or related as a G->property.
// TODO: Should a Basic method pick delta automatically?
int LAGraph_SingleSourceShortestPath    // returns 0 if successful, -1 if fail
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
    //      case 2: A only has positive entries (see FIXME below)
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

int LAGraph_VertexCentrality_Betweenness    // vertex betweenness-centrality
(
    // output:
    GrB_Vector *centrality,     // centrality(i): betweeness centrality of i
    // inputs:
    LAGraph_Graph G,            // input graph
    const GrB_Index *sources,   // source vertices to compute shortest paths
    int32_t ns,                 // number of source vertices
    char *msg
) ;

int LAGraph_VertexCentrality_PageRankGAP // returns -1 on failure, 0 on success
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

int LAGraph_TriangleCount_Methods   // returns 0 if successful, < 0 if failure
(
    uint64_t *ntriangles,   // # of triangles
    // input:
    LAGraph_Graph G,
    int method,             // selects the method to use (TODO: enum)
    // input/output:
    int *presort,           // controls the presort of the graph (TODO: enum)
        //  0: no sort
        //  1: sort by degree, ascending order
        // -1: sort by degree, descending order
        //  2: auto selection: no sort if rule is not triggered.  Otherise:
        //  sort in ascending order for methods 3 and 5, descending ordering
        //  for methods 4 and 6.  On output, presort is modified to reflect the
        //  sorting method used (0, -1, or 1).  If presort is NULL on input, no
        //  sort is performed.
    char *msg
) ;

int LAGraph_TriangleCount_vanilla   // returns 0 if successful, < 0 if failure
(
    uint64_t *ntriangles,   // # of triangles
    // input:
    LAGraph_Graph G,
    int method,             // selects the method to use (TODO: enum)
    // input/output:
    int *presort,           // controls the presort of the graph (TODO: enum)
        //  0: no sort
        //  1: sort by degree, ascending order
        // -1: sort by degree, descending order
        //  2: auto selection: no sort if rule is not triggered.  Otherise:
        //  sort in ascending order for methods 3 and 5, descending ordering
        //  for methods 4 and 6.  On output, presort is modified to reflect the
        //  sorting method used (0, -1, or 1).  If presort is NULL on input, no
        //  sort is performed.
    char *msg
) ;

#endif
