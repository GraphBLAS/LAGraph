

//==============================================================================
// working notes only
//==============================================================================

#include "LAGraph.h"
#define TRY(method) LAGraph_TRY(method)


// TODO: discuss sanity checking.
// what if A is unsymmetric, but the graph is labeled "undirected"?
// error?  use A+A' ?

// LAGraph_Algorithm_Variant ?
GrB_Info LAGraph_BreadthFirstSearch_variant

// rule: NULL vs GrB_NULL:  use NULL (discuss this later)

// rule: if any output (*arg) is NULL, do not compute

// rule: if any input arg is NULL: use default, don't use it, etc

// rule:  output vectors and matrices are "new'd" by the LAGraph function:
    // usage:
    GrB_Vector level ;
    LAGraph_BreadthFirstSearch (&level, NULL, ...)

//------------------------------------------------------------------------------
// ALGO: breadth first search
//------------------------------------------------------------------------------

// Policy: in general, if an input *x is NULL, do not compute

// TODO: hops for all methods?

GrB_Vector level ;
GrB_Vector parent ;
int LAGraph_BreadthFirstSearch     // no _Variant suffix
(
    // outputs:
    GrB_Vector *level,              // if NULL do not compute
    GrB_Vector *parent,             // if NULL do not compute
    // inputs:
    LAGraph_Graph G,
    GrB_Index source,
    char *message
) ;

GrB_Vector parent ; // undefined
result = LAGraph_BreadthFirstSearch (NULL, &parent, G, 0, NULL) ;
if result is bad, parent = NULL

int LAGraph_BreadthFirstSearch_MultiSource
(
    // outputs:
    GrB_Matrix *level,              // if NULL do not compute
    GrB_Matrix *parent,             // if NULL do not compute
    // inputs:
    LAGraph_Graph G,
    GrB_Index *sources, size_t nsources     // or LAGraph_array?  GrB_array?
    char *message
) ;

info = LAGraph_whatever (&level, NULL, G, Src, 5) ;
if info is 'error' then level is returned as level == NULL.

// expert function:
GrB_Info LAGraph_BreadthFirstSearch_Frontier
(
    // input/output:
    GrB_Vector level,
    GrB_Vector parent,
    GrB_Vector frontier,
    // inputs:
    LAGraph_Graph G,
    GrB_Vector vertex_mask,         // filter the vertices
    bool vertex_mask_complement,
    bool vertex_mask_structural,
//  or:
//  GrB_Descriptor desc,    // mask complement, mask structural:
//                          // NULL, GrB_DESC_C, GrB_DESC_S, or GrB_DESC_SC
    uint64_t hops           // if hops=0, no limit (or INT64_MAX), or "NONE"
) ;

//------------------------------------------------------------------------------
// ALGO: connected components:
//------------------------------------------------------------------------------

// TODO: utility to do SCC = Comm * G * Comm'
// TODO: utility to do Gnew = Perm * G * Perm' ; or Gnew = G (P,P)
// TODO: utility to convert Perm matrix <==> perm "vector" (array?) GrB_Index*

// SCC graph = Community * G * Community', a contracted graph

    uint64_t nc ;
    GrB_Vector c ;
    LAGraph_Community_Map (&nc, &c, Commun, NULL) ;     // constructor, like dup
    ...
    GrB_free (&c) ;


GrB_Info LAGraph_Community_Map
(
    // output:
    uint64_t *ncomponents,
    GrB_Vector *component,  // vertex i is in component(i), a dense vector.
    // input:
    GrB_Matrix Community,   // Community(c,i)=1 becomes c=component(i)
    char *msg
) ;

// Gnew = Perm * G *Perm', then Gnew is block diagonal
GrB_Info LAGraph_Community_Permutation
(
    // output:
    GrB_Matrix *Perm,        // Perm(...), n-by-n
    // input:
    GrB_Matrix Community    // Community(c,i)=1 becomes c=component(i)
) ;

GrB_Info LAGraph_ConnectedComponents
(
    // output:
    GrB_Matrix *Community,  // k-by-n if there are k components.  GrB_BOOL.
                            // Community(c,i)=1 if c=component(i)
                            // TODO: or n-by-k?
    // input:
    LAGraph_Graph G         // if directed: strongly cc
) ;

GrB_Info LAGraph_ConnectedComponents_Weakly
(
    // output:
    GrB_Matrix *Community,  // k-by-n if there are k components.  GrB_BOOL.
                            // Community(c,i)=1 if c=component(i)
    // input:
    LAGraph_Graph G
) ;

//------------------------------------------------------------------------------
// ALGO: centrality:
//------------------------------------------------------------------------------

// TODO: utility to remove self loops, or add all self-loops (G=G+I), ...
// TODO: utility like trace(G), trace-like methods
// TODO: G =G+G', G = spones(G)

// random number generator?  vertex sampler?
// exact vs approx

// TODO make 2 different typedef enums, for VertexC. and EdgeC.

typedef enum
{
    LAGR_BETWEENNESS_FP32 = 4,
    LAGR_BETWEENNESS_FP64 = 5, 
    ...
}
LAGraph_VertexCentrality_Type ;

// this: easy mode
GrB_Vector centrality ; // uninit
LAGraph_something (&centrality, ...)
// then for 'expert' mode
GrB_Vector x ; // 
LAGraph_something_init (&x, ...)     or GrB_Vector_new (&x, n, 1, INT32)
    for a billion times:
        LAGraph_something_tiny (x, ...)        // reuse it, might be faster

GrB_Info LAGraph_VertexCentrality
(
    // output:
    GrB_Vector *centrality,

    // input:
    LAGraph_Graph G,
    LAGraph_VertexCentrality_Type kind
                            // betweenness, eigenvector, degree, pagerank, ...
                            // which pagerank? parameters to pagerank?
                            // LAGR_BETWEENNESS_FP32,
                            // LAGR_BETWEENNESS_FP64
) ;

LAGraph_VertexCentrality (&degree, G, LAGR_OUTDEGREE) ;

// utility function:
LAGraph_Degree (&degree, G, kind)
// or:
LAGraph_Degree_in (&degree, G)
LAGraph_Degree_out (&degree, G) // ... etc
LAGraph_Degree_etc (&degree, G) // ... etc

GrB_Info LAGraph_EdgeCentrality
(
    // output:
    GrB_Matrix *Centrality,

    // input:
    LAGraph_Graph G,
    LAGraph_EdgeCentrality_Type kind
                            // betweeness, eigenvector, degree, pagerank, ...
                            // LAGR_BETWEENNESS_FP32,
                            // LAGR_BETWEENNESS_FP64
) ;


//------------------------------------------------------------------------------
// ALGO: shortest paths:
//------------------------------------------------------------------------------

GrB_Info LAGraph_ShortestPath_[...]

GrB_Info LAGraph_ShortestPath_SingleSource
(
    // output: if NULL do not compute
    GrB_Vector *distance,       // type: INT64, FP32, or FP64
                                // INT64: if G is int*
                                // UINT64: if G is bool, uint*
                                // FP32: if G is FP32
                                // FP64: if G is FP64
    GrB_Vector *parent,         // tree
    GrB_Vector *hops,           // # of hops, level from source (call it "level")?

    // input:
    GrB_Index source,
    LAGraph_Graph G,
    // error handling (every LAGraph function has this last):
    char *message
) ;

    // negative-weight-cycle: result is not defined, must stop!
    // report info < 0 (an error)




GrB_Info LAGraph_ShortestPath_AllPairs
(
    // output:
    GrB_Matrix *Distance,
    GrB_Matrix *Parent,
    GrB_Matrix *Hops,
    // input:
    LAGraph_Graph G,
    // input/output
    char *message
) ;

// Dist, Parent undefined on input, created on output, just like GrB_Matrix_new
info = LAGraph_something (&Dist, &Parent, NULL, G) ;

//------------------------------------------------------------------------------
// what's next?
//------------------------------------------------------------------------------

GAP: BFS, SSSP, TriangleCount, Conn.Components, PageRank, BetweennessCentrality

LDBC: Local Clustering Coef, CDLP, different PageRank

Luby's MIS
minimum spanning forest
k-truss
max flow
Louvain community detection

Notation

Doxygen

//------------------------------------------------------------------------------
// ALGO: triangle counting
//------------------------------------------------------------------------------

// LAGraph_Graph G: is a read-only object if "input"

int LAGraph_TriangleCount
(
    uint64_t *ntriangles,   // # of triangles
    // input:
    LAGraph_Graph G,
    // input/output:
    char *message
) ;

int LAGraph_TriangleCount_expert
(
    uint64_t *ntriangles,   // # of triangles
    // input:
    int method,
    LAGraph_Graph G,
    // input/output:
    char *message
) ;

// TODO:  # triangles incident on each node/edge

//------------------------------------------------------------------------------
// ALGO: k-truss
//------------------------------------------------------------------------------

    // this should OK: G = func(G); we know G is aliased in/out
    LAGraph_anything (/* out: */ &G, k, /* in: */ G, &message) ;

    G = NULL ;      // make a new G
    LAGraph_anything (/* out: */ &G, k, /* in: */ H, &message) ;


    LAGraph_Graph G ;
    LAGraph_anything (/* out: */ &G, k, /* in: */ G, &message) ;

LAGraph_Graph Gnew = LAGraph_create ( ) ;
LAGraph_anything (/* out: */ Gnew, k, /* in: */ G, &message) ;

int LAGraph_Ktruss_next
(
    // input/output:
    LAGraph_Graph Gnext,
    // input:
    uint64_t k,
    LAGraph_Graph Gprev,        // (k-1)-truss on input, with support A(i,j)
    // input/output:
    char *message
) ;

LAGraph_Graph K = NULL ;

int LAGraph_Ktruss_all
(
    // output:
    LAGraph_Graph *K,           // an array of LAGraph_Graph, each with support
    // output:
    uint64_t *kmax,                
    LAGraph_Graph G,            // just a graph
    // input/output:
    char *message
) ;

    // k = 12 on output:
    // K[3], ... K[11].

