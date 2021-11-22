//------------------------------------------------------------------------------
// LAGraph_BF_full1.c: Bellman-Ford single-source shortest paths, returns tree,
// while diagonal of input matrix A needs not to be explicit 0
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// LAGraph_BF_full1: Bellman-Ford single source shortest paths, returning both
// the path lengths and the shortest-path tree.  contributed by Jinhao Chen and
// Tim Davis, Texas A&M.

// LAGraph_BF_full1 performs a Bellman-Ford to find out shortest path, parent
// nodes along the path and the hops (number of edges) in the path from given
// source vertex s in the range of [0, n) on graph given as matrix A with size
// n*n. The sparse matrix A has entry A(i, j) if there is an edge from vertex i
// to vertex j with weight w, then A(i, j) = w.

// LAGraph_BF_full1 returns GrB_SUCCESS if it succeeds.  In this case, there
// are no negative-weight cycles in the graph, and d, pi, and h are returned.
// The vector d has d(k) as the shortest distance from s to k. pi(k) = p+1,
// where p is the parent node of k-th node in the shortest path. In particular,
// pi(s) = 0. h(k) = hop(s, k), the number of edges from s to k in the shortest
// path.

// If the graph has a negative-weight cycle, GrB_NO_VALUE is returned, and the
// GrB_Vectors d(k), pi(k) and h(k)  (i.e., *pd_output, *ppi_output and
// *ph_output respectively) will be NULL when negative-weight cycle detected.

// Otherwise, other errors such as GrB_OUT_OF_MEMORY, GrB_INVALID_OBJECT, and
// so on, can be returned, if these errors are found by the underlying
// GrB_* functions.

//------------------------------------------------------------------------------

#define LAGraph_FREE_WORK              \
{                                      \
    GrB_free(&d);                      \
    GrB_free(&dmasked);                \
    GrB_free(&dless);                  \
    GrB_free(&Atmp);                   \
    GrB_free(&BF_Tuple3);              \
    GrB_free(&BF_lMIN_Tuple3);         \
    GrB_free(&BF_PLUSrhs_Tuple3);      \
    GrB_free(&BF_LT_Tuple3);           \
    GrB_free(&BF_lMIN_Tuple3_Monoid);  \
    GrB_free(&BF_lMIN_PLUSrhs_Tuple3); \
    LAGraph_Free ((void**)&I);                  \
    LAGraph_Free ((void**)&J);                  \
    LAGraph_Free ((void**)&w);                  \
    LAGraph_Free ((void**)&W);                  \
    LAGraph_Free ((void**)&h);                  \
    LAGraph_Free ((void**)&pi);                 \
}

#define LAGraph_FREE_ALL               \
{                                      \
    LAGraph_FREE_WORK                  \
    GrB_free (pd_output);              \
    GrB_free (ppi_output);             \
    GrB_free (ph_output);              \
}

#include <LAGraph.h>
#include <LAGraphX.h>
#include <LG_internal.h>  // from src/utility


typedef void (*LAGraph_binary_function) (void *, const void *, const void *) ;

//------------------------------------------------------------------------------
// data type for each entry of the adjacent matrix A and "distance" vector d;
// <INFINITY,INFINITY,INFINITY> corresponds to nonexistence of a path, and
// the value  <0, 0, NULL> corresponds to a path from a vertex to itself
//------------------------------------------------------------------------------
typedef struct
{
    double w;    // w  corresponds to a path weight.
    GrB_Index h; // h  corresponds to a path size or number of hops.
    GrB_Index pi;// pi corresponds to the penultimate vertex along a path.
                 // vertex indexed as 1, 2, 3, ... , V, and pi = 0 (as nil)
                 // for u=v, and pi = UINT64_MAX (as inf) for (u,v) not in E
}
BF1_Tuple3_struct;

//------------------------------------------------------------------------------
// 2 binary functions, z=f(x,y), where Tuple3xTuple3 -> Tuple3
//------------------------------------------------------------------------------
void BF1_lMIN
(
    BF1_Tuple3_struct *z,
    const BF1_Tuple3_struct *x,
    const BF1_Tuple3_struct *y
)
{
    if (x->w < y->w
        || (x->w == y->w && x->h < y->h)
        || (x->w == y->w && x->h == y->h && x->pi < y->pi))
    {
        if (z != x) { *z = *x; }
    }
    else
    {
        *z = *y;
    }
}

void BF1_PLUSrhs
(
    BF1_Tuple3_struct *z,
    const BF1_Tuple3_struct *x,
    const BF1_Tuple3_struct *y
)
{
    z->w = x->w + y->w;
    z->h = x->h + y->h;
    z->pi = (x->pi != UINT64_MAX && y->pi != 0) ?  y->pi : x->pi ;
}

void BF1_Identity
(
    BF1_Tuple3_struct *z,
    const BF1_Tuple3_struct *x
)
{
    *z = *x;
}

void BF1_LT
(
    bool *z,
    const BF1_Tuple3_struct *x,
    const BF1_Tuple3_struct *y
)
{
    (*z) = (x->w < y->w
        || (x->w == y->w && x->h < y->h)
        || (x->w == y->w && x->h == y->h && x->pi < y->pi)) ;
}

// Given a n-by-n adjacency matrix A and a source vertex s.
// If there is no negative-weight cycle reachable from s, return the distances
// of shortest paths from s and parents along the paths as vector d. Otherwise,
// returns d=NULL if there is a negtive-weight cycle.
// pd_output is pointer to a GrB_Vector, where the i-th entry is d(s,i), the
//   sum of edges length in the shortest path
// ppi_output is pointer to a GrB_Vector, where the i-th entry is pi(i), the
//   parent of i-th vertex in the shortest path
// ph_output is pointer to a GrB_Vector, where the i-th entry is h(s,i), the
//   number of edges from s to i in the shortest path
// A has weights on corresponding entries of edges
// s is given index for source vertex
GrB_Info LAGraph_BF_full1
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    GrB_Vector *ppi_output,     //the pointer to the vector of parent
    GrB_Vector *ph_output,      //the pointer to the vector of hops
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Index s           //given index of the source
)
{
    GrB_Info info;
    char *msg = NULL ;
    // tmp vector to store distance vector after n (i.e., V) loops
    GrB_Vector d = NULL, dmasked = NULL, dless = NULL;
    GrB_Matrix Atmp = NULL;
    GrB_Type BF_Tuple3;

    GrB_BinaryOp BF_lMIN_Tuple3;
    GrB_BinaryOp BF_PLUSrhs_Tuple3;
    GrB_UnaryOp BF_Identity_Tuple3;
    GrB_BinaryOp BF_LT_Tuple3;

    GrB_Monoid BF_lMIN_Tuple3_Monoid;
    GrB_Semiring BF_lMIN_PLUSrhs_Tuple3;

    GrB_Index nrows, ncols, n, nz;  // n = # of row/col, nz = # of nnz in graph
    GrB_Index *I = NULL, *J = NULL; // for col/row indices of entries from A
    GrB_Index *h = NULL, *pi = NULL;
    double *w = NULL;
    BF1_Tuple3_struct *W = NULL;

    LG_CHECK (A == NULL || pd_output == NULL ||
        ppi_output == NULL || ph_output == NULL, -1001, "inputs are NULL") ;

    *pd_output  = NULL;
    *ppi_output = NULL;
    *ph_output  = NULL;
    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_TRY (GrB_Matrix_nvals (&nz, A));
    LG_CHECK (nrows != ncols, -1002, "A must be square") ;
    n = nrows;

    LG_CHECK (s >= n || s < 0, -1003, "invalid source node") ;

    //--------------------------------------------------------------------------
    // create all GrB_Type GrB_BinaryOp GrB_Monoid and GrB_Semiring
    //--------------------------------------------------------------------------
    // GrB_Type
    GrB_TRY (GrB_Type_new(&BF_Tuple3, sizeof(BF1_Tuple3_struct)));

    // GrB_BinaryOp
    GrB_TRY (GrB_UnaryOp_new(&BF_Identity_Tuple3,
        (void*) (&BF1_Identity), BF_Tuple3, BF_Tuple3));
    GrB_TRY (GrB_BinaryOp_new(&BF_LT_Tuple3,
        (LAGraph_binary_function) (&BF1_LT), GrB_BOOL, BF_Tuple3, BF_Tuple3));
    GrB_TRY (GrB_BinaryOp_new(&BF_lMIN_Tuple3,
        (LAGraph_binary_function) (&BF1_lMIN), BF_Tuple3, BF_Tuple3, BF_Tuple3));
    GrB_TRY (GrB_BinaryOp_new(&BF_PLUSrhs_Tuple3,
        (LAGraph_binary_function)(&BF1_PLUSrhs),
        BF_Tuple3, BF_Tuple3, BF_Tuple3));

    // GrB_Monoid
    BF1_Tuple3_struct BF_identity = (BF1_Tuple3_struct) { .w = INFINITY,
        .h = UINT64_MAX, .pi = UINT64_MAX };
    LAGRAPH_OK(GrB_Monoid_new_UDT(&BF_lMIN_Tuple3_Monoid, BF_lMIN_Tuple3,
        &BF_identity));

    //GrB_Semiring
    GrB_TRY (GrB_Semiring_new(&BF_lMIN_PLUSrhs_Tuple3,
        BF_lMIN_Tuple3_Monoid, BF_PLUSrhs_Tuple3));

    //--------------------------------------------------------------------------
    // allocate arrays used for tuplets
    //--------------------------------------------------------------------------
    I = LAGraph_Malloc (nz, sizeof(GrB_Index)) ;
    J = LAGraph_Malloc (nz, sizeof(GrB_Index)) ;
    w = LAGraph_Malloc (nz, sizeof(double)) ;
    W = LAGraph_Malloc (nz, sizeof(BF1_Tuple3_struct)) ;
    LG_CHECK (I == NULL || J == NULL || w == NULL || W == NULL, -1004,
        "out of memory") ;

    //--------------------------------------------------------------------------
    // create matrix Atmp based on A, while its entries become BF_Tuple3 type
    //--------------------------------------------------------------------------
    LAGRAPH_OK(GrB_Matrix_extractTuples_FP64(I, J, w, &nz, A));
    int nthreads;
    LAGRAPH_OK(LAGraph_GetNumThreads (&nthreads, NULL)) ;
    printf ("nthreads %d\n", nthreads) ;
    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (GrB_Index k = 0; k < nz; k++)
    {
        W[k] = (BF1_Tuple3_struct) { .w = w[k], .h = 1, .pi = I[k] + 1 };
    }
    GrB_TRY (GrB_Matrix_new(&Atmp, BF_Tuple3, n, n));
    LAGRAPH_OK(GrB_Matrix_build_UDT(Atmp, I, J, W, nz, BF_lMIN_Tuple3));
    LAGraph_Free ((void**)&I);
    LAGraph_Free ((void**)&J);
    LAGraph_Free ((void**)&W);
    LAGraph_Free ((void**)&w);

    //--------------------------------------------------------------------------
    // create and initialize "distance" vector d, dmasked and dless
    //--------------------------------------------------------------------------
    GrB_TRY (GrB_Vector_new(&d, BF_Tuple3, n));
    // make d dense
    LAGRAPH_OK(GrB_Vector_assign_UDT(d, NULL, NULL, (void*)&BF_identity,
        GrB_ALL, n, NULL));
    // initial distance from s to itself
    BF1_Tuple3_struct d0 = (BF1_Tuple3_struct) { .w = 0, .h = 0, .pi = 0 };
    LAGRAPH_OK(GrB_Vector_setElement_UDT(d, &d0, s));

    // creat dmasked as a sparse vector with only one entry at s
    GrB_TRY (GrB_Vector_new(&dmasked, BF_Tuple3, n));
    LAGRAPH_OK(GrB_Vector_setElement_UDT(dmasked, &d0, s));

    // create dless
    GrB_TRY (GrB_Vector_new(&dless, GrB_BOOL, n));

    //--------------------------------------------------------------------------
    // start the Bellman Ford process
    //--------------------------------------------------------------------------
    bool any_dless= true;      // if there is any newly found shortest path
    int64_t iter = 0;          // number of iterations

    // terminate when no new path is found or more than V-1 loops
    while (any_dless && iter < n - 1)
    {
        // execute semiring on d and A, and save the result to dtmp
        GrB_TRY (GrB_vxm(dmasked, GrB_NULL, GrB_NULL,
            BF_lMIN_PLUSrhs_Tuple3, dmasked, Atmp, GrB_NULL));

        // dless = d .< dtmp
        //GrB_TRY (GrB_Vector_clear(dless));
        GrB_TRY (GrB_eWiseMult(dless, NULL, NULL, BF_LT_Tuple3, dmasked, d,
            NULL));

        // if there is no entry with smaller distance then all shortest paths
        // are found
        GrB_TRY (GrB_reduce (&any_dless, NULL, GrB_LOR_MONOID_BOOL, dless,
            NULL)) ;
        if(any_dless)
        {
            // update all entries with smaller distances
            GrB_TRY (GrB_apply(d, dless, NULL, BF_Identity_Tuple3, dmasked, NULL));

            // only use entries that were just updated
            GrB_TRY (GrB_Vector_clear(dmasked));
            GrB_TRY (GrB_apply(dmasked, dless, NULL, BF_Identity_Tuple3, d, NULL));
            //try:
            //GrB_TRY (GrB_assign(dmasked, dless, NULL, d, GrB_ALL, n, GrB_DESC_R);
        }
        iter ++;
    }

    // check for negative-weight cycle only when there was a new path in the
    // last loop, otherwise, there can't be a negative-weight cycle.
    if (any_dless)
    {
        // execute semiring again to check for negative-weight cycle
        GrB_TRY (GrB_vxm(dmasked, GrB_NULL, GrB_NULL,
            BF_lMIN_PLUSrhs_Tuple3, dmasked, Atmp, GrB_NULL));

        // dless = d .< dtmp
        //GrB_TRY (GrB_Vector_clear(dless));
        GrB_TRY (GrB_eWiseMult(dless, NULL, NULL, BF_LT_Tuple3, dmasked, d, NULL));

        // if there is no entry with smaller distance then all shortest paths
        // are found
        GrB_TRY (GrB_reduce (&any_dless, NULL, GrB_LOR_MONOID_BOOL, dless, NULL)) ;
        if(any_dless)
        {
            // printf("A negative-weight cycle found. \n");
            LAGraph_FREE_ALL;
            return (GrB_NO_VALUE) ;
        }
    }

    //--------------------------------------------------------------------------
    // extract tuple from "distance" vector d and create GrB_Vectors for output
    //--------------------------------------------------------------------------
    I = LAGraph_Malloc (n, sizeof(GrB_Index)) ;
    W = LAGraph_Malloc (n, sizeof(BF1_Tuple3_struct)) ;
    w = LAGraph_Malloc (n, sizeof(double)) ;
    h  = LAGraph_Malloc (n, sizeof(GrB_Index)) ;
    pi = LAGraph_Malloc (n, sizeof(GrB_Index)) ;
    LG_CHECK (I == NULL || W == NULL || w == NULL || h == NULL || pi == NULL,
        -1004, "out of memory") ;

    LAGRAPH_OK(GrB_Vector_extractTuples_UDT (I, (void *) W, &n, d));

    for (GrB_Index k = 0; k < n; k++)
    {
        w [k] = W[k].w ;
        h [k] = W[k].h ;
        pi[k] = W[k].pi;
    }
    GrB_TRY (GrB_Vector_new(pd_output,  GrB_FP64,   n));
    GrB_TRY (GrB_Vector_new(ppi_output, GrB_UINT64, n));
    GrB_TRY (GrB_Vector_new(ph_output,  GrB_UINT64, n));
    GrB_TRY (GrB_Vector_build (*pd_output , I, w , n, GrB_MIN_FP64  ));
    GrB_TRY (GrB_Vector_build (*ppi_output, I, pi, n, GrB_MIN_UINT64));
    GrB_TRY (GrB_Vector_build (*ph_output , I, h , n, GrB_MIN_UINT64));
    LAGraph_FREE_WORK;
    return (GrB_SUCCESS) ;
}
