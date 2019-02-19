//------------------------------------------------------------------------------
// LAGraph_BF_full.c: implementation of Bellman-Ford method for shortest 
// paths in given graph using GraphBLAS with both distance and parent returned
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// LAGraph_BF performs a Bellman-Ford to find out shortest from given source
// vertex s. The result is vector d where d(k) is the shortest distance from s
// to k.
//------------------------------------------------------------------------------
//TODO: Think about how
#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL               \
{                                      \
    GrB_free(&d);                      \
    GrB_free(&dtmp);                   \
    GrB_free(&Atmp);                   \
    GrB_free(&BF_Tuple3);              \
    GrB_free(&BF_lMIN_Tuple3);         \
    GrB_free(&BF_PLUSrhs_Tuple3;)      \
    GrB_free(&BF_EQ_Tuple3);           \
    GrB_free(&BF_lMIN_Tuple3_Monoid);  \
    GrB_free(&BF_lMIN_PLUSrhs_Tuple3); \
    free(I);                           \
    free(J);                           \
    free(w);                           \
    free(W);                           \
    free(h);                           \
    free(pi);                          \
}

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
BF_Tuple3_struct;

//------------------------------------------------------------------------------
// 2 binary functions, z=f(x,y), where Tuple3xTuple3 -> Tuple3
//------------------------------------------------------------------------------
void BF_lMIN
(
    BF_Tuple3_struct *z,
    const BF_Tuple3_struct *y,
    const BF_Tuple3_struct *x
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

void BF_PLUSrhs
(
    BF_Tuple3_struct *z,
    const BF_Tuple3_struct *y,
    const BF_Tuple3_struct *x
)
{
    z->w = x->w + y->w;
    z->h = x->h + y->h;
    if (x->pi != UINT64_MAX && y->pi != 0)
    {
        z->pi = y->pi;
    }
    else
    {
        z->pi = x->pi;
    }
}

void BF_EQ
(
    BF_Tuple3_struct *z,
    const BF_Tuple3_struct *y,
    const BF_Tuple3_struct *x
)
{
    if (x->w == y->w && x->h == y->h && x->pi == y->pi)
    {
	z->w  = 1;
	z->h  = 1;
	z->pi = 1;
    }
    else
    {
        z->w  = 0;
        z->h  = 0;
        z->pi = 0;
    }
}

/*
* Given a n-by-n adjacency matrix A and a source vertex s.
* If there is no negative-weight cycle reachable from s, return the distances of
* shortest paths from s and parents along the paths as vector d. Otherwise,
* returns d=NULL if there is a negtive-weight cycle.
* pd is given pointer to a GrB_Vector, where the i-th entry is d(s,i), the
*   sum of edges length in the shortest path
* ppi is given pointer to a GrB_Vector, where the i-th entry is pi(i), the
*   parent of i-th vertex in the shortest path
* ph is given pointer to a GrB_Vector, where the i-th entry is h(s,i), the 
*   number of edges from s to i in the shortest path
* A has zeros on diagonal and weights on corresponding entries of edges
* s is given index for source vertex
*/
GrB_Info LAGraph_BF_full
(
    GrB_Vector *pd,             //the pointer to the vector of distance
    GrB_Vector *ppi,            //the pointer to the vector of parent
    GrB_Vector *ph,             //the pointer to the vector of hops
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Index s,          //given index of the source
)
{
    GrB_Info info;
    // tmp vector to store distance vector after n (i.e., V) loops
    GrB_Vector d = NULL; dtmp = NULL;
    GrB_Matrix Atmp = NULL;
    GrB_Type BF_Tuple3;
    
    GrB_BinaryOp BF_lMIN_Tuple3;
    GrB_BinaryOp BF_PLUSrhs_Tuple3;
    GrB_BinaryOp BF_EQ_Tuple3;

    GrB_Monoid BF_lMIN_Tuple3_Monoid;
    GrB_Semiring BF_lMIN_PLUSrhs_Tuple3;
 
    GrB_Index n, nz, *I, *J;    // n = # of row/col, nz = # of nnz in graph
    GrB_Index *h, *pi;
    double *w;
    BF_Tuple3_struct *W;

    if (A == NULL || pd == NULL)
    {
        // required argument is missing
        LAGRAPH_ERROR ("required arguments are NULL", GrB_NULL_POINTER) ;
    }

    LAGRAPH_OK (GrB_Matrix_nrows(&n, A));
    LAGRAPH_OK (GrB_Matrix_nvals(&nz, A));

    if (s >= n || s < 0)
    {
        LAGRAPH_ERROR ("invalid value for source vertex s", GrB_INVALID_VALUE);
    }
    //--------------------------------------------------------------------------
    // create all GrB_Type GrB_BinaryOp GrB_Monoid and GrB_Semiring
    //--------------------------------------------------------------------------
    LAGRAPH_OK (GrB_Type_new(&BF_Tuple3, sizeof(BF_Tuple3_struct)));
    
    LAGRAPH_OK (GrB_BinaryOp_new(&BF_EQ_Tuple3, &BF_EQ, BF_Tuple3,
        BF_Tuple3, BF_Tuple3));
    LAGRAPH_OK (GrB_BinaryOp_new(&BF_lMIN_Tuple3, &BF_lMIN, BF_Tuple3,
        BF_Tuple3, BF_Tuple3));
    LAGRAPH_OK (GrB_BinaryOp_new(&BF_PLUSrhs_Tuple3, &BF_PLUSrhs, BF_Tuple3,
        BF_Tuple3, BF_Tuple3)); 

    BF_Tuple3_struct BF_identity = (BF_Tuple3_struct) { .w = INFINITY,
        .h = UINT64_MAX, .pi = UINT64_MAX };
    LAGRAPH_OK (GrB_Monoid_new_UDT(&BF_lMIN_Tuple3_Monoid, BF_lMIN_Tuple3,
        &BF_identity));
    LAGRAPH_OK (GrB_Semiring_new(&BF_lMIN_PLUSrhs_Tuple3,
        BF_lMIN_Tuple3_Monoid, BF_PLUSrhs_Tuple3));

    //--------------------------------------------------------------------------
    // malloc for arrays used for tuplets
    //--------------------------------------------------------------------------
    I = LAGraph_malloc(nz, sizeof(GrB_Index));
    J = LAGraph_malloc(nz, sizeof(GrB_Index));
    w = LAGraph_malloc(nz, sizeof(double));
    W = LAGraph_malloc(nz, sizeof(BF_Tuple3_struct));
    LAGRAPH_OK (GrB_Matrix_extractTuples_FP64(I, J, w, &nz, A));
    for (GrB_Index k = 0; k < nz; k++)
    {
        if (w == 0)             //diagonal entries
        {   
            W[k] = (Tuple3) { .w = 0, .h = 0, .pi = 0 };
        }
        else
        {   
            W[k] = (Tuple3) { .w = w, .h = 1, .pi = I[k] + 1 };
        }
    }
    LAGRAPH_OK (GrB_Matrix_new(&Atmp, BF_Tuple3, n, n));
    LAGRAPH_OK (GrB_Matrix_build_UDT(Atmp, I, J, W, nz, BF_lMIN_Tuple3));

    //initialize "distance" vector
    LAGRAPH_OK (GrB_Vector_new(&d, BF_Tuple3, n));
    // initial distance from s to itself
    BF_Tuple3_struct d0 = (BF_Tuple3_struct) { .w = 0, .h = 0, .pi = 0 };
    LAGRAPH_OK (GrB_Vector_setElement_UDT(d, &d0, s));

    // copy d to dtmp in order to create a same size of vector
    LAGRAPH_OK (GrB_Vector_dup(&dtmp, d));
    bool same= false;           // variable indicating if d=dtmp
    int32_t count = 0;          // number of iterations

    // terminate when no new path is found or more than V-1 loops
    while (!same && count < n - 1)
    {
        // execute semiring on d and A, and save the result to dtmp
        LAGRAPH_OK (GrB_mxv(dtmp, GrB_NULL, GrB_NULL, GxB_lMIN_PLUSrhs_Tuple3, 
            Atmp, d, GrB_NULL));
	LAGRAPH_OK (LAGraph_Vector_isequal(&same, dtmp, d, BF_EQ_Tuple3));
        if (!same)
        {
            GrB_Vector ttmp = dtmp;
            dtmp = d;
            d = ttmp;
        }
        count++;
    }

    // check for negative-weight cycle only when there was a new path in the  
    // last loop, otherwise, there can't be a negative-weight cycle.
    if (!same)
    {
        // execute semiring again to check for negative-weight cycle
        LAGRAPH_OK (GrB_mxv(dtmp, GrB_NULL, GrB_NULL, GxB_lMIN_PLUSrhs_Tuple3, 
            Atmp, d, GrB_NULL));
	LAGRAPH_OK (LAGraph_Vector_isequal(&same, dtmp, d, BF_EQ_Tuple3));
        
	// if d != dtmp, then there is a negative-weight cycle in the graph
        if (!same)
        {
            printf("A negative-weight cycle found. \n");
            *pd = NULL;
	    *ppi = NULL;
	    *ph = NULL;
	    LAGRAPH_FREE_ALL;
            return (GrB_SUCCESS) ;
        }
    }
    LAGRAPH_OK (GrB_Vector_extractTuples_FP64(I, W, &nz, d));
    h = LAGraph_malloc(nz, sizeof(GrB_Index));
    pi = LAGraph_malloc(nz, sizeof(GrB_Index));

    for (GrB_Index k = 0; k < nz; k++)
    {
	w [k] = W[k].w ;
	h [k] = W[k].h ;
	pi[k] = W[k].pi;
    }
    LAGRAPH_OK (GrB_Vector_new(pd,  GrB_FP64,   n));
    LAGRAPH_OK (GrB_Vector_new(ppi, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new(ph,  GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_build_FP64  (*pd , I, w , nz, GrB_MIN_FP64  ));
    LAGRAPH_OK (GrB_Vector_build_UINT64(*ppi, I, pi, nz, GrB_MIN_UINT64));
    LAGRAPH_OK (GrB_Vector_build_UINT64(*ph , I, h , nz, GrB_MIN_UINT64));
    LAGRAPH_FREE_ALL;
    return (GrB_SUCCESS) ;
}
