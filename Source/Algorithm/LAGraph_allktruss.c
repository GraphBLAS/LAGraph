//------------------------------------------------------------------------------
// LAGraph_allktruss.c: find all k-trusses of a graph via GraphBLAS
//------------------------------------------------------------------------------

// Given a symmetric graph A with no-self edges, LAGraph_allktruss finds all
// k-trusses of A.

// The edge weights of A are treated as binary.  Explicit zero entries in A are
// treated as non-edges.  Any type will work, but uint32 is recommended for
// fastest results since that is the type used here for the semiring.
// GraphBLAS will do typecasting internally, but that takes extra time. 

// The optional output matrices Cset [3..kmax-1] are the k-trusses of A.  Their
// edges are a subset of A.  Each edge in C = Cset [k] is part of at least k-2
// triangles in C.  The pattern of C is the adjacency matrix of the k-truss
// subgraph of A.  The edge weights of C are the support of each edge.  That
// is, C(i,j)=nt if the edge (i,j) is part of nt triangles in C.  All edges in
// C have support of at least k-2.  The total number of triangles in C is
// sum(C)/6.  The number of edges in C is nnz(C)/2.  C = Cset [k] is returned
// as symmetric with a zero-free diagonal, if Cset is not NULL on input.  The
// k-trusses are not returned if Cset is NULL.  Cset [kmax] is NULL since the
// kmax-truss is empty.

// The arrays ntris, nedges, and nstepss hold the output statistics.
// ntris   [k] = # of triangles in the k-truss
// nedges  [k] = # of edges in the k-truss
// nstepss [k] = # of steps required to compute the k-truss

// Usage: constructs k-trusses of A, for k = 3:kmax

//      GrB_Matrix_nrows (&n, A) ;
//      GrB_Matrix *Cset = LAGraph_malloc (n, sizeof (GrB_Matrix)) ;
//      int64_t *ntris   = LAGraph_malloc (n, sizeof (int64_t)) ;
//      int64_t *nedges  = LAGraph_malloc (n, sizeof (int64_t)) ;
//      int64_t *nstepss = LAGraph_malloc (n, sizeof (int64_t)) ;
//      GrB_Info info = LAGraph_allktruss (&Cset, A, &kmax,
//          ntris, nedges, nstepss) ;

// Compare this function with the MATLAB equivalent, allktruss.m, in
// LAGraph/Test/AllKTruss.  This function is derived from SuiteSparse/
// GraphBLAS/Extras/ktruss/allktruss_graphblas.c

#define LAGRAPH_FREE_ALL                        \
    if (keep_all_ktrusses)                      \
    {                                           \
        for (int64_t kk = 3 ; kk <= k ; kk++)   \
        {                                       \
            GrB_free (&(Cset [kk])) ;           \
        }                                       \
    }                                           \
    GrB_free (&C) ;

#include "LAGraph_internal.h"

//------------------------------------------------------------------------------
// C = LAGraph_allktruss (A,k): find all k-trusses a graph
//------------------------------------------------------------------------------

GrB_Info LAGraph_allktruss      // compute all k-trusses of a graph
(
    GrB_Matrix *Cset,           // size n, output k-truss subgraphs (optional)
    GrB_Matrix A,               // input adjacency matrix, A, not modified
    // output statistics
    int64_t *kmax,              // smallest k where k-truss is empty
    int64_t *ntris,             // size n, ntris [k] is #triangles in k-truss
    int64_t *nedges,            // size n, nedges [k] is #edges in k-truss
    int64_t *nstepss            // size n, nstepss [k] is #steps for k-truss
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (nstepss == NULL || kmax == NULL || ntris == NULL || nedges == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    bool keep_all_ktrusses = (Cset != NULL) ;

    int64_t k ;
    for (k = 0 ; k < 3 ; k++)
    {
        if (keep_all_ktrusses)
        {
            Cset [k] = NULL ;
        }
        ntris   [k] = 0 ;
        nedges  [k] = 0 ;
        nstepss [k] = 0 ;
    }
    (*kmax) = 0 ;
    k = 0 ;

    GrB_Info info ;

    // the current k-truss
    GrB_Matrix C = NULL ;

    // get the size of A
    GrB_Index n ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;

    //--------------------------------------------------------------------------
    // C<A> = A*A
    //--------------------------------------------------------------------------

    GrB_Index last_cnz ;
    LAGRAPH_OK (GrB_Matrix_nvals (&last_cnz, A)) ;       // last_cnz = nnz (A)
    LAGRAPH_OK (GrB_Matrix_new (&C, GrB_UINT32, n, n)) ;
    LAGRAPH_OK (GrB_mxm (C, A, NULL, GxB_PLUS_LAND_UINT32, A, A, NULL)) ;
    int64_t nsteps = 1 ;

    //--------------------------------------------------------------------------
    // find all k-trusses
    //--------------------------------------------------------------------------

    for (k = 3 ; ; k++)
    {

        //----------------------------------------------------------------------
        // find the k-truss
        //----------------------------------------------------------------------

        int64_t support = (k-2) ;

        while (1)
        {

            //------------------------------------------------------------------
            // C = C .* (C >= support)
            //------------------------------------------------------------------

            LAGRAPH_OK (GxB_select (C, NULL, NULL, LAGraph_support, C,
                &support, NULL)) ;

            //------------------------------------------------------------------
            // check if k-truss has been found
            //------------------------------------------------------------------

            GrB_Index cnz ;
            LAGRAPH_OK (GrB_Matrix_nvals (&cnz, C)) ;
            if (cnz == last_cnz)
            {
                // k-truss has been found
                int64_t nt = 0 ;
                LAGRAPH_OK (GrB_reduce (&nt, NULL, GxB_PLUS_INT64_MONOID,
                    C, NULL)) ;
                ntris   [k] = nt / 6 ;
                nedges  [k] = cnz / 2 ;
                nstepss [k] = nsteps ;
                nsteps = 0 ;
                if (cnz == 0)
                {
                    // this is the last k-truss
                    LAGRAPH_OK (GrB_free (&C)) ;    // free last empty k-truss
                    (*kmax) = k ;
                    if (keep_all_ktrusses)
                    {
                        Cset [k] = NULL ;
                    }
                    return (GrB_SUCCESS) ;
                }
                else if (keep_all_ktrusses)
                {
                    // save the k-truss in the list of output k-trusses
                    // TODO: if Cset [k] == Cset [k-1], then do not save it.
                    // Set it to NULLL to denote that the k-truss is the
                    // same as the (k-1)-truss.  Also, advance quickly to
                    // the next k, setting k = min (C).
                    LAGRAPH_OK (GrB_Matrix_dup (&(Cset [k]), C)) ;
                }
                // start finding the next k-truss
                break ;
            }

            // continue searching for this k-truss
            last_cnz = cnz ;
            nsteps++ ;

            //------------------------------------------------------------------
            // C<C> = C*C
            //------------------------------------------------------------------

            LAGRAPH_OK (GrB_mxm (C, C, NULL, GxB_PLUS_LAND_UINT32, C, C, NULL));
        }
    }
}

