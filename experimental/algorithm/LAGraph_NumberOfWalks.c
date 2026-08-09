#include "GraphBLAS.h"
#include "LAGraphX.h"

// Inner recursive helper: computes C = A^k using binary exponentiation.
// Uses LAGraph_plus_one_int64 (structural semiring) only for A*A (k==2),
// where both operands are still the original 0/1 adjacency matrix
// All other squarings use PLUS_TIMES because T = A^(k/2) has walk counts
static GrB_Info NumberOfWalks_inner(GrB_Matrix *C, GrB_Matrix A, int64_t k)
{
    if (C == NULL || A == NULL || k < 0) return GrB_INVALID_VALUE ;

    GrB_Info info ;
    GrB_Index n ;
    GrB_Matrix_nrows (&n, A) ;

    // --- BASE CASES ---
    if (k == 0)
    {
        GrB_Matrix_new (C, GrB_INT64, n, n) ;
        for (GrB_Index i = 0 ; i < n ; i++)
            GrB_Matrix_setElement_INT64 (*C, 1, i, i) ;
        return GrB_SUCCESS ;
    }

    if (k == 1) return GrB_Matrix_dup (C, A) ;

    // k==2: adjacency matrix, 
    if (k == 2)
    {
        GrB_Matrix_new (C, GrB_INT64, n, n) ;
        return GrB_mxm (*C, NULL, NULL, LAGraph_plus_one_int64, A, A, NULL) ;
    }

    // --- RECURSION (binary exponentiation) ---
    GrB_Matrix T = NULL ;
    info = NumberOfWalks_inner (&T, A, k / 2) ;
    if (info != GrB_SUCCESS) return info ;

    // T = A^(k/2) 
    GrB_Matrix_new (C, GrB_INT64, n, n) ;
    GrB_mxm (*C, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, T, T, NULL) ;
    GrB_free (&T) ;

    // If k is odd, multiply by one more A: *C = *C * A
    if (k % 2 != 0)
    {
        GrB_mxm (*C, NULL, NULL, LAGraph_plus_first_int64, *C, A, NULL) ;
    }

    return GrB_SUCCESS ;
}

/**
 * LAGraph_NumberOfWalks: compute number of walks of length k.
 *
 * If src is NULL, computes the full n×n matrix A^k where C(i,j) is the
 * number of distinct walks of length k from node i to node j.
 *
 * If src is a non-NULL indicator vector (1 at the source node index),
 * computes only the walks from that source to all destinations.
 * Result is stored in *C as a 1×n matrix where C(0,j) = walks from src to j.
 */
GrB_Info LAGraph_NumberOfWalks
(
    GrB_Matrix *C,      // output: A^k (n×n) or single-source walks (1×n)
    GrB_Matrix  A,      // input: adjacency matrix
    GrB_Vector  src,    // input: source indicator vector (NULL = all pairs)
    int64_t     k       // input: walk length
)
{
    if (C == NULL || A == NULL || k < 0) return GrB_INVALID_VALUE ;

    // All-pairs case: compute full A^k
    if (src == NULL)
    {
        return NumberOfWalks_inner (C, A, k) ;
    }

    // Single-source case: compute A^k, then w = src * A^k via vxm
    GrB_Info info ;
    GrB_Index n ;
    GrB_Matrix_nrows (&n, A) ;

    // w represents the counts of walks of length 't' starting from 'src'
    GrB_Vector w = NULL ;
    GrB_Vector_dup (&w, src) ;

    for (int64_t t = 0 ; t < k ; t++)
    {
        // w^T = w^T * A: walks of length t+1 from source to each vertex
        info = GrB_vxm (w, NULL, NULL, LAGraph_plus_first_int64, w, A, NULL) ;
        if (info != GrB_SUCCESS) { GrB_free (&w) ; return info ; }
    }

    GrB_Matrix_new (C, GrB_INT64, 1, n) ;
    info = GrB_assign (*C, NULL, NULL, w, 0, GrB_ALL, n, NULL) ;

    GrB_free (&w) ;

    return info;
}
