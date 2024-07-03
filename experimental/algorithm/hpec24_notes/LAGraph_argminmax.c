// A simplified algorithm for HPEC'24
// assume the matrix type is FP64
// assume argmax
// use mxv where appropriate
// don't use the ANY monoid.

//------------------------------------------------------------------------------
// argmax: compute argmax of each row of A
//------------------------------------------------------------------------------

int argmax
(
    // output
    GrB_Vector *x_handle,   // max value in each row of A
    GrB_Vector *p_handle,   // index of max value in each row of A
    // input
    GrB_Matrix A            // assumed to be GrB_FP64
)
{

    //--------------------------------------------------------------------------
    // create outputs x and p, and the iso full vector y
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols ;
    GrB_Matrix_nrows (&nrows, A) ;
    GrB_Matrix_ncols (&ncols, A) ;
    GrB_Vector y = NULL, x = NULL, p = NULL ;
    GrB_Matrix G = NULL, D = NULL ;
    GrB_Vector_new (&x, GrB_FP64, nrows) ;
    GrB_Vector_new (&y, GrB_FP64, ncols) ;
    GrB_Vector_new (&p, GrB_INT64, nrows) ;

    // y (:) = 1, an full vector with all entries equal to 1
    GrB_Matrix_assign_INT64 (y, NULL, NULL, 1, GrB_ALL, ncols, NULL) ;

    //--------------------------------------------------------------------------
    // compute x = max(A)
    //--------------------------------------------------------------------------

    // x = max (A) where x(i) = max (A (i,:))
    GrB_mxv (x, NULL, NULL, GrB_MAX_FIRST_SEMIRING_FP64, A, y, NULL) ;

    //--------------------------------------------------------------------------
    // compute G, where G(i,j)=1 if A(i,j) is the max in its row
    //--------------------------------------------------------------------------

    // D = diag (x)
    GrB_Matrix_diag (&D, x, 0) ;
    GrB_Matrix_new (&G, GrB_BOOL, nrows, ncols) ;
    // G = D*A using the EQ_EQ_FP64 semiring
    GrB_mxm (G, NULL, NULL, GxB_EQ_EQ_FP64, D, A, NULL) ;
    // drop explicit zeros from G
    GrB_Matrix_select_BOOL (G, NULL, NULL, GrB_VALUENE_BOOL, G, 0, NULL) ;

    //--------------------------------------------------------------------------
    // extract the positions of the entries in G
    //--------------------------------------------------------------------------

    // find the position of the max entry in each row:
    // p = G*y, so that p(i) = j if x(i) = A(i,j) = max (A (i,:)).

    if (no 2ndI op)
    {
        // H = rowindex (G)
        GrB_Matrix H = NULL ;
        GrB_Matrix_new (&H, nrows, ncols) ;
        GrB_apply (H, NULL, NULL, GrB_ROWINDEX_INT64, G, NULL) ;
        // p = H*y
        GrB_mxv (p, NULL, NULL, GrB_MIN_FIRST_SEMIRING_INT64, H, y, NULL) ;
        GrB_free (&H) ;
    }
    else
    {
        // using the SECONDI operator
        GrB_mxm (p, NULL, NULL, GxB_MIN_SECONDI_INT64, G, y, NULL) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    GrB_Matrix_free (&D) ;
    GrB_Matrix_free (&G) ;
    GrB_Matrix_free (&y) ;
    (*x_handle) = x ;
    (*p_handle) = p ;
}

