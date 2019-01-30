//------------------------------------------------------------------------------
// LAGraph_mmwrite:  write a matrix to a Matrix Market file
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Writes a matrix to a file in the Matrix Market format.  See LAGraph_mmread
// for a description of the format.

// Parts of this code are from SuiteSparse/CHOLMOD/Check/cholmod_write.c, and
// are used here by permission of the author of CHOLMOD/Check (T. A. Davis).

#include "LAGraph_internal.h"

//------------------------------------------------------------------------------
// TODO: include_comments
//------------------------------------------------------------------------------

// Read in the comments file, if it exists, and copy it to the Matrix Market
// file.  A "%" is prepended to each line.  Returns true if successful, false
// if an I/O error occu.

#if 0
static bool include_comments
(
    FILE *f,
    const char *comments
)
{
    FILE *cf = NULL ;
    char buffer [MAXLINE] ;
    int ok = TRUE ;
    if (comments != NULL && comments [0] != '\0')
    {
        cf = fopen (comments, "r") ;
        if (cf == NULL)
        {
            return (FALSE) ;
        }
        while (ok && fgets (buffer, MAXLINE, cf) != NULL)
        {
            // ensure the line is not too long
            buffer [MMLEN-1] = '\0' ;
            buffer [MMLEN-2] = '\n' ;
            ok = ok && (fprintf (f, "%%%s", buffer) > 0) ;
        }
        fclose (cf) ;
    }
    return (ok) ;
}
#endif

//------------------------------------------------------------------------------
// print_double
//------------------------------------------------------------------------------

// Print a double value to the file, using the shortest format that ensures the
// value is written precisely.  Returns true if successful, false if an I/O
// error occurred.

static bool print_double
(
    FILE *f,        // file to print to
    double x        // value to print
)
{

    char s [MAXLINE], *p ;
    Int i, dest = 0, src = 0 ;
    int width, ok ;

    //--------------------------------------------------------------------------
    // handle Inf and NaN
    //--------------------------------------------------------------------------

    if (isnan (x))
    {
        return (fprintf (f, "nan") > 0) ;
    }
    if (isinf (x))
    {
        return (fprintf (f, (x < 0) ? "-inf" : "inf") > 0) ;
    }

    //--------------------------------------------------------------------------
    // find the smallest acceptable precision
    //--------------------------------------------------------------------------

    for (width = 6 ; width < 20 ; width++)
    {
        double y ;
        sprintf (s, "%.*g", width, x) ;
        sscanf (s, "%lg", &y) ;
        if (x == y) break ;
    }

    //--------------------------------------------------------------------------
    // shorten the string
    //--------------------------------------------------------------------------

    // change "e+0" to "e", change "e+" to "e", and change "e-0" to "e-"
    for (i = 0 ; i < MAXLINE && s [i] != '\0' ; i++)
    {
        if (s [i] == 'e')
        {
            if (s [i+1] == '+')
            {
                dest = i+1 ;
                if (s [i+2] == '0')
                {
                    // delete characters s[i+1] and s[i+2]
                    src = i+3 ;
                }
                else
                {
                    // delete characters s[i+1]
                    src = i+2 ;
                }
            }
            else if (s [i+1] == '-')
            {
                dest = i+2 ;
                if (s [i+2] == '0')
                {
                    // delete character s[i+2]
                    src = i+3 ;
                }
                else
                {
                    // no change
                    break ;
                }
            }
            while (s [src] != '\0')
            {
                s [dest++] = s [src++] ;
            }
            s [dest] = '\0' ;
            break ;
        }
    }

    // delete the leading "0" if present and not necessary
    p = s ;
    s [MAXLINE-1] = '\0' ;
    i = strlen (s) ;
    if (i > 2 && s [0] == '0' && s [1] == '.')
    {
        // change "0.x" to ".x"
        p = s + 1 ;
    }
    else if (i > 3 && s [0] == '-' && s [1] == '0' && s [2] == '.')
    {
        // change "-0.x" to "-.x"
        s [1] = '-' ;
        p = s + 1 ;
    }

#if 0
    // double-check
    i = sscanf (p, "%lg", &z) ;
    if (i != 1 || y != z)
    {
        // oops! something went wrong in the "e+0" edit, above.
        // this "cannot" happen
        sprintf (s, "%.*g", width, x) ;
        p = s ;
    }
#endif

    //--------------------------------------------------------------------------
    // print the value to the file
    //--------------------------------------------------------------------------

    return (fprintf (f, "%s", p) > 0) ;
}

//------------------------------------------------------------------------------
// FPRINTF: fprintf and check result
//------------------------------------------------------------------------------

#define FPRINTF(f,...)                  \
{                                       \
    if (fprintf (f, __VA_ARGS__) < 0)   \
    {                                   \
        /* file I/O error */            \
        LAGRAPH_FREE_ALL ;              \
        return (GrB_INVALID_VALUE) ;    \
    }                                   \
}

//------------------------------------------------------------------------------
// LAGraph_mmwrite
//------------------------------------------------------------------------------

GrB_info LAGraph_mmwrite
(
    GrB_Matrix A,           // matrix to write to the file
    FILE *f                 // file to write it to
    // TODO , FILE *fcomments         // optional file with extra comments
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (A == NULL || f == NULL)
    {
        // input arguments invalid
        return (GrB_NULL_POINTER) ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Index *I = NULL ;
    GrB_Index *J = NULL ;
    void      *X = NULL ;
    GrB_Matrix M = NULL, AT = NULL ;

    #define LAGRAPH_FREE_ALL        \
    {                               \
        LAGraph_free (&I) ;         \
        LAGraph_free (&J) ;         \
        LAGraph_free (&X) ;         \
        GrB_free (&AT) ;            \
        GrB_free (&M) ;             \
    }

    //--------------------------------------------------------------------------
    // determine the basic matrix properties
    //--------------------------------------------------------------------------

    GrB_Type type ;
    GrB_Index nrows, ncols, nvals ;

    LAGRAPH_OK (GxB_Matrix_type  (&type,  A)) ;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&colss, A)) ;
    LAGRAPH_OK (GrB_Matrix_nvals (&nvals, A)) ;
    GrB_Index n = nrows ;

    //--------------------------------------------------------------------------
    // determine if the matrix is dense
    //--------------------------------------------------------------------------

    MM_fmt_enum MM_fmt = MM_coordinate ;

    // guard against integer overflow
    if (((double) nrows * (double) ncols < (double) INT64_MAX) &&
        (nvals == nrows * ncols))
    {
        MM_fmt = MM_array ;
    }

    //--------------------------------------------------------------------------
    // determine the entry type
    //--------------------------------------------------------------------------

    MM_type_enum MM_type = MM_integer ;

    if (type == GrB_BOOL   || type == GrB_INT8   || type == GrB_INT16  ||
        type == GrB_INT32  || type == GrB_INT64  || type == GrB_UINT8  ||
        type == GrB_UINT16 || type == GrB_UINT32 || type == GrB_UINT64)
    {
        MM_type = MM_integer ;
    }
    if (type == GrB_FP32 || type == GrB_FP64)
    {
        MM_type = MM_real ;
    }
    else if (type == LAGraph_Complex)
    {
        MM_type = MM_complex ;
    }
    else
    {
        // type not supported
        return (GrB_INVALID_VALUE) ;
    }

    //--------------------------------------------------------------------------
    // determine symmetry
    //--------------------------------------------------------------------------

    MM_storage_enum MM_storage = MM_general ;

    if (nrows == ncols)
    {
        // AT = A'
        LAGRAPH_OK (GrB_Matrix_new (&AT, type, n, n)) ;
        LAGRAPH_OK (GrB_transpose (AT, NULL, NULL, A, NULL)) ;

        //----------------------------------------------------------------------
        // check for symmetry
        //----------------------------------------------------------------------

        bool isequal = false ;
        LAGRAPH_OK (LAGraph_isequal (&isequal, A, AT, LAGraph_EQ_Complex)) ;
        if (isequal)
        {
            MM_storage = MM_symmetric ;
        }

        //----------------------------------------------------------------------
        // check for skew-symmetry
        //----------------------------------------------------------------------

        // for signed types only

        if (MM_storage == MM_general)
        {

            // select the operator
            GrB_BinaryOp op = NULL ;
            if      (type == GrB_INT8       ) op = LAGraph_SKEW_INT8 ;
            else if (type == GrB_INT16      ) op = LAGraph_SKEW_INT16 ;
            else if (type == GrB_INT32      ) op = LAGraph_SKEW_INT32 ;
            else if (type == GrB_INT64      ) op = LAGraph_SKEW_INT64 ;
            else if (type == GrB_FP32       ) op = LAGraph_SKEW_FP32   ;
            else if (type == GrB_FP64       ) op = LAGraph_SKEW_FP64   ;
            else if (type == LAGraph_Complex) op = LAGraph_SKEW_Complex ;

            if (op != NULL)
            {
                LAGRAPH_OK (LAGraph_isall (&isequal, A, AT, op)) ;
                if (isequal)
                {
                    MM_storage = MM_skew_symmetric ;
                }
            }
        }

        //----------------------------------------------------------------------
        // check for Hermitian
        //----------------------------------------------------------------------

        if (MM_type == MM_complex && MM_storage == MM_general)
        {
            LAGRAPH_OK (LAGraph_isall (&isequal, A, AT, LAGraph_Hermitian)) ;
            if (isequal)
            {
                MM_storage = MM_hermitian ;
            }
        }

        GrB_free (&AT) ;
    }

    //--------------------------------------------------------------------------
    // determine if the matrix is pattern-only
    //--------------------------------------------------------------------------

    bool is_pattern = false ;
    if (! (MM_storage == MM_skew_symmetric || MM_storage == MM_hermitian))
    {
        LAGRAPH_OK (LAGraph_ispattern (&is_pattern, A, LAGraph_ISONE_Complex)) ;
        if (is_pattern)
        {
            MM_type = MM_pattern ;
        }
    }

    //--------------------------------------------------------------------------
    // write the Matrix Market header
    //--------------------------------------------------------------------------

    FPRINTF (f, "%%%%MatrixMarket matrix") ;

    switch (MM_fmt)
    {
        case MM_coordinate      : FPRINTF (f, " coordinate")        ; break ;
        case MM_array           : FPRINTF (f, " array")             ; break ;
    }

    switch (MM_type)
    {
        case MM_real            : FPRINTF (f, " real")              ; break ;
        case MM_integer         : FPRINTF (f, " integer")           ; break ;
        case MM_complex         : FPRINTF (f, " complex")           ; break ;
        case MM_pattern         : FPRINTF (f, " pattern")           ; break ;
    }

    switch (MM_storage)
    {
        case MM_general         : FPRINTF (f, " general\n")         ; break ;
        case MM_symmetric       : FPRINTF (f, " symmetric\n")       ; break ;
        case MM_skew_symmetric  : FPRINTF (f, " skew-symmetric\n")  ; break ;
        case MM_hermitian       : FPRINTF (f, " Hermitian\n")       ; break ;
    }

    FPRINTF (f, "%%%%GraphBLAS ") ;
    if      (type == GrB_BOOL  ) FPRINTF (f, "GrB_BOOL\n") ;
    else if (type == GrB_INT8  ) FPRINTF (f, "GrB_INT8\n") ;
    else if (type == GrB_INT16 ) FPRINTF (f, "GrB_INT16\n") ;
    else if (type == GrB_INT32 ) FPRINTF (f, "GrB_INT32\n") ;
    else if (type == GrB_INT64 ) FPRINTF (f, "GrB_INT64\n") ;
    else if (type == GrB_UINT8 ) FPRINTF (f, "GrB_UINT8\n") ;
    else if (type == GrB_UINT16) FPRINTF (f, "GrB_UINT16\n") ;
    else if (type == GrB_UINT32) FPRINTF (f, "GrB_UINT32\n") ;
    else if (type == GrB_UINT64) FPRINTF (f, "GrB_UINT64\n") ;
    else if (type == GrB_FP32  ) FPRINTF (f, "GrB_FP32\n") ;
    else if (type == GrB_FP64  ) FPRINTF (f, "GrB_FP64\n") ;
    else                         FPRINTF (f, "LAGraph_Complex\n") ;

    //--------------------------------------------------------------------------
    // include any additional comments
    //--------------------------------------------------------------------------

    // TODO: read comments from the file fcomments, until reaching EOF

    //--------------------------------------------------------------------------
    // print the first line
    //--------------------------------------------------------------------------

    bool is_general = (MM_storage == MM_general) ;
    GrB_Index nvals_to_print = nvals ;

    if (!is_general)
    {
        // count the entries on the diagonal
        LAGRAPH_OK (GrB_Matrix_new (&M, GrB_BOOL, n, n)) ;
        for (int64_t k = 0 ; k < n ; k++)
        {
            // M = diagonal matrix, all ones
            LAGRAPH_OK (GrB_Matrix_setElement_BOOL (M, true, k, k)) ;
        }
        // M<M> = A
        LAGRAPH_OK (GrB_assign (&M, M, NULL, A, GrB_ALL, n, GrB_ALL, n, NULL)) ;
        GrB_Index ndiag = 0 ;
        LAGRAPH_OK (GrB_Matrix_nvals (&ndiag, M)) ;
        GrB_free (&M) ;
        // nvals_to_print = # of entries in tril(A), including diagonal 
        nvals_to_print = ndiag + (nvals - ndiag) / 2 ;
    }

    FPRINTF (f, "%" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
        nrows, ncols, nvals_to_print) ;

    if (nvals_to_print == 0)
    {
        // quick return if nothing more to do
        LAGRAPH_FREE_ALL ;
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // extract the and print tuples
    //--------------------------------------------------------------------------

    LAGRAPH_OK (LAGraph_malloc (&I, nvals, sizeof (GrB_Index))) ;
    LAGRAPH_OK (LAGraph_malloc (&J, nvals, sizeof (GrB_Index))) ;

    // TODO need to sort them as well, into column-major order.  This can be
    // done with a few lines in SuiteSparse:GraphBLAS, but requires an explicit
    // sort otherwise.  The ANSI C11 sort could be used, but it would required
    // the [I J X] arrays to be concatenated.

    GrB_Index nvals_printed = 0 ;

    #define WRITE_TUPLES(ctype,is_unsigned,is_signed,is_real,is_complex)    \
    {                                                                       \
        ctype *X = NULL ;                                                   \
        LAGRAPH_OK (LAGraph_malloc (&X, nvals, sizeof (ctype))) ;           \
        LAGRAPH_OK (GrB_extractTuples (&I, &J, &X, &nvals, A)) ;            \
        for (int64_t k = 0 ; k < nvals ; k++)                               \
        {                                                                   \
            /* convert the row and column index to 1-based */               \
            GrB_Index i = I [k] + 1 ;                                       \
            GrB_Index j = J [k] + 1 ;                                       \
            ctype     x = X [k] ;                                           \
            if (is_general || i >= j)                                       \
            {                                                               \
                /* print the row and column index of the tuple */           \
                FPRINTF (f, "%" PRIu64 " %" PRu64, i, j) ;                  \
                /* print the value of the tuple */                          \
                if (is_pattern)                                             \
                {                                                           \
                    /* print nothing */ ;                                   \
                }                                                           \
                else if (is_unsigned)                                       \
                {                                                           \
                    FPRINTF (f, " %" PRIu64, (uint64_t) x) ;                \
                }                                                           \
                else if (is_signed)                                         \
                {                                                           \
                    FPRINTF (f, " %" PRId64, (int64_t) x) ;                 \
                }                                                           \
                else if (is_real)                                           \
                {                                                           \
                    FPRINTF (f, " ") ;                                      \
                    print_double (f, (double) x) ;                          \
                }                                                           \
                else if (is_complex)                                        \
                {                                                           \
                    FPRINTF (f, " ") ;                                      \
                    print_double (f, creal ((double) x)) ;                  \
                    FPRINTF (f, " ") ;                                      \
                    print_double (f, cimag ((double) x)) ;                  \
                }                                                           \
                FPRINTF (f, "\n") ;                                         \
            }                                                               \
            nvals_printed++ ;                                               \
        }                                                                   \
    }

    if      (type == GrB_BOOL   ) WRITE_TUPLES (bool            , 1, 0, 0, 0) ;
    else if (type == GrB_INT8   ) WRITE_TUPLES (int8_t          , 0, 1, 0, 0) ;
    else if (type == GrB_INT16  ) WRITE_TUPLES (int16_t         , 0, 1, 0, 0) ;
    else if (type == GrB_INT32  ) WRITE_TUPLES (int32_t         , 0, 1, 0, 0) ;
    else if (type == GrB_INT64  ) WRITE_TUPLES (int64_t         , 0, 1, 0, 0) ;
    else if (type == GrB_UINT8  ) WRITE_TUPLES (uint8_t         , 1, 0, 0, 0) ;
    else if (type == GrB_UINT16 ) WRITE_TUPLES (uint16_t        , 1, 0, 0, 0) ;
    else if (type == GrB_UINT32 ) WRITE_TUPLES (uint32_t        , 1, 0, 0, 0) ;
    else if (type == GrB_UINT64 ) WRITE_TUPLES (uint64_t        , 1, 0, 0, 0) ;
    else if (type == GrB_FP32   ) WRITE_TUPLES (float           , 0, 0, 1, 0) ;
    else if (type == GrB_FP64   ) WRITE_TUPLES (double          , 0, 0, 1, 0) ;
    else                          WRITE_TUPLES (double complex  , 0, 0, 1, 0) ;

    ASSERT (nvals_to_print == nvals_printed) ;

    //--------------------------------------------------------------------------
    // free workspace and return
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS) ;
}

