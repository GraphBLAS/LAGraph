//------------------------------------------------------------------------------
// LAGraph_MMWrite:  write a matrix to a Matrix Market file
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// LAGraph_MMWrite:  write a matrix to a Matrix Market file.
// Contributed by Tim Davis, Texas A&M.

// Writes a matrix to a file in the Matrix Market format.  See LAGraph_mmread
// for a description of the format.

// TODO output is yet not sorted, as required by the Matrix Market format.
// The LAGraph_MMRead can read in an unsorted file, however.

// Parts of this code are from SuiteSparse/CHOLMOD/Check/cholmod_write.c, and
// are used here by permission of the author of CHOLMOD/Check (T. A. Davis).

#include "LG_internal.h"

#undef  LAGraph_FREE_WORK
#define LAGraph_FREE_WORK           \
{                                   \
    LAGraph_Free ((void **) &I) ;   \
    LAGraph_Free ((void **) &J) ;   \
    LAGraph_Free ((void **) &K) ;   \
    LAGraph_Free ((void **) &X) ;   \
    GrB_free (&AT) ;                \
    GrB_free (&M) ;                 \
    GrB_free (&C) ;                 \
}

#undef  LAGraph_FREE_ALL
#define LAGraph_FREE_ALL LAGraph_FREE_WORK

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
    int64_t i, dest = 0, src = 0 ;
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

#if 0
#define FPRINTF(f,...)                  \
{                                       \
    if (fprintf (f, __VA_ARGS__) < 0)   \
    {                                   \
        /* file I/O error */            \
        printf ("LAGraph_MMWrite: file I/O error\n") ; \
        LAGraph_FREE_ALL ;              \
        return (GrB_INVALID_VALUE) ;    \
    }                                   \
}
#endif

#define FPRINTF(f,...)                  \
{                                       \
    LG_CHECK (fprintf (f, __VA_ARGS__) < 0, -2, "Unable to write to file") ; \
}

//------------------------------------------------------------------------------
// LAGraph_MMWrite
//------------------------------------------------------------------------------

int LAGraph_MMWrite
(
    GrB_Matrix A,           // matrix to write to the file
    FILE *f,                // file to write it to, must be already open
    // TODO , FILE *fcomments         // optional file with extra comments
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    void *X = NULL ;
    GrB_Index *I = NULL, *J = NULL, *K = NULL ;
    GrB_Matrix M = NULL, AT = NULL, C = NULL ;
    LG_CHECK (A == NULL || f == NULL, -1005, "invalid inputs") ;

    //--------------------------------------------------------------------------
    // determine the basic matrix properties
    //--------------------------------------------------------------------------

    GrB_Type type ;
    GrB_Index nrows, ncols, nvals ;

    GrB_TRY (GxB_Matrix_type  (&type,  A)) ;
    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
    GrB_Index n = nrows ;

    //--------------------------------------------------------------------------
    // determine if the matrix is dense
    //--------------------------------------------------------------------------

    MM_fmt_enum MM_fmt = MM_coordinate ;

/*  TODO this requires column-major order
    // guard against integer overflow
    if (((double) nrows * (double) ncols < (double) INT64_MAX) &&
        (nvals == nrows * ncols))
    {
        MM_fmt = MM_array ;
    }
*/

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
    else if (type == GrB_FP32 || type == GrB_FP64)
    {
        MM_type = MM_real ;
    }
    #if 0
    else if (type == GxB_FC32 || type == GxB_FC64)
    {
        MM_type = MM_complex ;
    }
    #endif
    else
    {
        // type not supported
        LG_CHECK (false, -1006, "unsupported matrix type") ;
    }

    //--------------------------------------------------------------------------
    // determine symmetry
    //--------------------------------------------------------------------------

    MM_storage_enum MM_storage = MM_general ;

    if (nrows == ncols)
    {
        // AT = A'
        GrB_TRY (GrB_Matrix_new (&AT, type, n, n)) ;
        GrB_TRY (GrB_transpose (AT, NULL, NULL, A, NULL)) ;

        //----------------------------------------------------------------------
        // check for symmetry
        //----------------------------------------------------------------------

        bool isequal = false ;
        LAGraph_TRY (LAGraph_IsEqual (&isequal, A, AT, NULL, msg)) ;
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
            GrB_UnaryOp op = NULL ;
            if      (type == GrB_INT8 ) op = GrB_AINV_INT8  ;
            else if (type == GrB_INT16) op = GrB_AINV_INT16 ;
            else if (type == GrB_INT32) op = GrB_AINV_INT32 ;
            else if (type == GrB_INT64) op = GrB_AINV_INT64 ;
            else if (type == GrB_FP32 ) op = GrB_AINV_FP32  ;
            else if (type == GrB_FP64 ) op = GrB_AINV_FP64  ;
            #if 0
            else if (type == GxB_FC32 ) op = GxB_AINV_FC32 ;
            else if (type == GxB_FC64 ) op = GxB_AINV_FC64 ;
            #endif
            if (op != NULL)
            {
                GrB_TRY (GrB_apply (AT, NULL, NULL, op, AT, NULL)) ;
                LAGraph_TRY (LAGraph_IsEqual (&isequal, A, AT, NULL, msg)) ;
                if (isequal)
                {
                    MM_storage = MM_skew_symmetric ;
                }
            }
        }

        //----------------------------------------------------------------------
        // check for Hermitian
        //----------------------------------------------------------------------

        #if 0
        if (MM_type == MM_complex && MM_storage == MM_general)
        {
            LAGraph_TRY (LAGraph_isall (&isequal, A, AT,
                LAGraph_HERMITIAN_ComplexFP64)) ;
            if (isequal)
            {
                MM_storage = MM_hermitian ;
            }
        }
        #endif

        GrB_free (&AT) ;
    }

    //--------------------------------------------------------------------------
    // determine if the matrix is pattern-only
    //--------------------------------------------------------------------------

    bool is_pattern = false ;
    if (! (MM_storage == MM_skew_symmetric || MM_storage == MM_hermitian))
    {
        if (type == GrB_BOOL)
        {
            LAGraph_TRY (GrB_reduce (&is_pattern, NULL, GrB_LAND_MONOID_BOOL,
                A, NULL)) ;
        }
        else
        {
            GrB_TRY (GrB_Matrix_new (&C, GrB_BOOL, nrows, ncols)) ;
            GrB_BinaryOp op = NULL ;
            if      (type == GrB_INT8  ) op = GrB_EQ_INT8   ;
            else if (type == GrB_INT16 ) op = GrB_EQ_INT16  ;
            else if (type == GrB_INT32 ) op = GrB_EQ_INT32  ;
            else if (type == GrB_INT64 ) op = GrB_EQ_INT64  ;
            else if (type == GrB_UINT8 ) op = GrB_EQ_UINT8  ;
            else if (type == GrB_UINT16) op = GrB_EQ_UINT16 ;
            else if (type == GrB_UINT32) op = GrB_EQ_UINT32 ;
            else if (type == GrB_UINT64) op = GrB_EQ_UINT64 ;
            else if (type == GrB_FP32  ) op = GrB_EQ_FP32   ;
            else if (type == GrB_FP64  ) op = GrB_EQ_FP64   ;
            #if 0
            else if (type == GxB_FC32  ) op = GrB_EQ_FC32 ;
            else if (type == GxB_FC64  ) op = GrB_EQ_FC64 ;
            #endif
            GrB_TRY (GrB_apply (C, NULL, NULL, op, A, 1, NULL)) ;
            GrB_TRY (GrB_reduce (&is_pattern, NULL, GrB_LAND_MONOID_BOOL,
                C, NULL)) ;
            GrB_free (&C) ;
        }
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
    if      (type == GrB_BOOL  ) FPRINTF (f, "GrB_BOOL\n")
    else if (type == GrB_INT8  ) FPRINTF (f, "GrB_INT8\n")
    else if (type == GrB_INT16 ) FPRINTF (f, "GrB_INT16\n")
    else if (type == GrB_INT32 ) FPRINTF (f, "GrB_INT32\n")
    else if (type == GrB_INT64 ) FPRINTF (f, "GrB_INT64\n")
    else if (type == GrB_UINT8 ) FPRINTF (f, "GrB_UINT8\n")
    else if (type == GrB_UINT16) FPRINTF (f, "GrB_UINT16\n")
    else if (type == GrB_UINT32) FPRINTF (f, "GrB_UINT32\n")
    else if (type == GrB_UINT64) FPRINTF (f, "GrB_UINT64\n")
    else if (type == GrB_FP32  ) FPRINTF (f, "GrB_FP32\n")
    else if (type == GrB_FP64  ) FPRINTF (f, "GrB_FP64\n")
    #if 0
    else if (type == GxB_FC32  ) FPRINTF (f, "GxB_FC32\n")
    else if (type == GxB_FC64  ) FPRINTF (f, "GxB_FC64\n")
    #endif

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
        GrB_TRY (GrB_Matrix_new (&C, type, nrows, ncols)) ;
        GrB_TRY (GxB_select (C, NULL, NULL, GxB_DIAG, A, NULL, NULL)) ;
        GrB_Index ndiag = 0 ;
        GrB_TRY (GrB_Matrix_nvals (&ndiag, C)) ;
        GrB_free (&C) ;
        // nvals_to_print = # of entries in tril(A), including diagonal
        nvals_to_print = ndiag + (nvals - ndiag) / 2 ;
    }

    FPRINTF (f, "%" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
        nrows, ncols, nvals_to_print) ;

    if (nvals_to_print == 0)
    {
        // quick return if nothing more to do
        LAGraph_FREE_ALL ;
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // extract and print tuples
    //--------------------------------------------------------------------------

    I = LAGraph_Malloc (nvals, sizeof (GrB_Index)) ;
    J = LAGraph_Malloc (nvals, sizeof (GrB_Index)) ;
    K = LAGraph_Malloc (nvals, sizeof (GrB_Index)) ;
    LG_CHECK (I == NULL || J == NULL || K == NULL, -1004, "out of memory") ;
    for (int64_t k = 0 ; k < nvals ; k++)
    {
        K [k] = k ;
    }
    int nthreads ;
    LAGraph_TRY (LAGraph_GetNumThreads (&nthreads, msg)) ;

    GrB_Index nvals_printed = 0 ;
    bool coord = (MM_fmt == MM_coordinate) ;

    #define WRITE_TUPLES(ctype,is_unsigned,is_signed,is_real,is_complex)    \
    {                                                                       \
        ctype *X = NULL ;                                                   \
        X = LAGraph_Malloc (nvals, sizeof (ctype)) ;                        \
        LG_CHECK (X == NULL, -1010, "out of memory") ;                      \
        GrB_TRY (GrB_Matrix_extractTuples (I, J, ARG(X), &nvals, A)) ;      \
        LAGraph_TRY (LAGraph_Sort3 ((int64_t *) J, (int64_t *) I,           \
            (int64_t *) K, nvals, nthreads, msg)) ;                         \
        for (int64_t k = 0 ; k < nvals ; k++)                               \
        {                                                                   \
            /* convert the row and column index to 1-based */               \
            GrB_Index i = I [k] + 1 ;                                       \
            GrB_Index j = J [k] + 1 ;                                       \
            ctype     x = X [K [k]] ;                                       \
            if (is_general || i >= j)                                       \
            {                                                               \
                /* print the row and column index of the tuple */           \
                if (coord) FPRINTF (f, "%" PRIu64 " %" PRIu64 " ", i, j) ;  \
                /* print the value of the tuple */                          \
                if (is_pattern)                                             \
                {                                                           \
                    /* print nothing */ ;                                   \
                }                                                           \
                else if (is_unsigned)                                       \
                {                                                           \
                    FPRINTF (f, "%" PRIu64, (uint64_t) x) ;                 \
                }                                                           \
                else if (is_signed)                                         \
                {                                                           \
                    FPRINTF (f, "%" PRId64, (int64_t) x) ;                  \
                }                                                           \
                else if (is_real)                                           \
                {                                                           \
                    LG_CHECK (!print_double (f, (double) x), -1002,         \
                        "Unable to write to file") ;                        \
                }                                                           \
            /*  else if (is_complex)                                 */     \
            /*  {                                                    */     \
            /*      LG_CHECK (!print_double (f, creal (x)), -1002,   */     \
            /*          "Unable to write to file") ;                 */     \
            /*      FPRINTF (f, " ") ;                               */     \
            /*      LG_CHECK (!print_double (f, cimag (x)), -1002,   */     \
            /*          "Unable to write to file") ;                 */     \
            /*  }                                                    */     \
                FPRINTF (f, "\n") ;                                         \
            }                                                               \
            nvals_printed++ ;                                               \
        }                                                                   \
        LAGraph_Free ((void **) &X) ;                                       \
    }

    #define ARG(X) X

    if      (type == GrB_BOOL   ) WRITE_TUPLES (bool    , 1, 0, 0, 0)
    else if (type == GrB_INT8   ) WRITE_TUPLES (int8_t  , 0, 1, 0, 0)
    else if (type == GrB_INT16  ) WRITE_TUPLES (int16_t , 0, 1, 0, 0)
    else if (type == GrB_INT32  ) WRITE_TUPLES (int32_t , 0, 1, 0, 0)
    else if (type == GrB_INT64  ) WRITE_TUPLES (int64_t , 0, 1, 0, 0)
    else if (type == GrB_UINT8  ) WRITE_TUPLES (uint8_t , 1, 0, 0, 0)
    else if (type == GrB_UINT16 ) WRITE_TUPLES (uint16_t, 1, 0, 0, 0)
    else if (type == GrB_UINT32 ) WRITE_TUPLES (uint32_t, 1, 0, 0, 0)
    else if (type == GrB_UINT64 ) WRITE_TUPLES (uint64_t, 1, 0, 0, 0)
    else if (type == GrB_FP32   ) WRITE_TUPLES (float   , 0, 0, 1, 0)
    else if (type == GrB_FP64   ) WRITE_TUPLES (double  , 0, 0, 1, 0)
    #if 0
    else if (type == GxB_FC32   ) WRITE_TUPLES (GxB_FC32_t, 0, 0, 0, 1) ;
    else if (type == GxB_FC64   ) WRITE_TUPLES (GxB_FC64_t, 0, 0, 0, 1) ;
    #endif

    ASSERT (nvals_to_print == nvals_printed) ;

    //--------------------------------------------------------------------------
    // free workspace and return
    //--------------------------------------------------------------------------

    LAGraph_FREE_ALL ;
    return (0) ;
}

