//------------------------------------------------------------------------------
// LAGraph_MMRead: read a matrix from a Matrix Market file
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// LAGraph_MMRead: read a matrix from a Matrix Market file

// Parts of this code are from SuiteSparse/CHOLMOD/Check/cholmod_read.c, and
// are used here by permission of the author of CHOLMOD/Check (T. A. Davis).

// TODO: need to decide on error codes to return. Currently:
//  -1001:  an input parameter is NULL
//  -1002:  the contents of the file are invalid in some way
// It would be possible to return a wider range of error codes that explain
// why the contents of the input file are invalid.

#define LAGraph_FREE_ALL GrB_free (A) ;

#include "LG_internal.h"

//------------------------------------------------------------------------------
// get_line
//------------------------------------------------------------------------------

// Read one line of the file, return true if successful, false if EOF.
// The string is returned in buf, converted to lower case.

static inline bool get_line
(
    FILE *f,        // file open for reading
    char *buf       // size MAXLINE+1
)
{

    // check inputs
    ASSERT (f != NULL) ;
    ASSERT (buf != NULL) ;

    // read the line from the file
    buf [0] = '\0' ;
    buf [1] = '\0' ;
    if (fgets (buf, MAXLINE, f) == NULL)
    {
        // EOF or other I/O error
        return (false) ;
    }
    buf [MAXLINE] = '\0' ;

    // convert the string to lower case
    for (int k = 0 ; k < MAXLINE && buf [k] != '\0' ; k++)
    {
        buf [k] = tolower (buf [k]) ;
    }
    return (true) ;
}

//------------------------------------------------------------------------------
// is_blank_line
//------------------------------------------------------------------------------

// returns true if buf is a blank line or comment, false otherwise.

static inline bool is_blank_line
(
    char *buf       // size MAXLINE+1, never NULL
)
{

    // check inputs
    ASSERT (buf != NULL) ;

    // check if comment line
    if (buf [0] == '%')
    {
        // line is a comment
        return (true) ;
    }

    // check if blank line
    for (int k = 0 ; k <= MAXLINE ; k++)
    {
        int c = buf [k] ;
        if (c == '\0')
        {
            // end of line
            break ;
        }
        if (!isspace (c))
        {
            // non-space character; this is not an error
            return (false) ;
        }
    }

    // line is blank
    return (true) ;
}

//------------------------------------------------------------------------------
// read_double
//------------------------------------------------------------------------------

// Read a single double value from a string.  The string may be any string
// recognized by sscanf, or inf, -inf, +inf, or nan.  The token infinity is
// also OK instead of inf (only the first 3 letters of inf* or nan* are
// significant, and the rest are ignored).

static inline bool read_double      // true if successful, false if failure
(
    char *p,        // string containing the value
    double *rval    // value to read in
)
{
    while (*p && isspace (*p)) p++ ;   // skip any spaces

    if ((strncmp (p, "inf", 3) == 0) || (strncmp (p, "+inf", 4) == 0))
    {
        (*rval) = INFINITY ;
    }
    else if (strncmp (p, "-inf", 4) == 0)
    {
        (*rval) = -INFINITY ;
    }
    else if (strncmp (p, "nan", 3) == 0)
    {
        (*rval) = NAN ;
    }
    else
    {
        if (sscanf (p, "%lg", rval) != 1)
        {
            // invalid file format, EOF, or other I/O error
            return (false) ;
        }
    }
    return (true) ;
}

//------------------------------------------------------------------------------
// read_entry
//------------------------------------------------------------------------------

static inline bool read_entry   // returns true if successful, false if failure
(
    char *p,        // string containing the value
    GrB_Type type,  // type of value to read
    bool structural,   // if true, then the value is 1
    char *x         // value read in, a pointer to space of size of the type
)
{

    int64_t ival = 1 ;
    double rval = 1, zval = 0 ;

    while (*p && isspace (*p)) p++ ;   // skip any spaces

    if (type == GrB_BOOL)
    {
        if (!structural && sscanf (p, "%" SCNd64, &ival) != 1) return (false) ;
        if (ival < 0 || ival > 1)
        {
            // entry out of range
            return (false) ;
        }
        bool *result = (bool *) x ;
        result [0] = (bool) ival ;
    }
    else if (type == GrB_INT8)
    {
        if (!structural && sscanf (p, "%" SCNd64, &ival) != 1) return (false) ;
        if (ival < INT8_MIN || ival > INT8_MAX)
        {
            // entry out of range
            return (false) ;
        }
        int8_t *result = (int8_t *) x ;
        result [0] = (int8_t) ival ;
    }
    else if (type == GrB_INT16)
    {
        if (!structural && sscanf (p, "%" SCNd64, &ival) != 1) return (false) ;
        if (ival < INT16_MIN || ival > INT16_MAX)
        {
            // entry out of range
            return (false) ;
        }
        int16_t *result = (int16_t *) x ;
        result [0] = (int16_t) ival ;
    }
    else if (type == GrB_INT32)
    {
        if (!structural && sscanf (p, "%" SCNd64, &ival) != 1) return (false) ;
        if (ival < INT32_MIN || ival > INT32_MAX)
        {
            // entry out of range
            return (false) ;
        }
        int32_t *result = (int32_t *) x ;
        result [0] = (int32_t) ival ;
    }
    else if (type == GrB_INT64)
    {
        if (!structural && sscanf (p, "%" SCNd64, &ival) != 1) return (false) ;
        int64_t *result = (int64_t *) x ;
        result [0] = (int64_t) ival ;
    }
    else if (type == GrB_UINT8)
    {
        if (!structural && sscanf (p, "%" SCNd64, &ival) != 1) return (false) ;
        if (ival < 0 || ival > UINT8_MAX)
        {
            // entry out of range
            return (false) ;
        }
        uint8_t *result = (uint8_t *) x ;
        result [0] = (uint8_t) ival ;
    }
    else if (type == GrB_UINT16)
    {
        if (!structural && sscanf (p, "%" SCNd64, &ival) != 1) return (false) ;
        if (ival < 0 || ival > UINT16_MAX)
        {
            // entry out of range
            return (false) ;
        }
        uint16_t *result = (uint16_t *) x ;
        result [0] = (uint16_t) ival ;
    }
    else if (type == GrB_UINT32)
    {
        if (!structural && sscanf (p, "%" SCNd64, &ival) != 1) return (false) ;
        if (ival < 0 || ival > UINT32_MAX)
        {
            // entry out of range
            return (false) ;
        }
        uint32_t *result = (uint32_t *) x ;
        result [0] = (uint32_t) ival ;
    }
    else if (type == GrB_UINT64)
    {
        uint64_t uval = 1 ;
        if (!structural && sscanf (p, "%" SCNu64, &uval) != 1) return (false) ;
        uint64_t *result = (uint64_t *) x ;
        result [0] = (uint64_t) uval ;
    }
    else if (type == GrB_FP32)
    {
        if (!structural && !read_double (p, &rval)) return (false) ;
        float *result = (float *) x ;
        result [0] = (float) rval ;
    }
    else if (type == GrB_FP64)
    {
        if (!structural && !read_double (p, &rval)) return (false) ;
        double *result = (double *) x ;
        result [0] = rval ;
    }
#if 0
    else if (type == GxB_FC32)
    {
        if (!structural && !read_double (p, &rval)) return (false) ;
        while (*p && !isspace (*p)) p++ ;   // skip real part
        if (!structural && !read_double (p, &zval)) return (false) ;
        float *result = (float *) x ;
        result [0] = (float) rval ;     // real part
        result [1] = (float) zval ;     // imaginary part
    }
    else if (type == GxB_FC64)
    {
        if (!structural && !read_double (p, &rval)) return (false) ;
        while (*p && !isspace (*p)) p++ ;   // skip real part
        if (!structural && !read_double (p, &zval)) return (false) ;
        double *result = (double *) x ;
        result [0] = rval ;     // real part
        result [1] = zval ;     // imaginary part
    }
#endif

    return (true) ;
}

//------------------------------------------------------------------------------
// negate_scalar: negate a scalar value
//------------------------------------------------------------------------------

// negate the scalar x.  Do nothing for bool or uint*.

static inline void negate_scalar
(
    GrB_Type type,
    void *x
)
{

    if (type == GrB_INT8)
    {
        int8_t *value = (int8_t *) x ;
        (*value) = - (*value) ;
    }
    else if (type == GrB_INT16)
    {
        int16_t *value = (int16_t *) x ;
        (*value) = - (*value) ;
    }
    else if (type == GrB_INT32)
    {
        int32_t *value = (int32_t *) x ;
        (*value) = - (*value) ;
    }
    else if (type == GrB_INT64)
    {
        int64_t *value = (int64_t *) x ;
        (*value) = - (*value) ;
    }
    else if (type == GrB_FP32)
    {
        float *value = (float *) x ;
        (*value) = - (*value) ;
    }
    else if (type == GrB_FP64)
    {
        double *value = (double *) x ;
        (*value) = - (*value) ;
    }
#if 0
    else if (type == GxB_FC32)
    {
        float complex *value = (float complex *) x ;
        (*value) = - (*value) ;
    }
    else if (type == GxB_FC64)
    {
        double complex *value = (double complex *) x ;
        (*value) = - (*value) ;
    }
#endif
}

//------------------------------------------------------------------------------
// set_value
//------------------------------------------------------------------------------

// A(i,j) = x using GrB_Matrix_setElement_<type>.  No typecasting is done.

static inline GrB_Info set_value
(
    GrB_Matrix A,
    GrB_Type type,
    GrB_Index i,
    GrB_Index j,
    char *x
)
{

    if (type == GrB_BOOL)
    {
        bool *value = (bool *) x ;
        return (GrB_Matrix_setElement_BOOL (A, *value, i, j)) ;
    }
    else if (type == GrB_INT8)
    {
        int8_t *value = (int8_t *) x ;
        return (GrB_Matrix_setElement_INT8 (A, *value, i, j)) ;
    }
    else if (type == GrB_INT16)
    {
        int16_t *value = (int16_t *) x ;
        return (GrB_Matrix_setElement_INT16 (A, *value, i, j)) ;
    }
    else if (type == GrB_INT32)
    {
        int32_t *value = (int32_t *) x ;
        return (GrB_Matrix_setElement_INT32 (A, *value, i, j)) ;
    }
    else if (type == GrB_INT64)
    {
        int64_t *value = (int64_t *) x ;
        return (GrB_Matrix_setElement_INT64 (A, *value, i, j)) ;
    }
    else if (type == GrB_UINT8)
    {
        uint8_t *value = (uint8_t *) x ;
        return (GrB_Matrix_setElement_UINT8 (A, *value, i, j)) ;
    }
    else if (type == GrB_UINT16)
    {
        uint16_t *value = (uint16_t *) x ;
        return (GrB_Matrix_setElement_UINT16 (A, *value, i, j)) ;
    }
    else if (type == GrB_UINT32)
    {
        uint32_t *value = (uint32_t *) x ;
        return (GrB_Matrix_setElement_UINT32 (A, *value, i, j)) ;
    }
    else if (type == GrB_UINT64)
    {
        uint64_t *value = (uint64_t *) x ;
        return (GrB_Matrix_setElement_UINT64 (A, *value, i, j)) ;
    }
    else if (type == GrB_FP32)
    {
        float *value = (float *) x ;
        return (GrB_Matrix_setElement_FP32 (A, *value, i, j)) ;
    }
    else // if (type == GrB_FP64)
    {
        double *value = (double *) x ;
        return (GrB_Matrix_setElement_FP64 (A, *value, i, j)) ;
    }
#if 0
    else if (type == GxB_FC32)
    {
        float complex *value = (float complex *) x ;
        return (GxB_Matrix_setElement_FC32 (A, *value, i, j)) ;
    }
    else if (type == GxB_FC64)
    {
        double complex *value = (double complex *) x ;
        return (GxB_Matrix_setElement_FC64 (A, *value, i, j)) ;
    }
#endif
}

//------------------------------------------------------------------------------
// LAGraph_mmread
//------------------------------------------------------------------------------

int LAGraph_MMRead          // returns 0 if successful, -1 if faillure
(
    GrB_Matrix *A,          // handle of matrix to create
    GrB_Type   *A_type,     // type of the scalar stored in A
    FILE *f,                // file to read from, already open
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_CHECK (A == NULL || A_type == NULL || f == NULL, -1001,
        "inputs are NULL") ;
    (*A) = NULL ;
    (*A_type) = NULL;

    //--------------------------------------------------------------------------
    // set the default properties
    //--------------------------------------------------------------------------

    MM_fmt_enum     MM_fmt     = MM_coordinate ;
    MM_type_enum    MM_type    = MM_real ;
    MM_storage_enum MM_storage = MM_general ;
    GrB_Type type = GrB_FP64 ;
    GrB_Index nrows = 0 ;
    GrB_Index ncols = 0 ;
    GrB_Index nvals = 0 ;

    //--------------------------------------------------------------------------
    // read the Matrix Market header
    //--------------------------------------------------------------------------

    // Read the header.  This consists of zero or more comment lines (blank, or
    // starting with a "%" in the first column), followed by a single data line
    // containing two or three numerical values.  The first line is normally:
    //
    //          %%MatrixMarket matrix <fmt> <type> <storage>
    //
    // but this is optional.  The 2nd line is also optional (the %%MatrixMarket
    // line is required for this 2nd line to be recognized):
    //
    //          %%GraphBLAS <graphblastype>
    //
    // If the %%MatrixMarket line is not present, then the <fmt> <type> and
    // <storage> are implicit.  If the first data line contains 3 items,
    // then the implicit header is:
    //
    //          %%MatrixMarket matrix coordinate real general
    //          %%GraphBLAS GrB_FP64
    //
    // If the first data line contains 2 items (nrows ncols), then the implicit
    // header is:
    //
    //          %%MatrixMarket matrix array real general
    //          %%GraphBLAS GrB_FP64
    //
    // The implicit header is an extension of the Matrix Market format.

    char buf [MAXLINE+1] ;

    bool got_mm_header = false ;

    for (int64_t line = 1 ; get_line (f, buf) ; line++)
    {

        //----------------------------------------------------------------------
        // parse the line
        //----------------------------------------------------------------------

        if ((line == 1) && (strncmp (buf, "%%matrixmarket", 14) == 0))
        {

            //------------------------------------------------------------------
            // read a Matrix Market header
            //------------------------------------------------------------------

            //  %%MatrixMarket matrix <fmt> <type> <storage>
            //  if present, it must be the first line in the file.

            got_mm_header = true ;
            char *p = buf + 14 ;

            //------------------------------------------------------------------
            // get "matrix" token and discard it
            //------------------------------------------------------------------

            while (*p && isspace (*p)) p++ ;        // skip any leading spaces

            if (strncmp (p, "matrix", 6) != 0)
            {
                // invalid Matrix Market object
                LG_CHECK (true, -1002, "invalid object") ;
            }
            p += 6 ;                                // skip past token "matrix"

            //------------------------------------------------------------------
            // get the fmt token
            //------------------------------------------------------------------

            while (*p && isspace (*p)) p++ ;        // skip any leading spaces

            if (strncmp (p, "coordinate", 10) == 0)
            {
                MM_fmt = MM_coordinate ;
                p += 10 ;
            }
            else if (strncmp (p, "array", 5) == 0)
            {
                MM_fmt = MM_array ;
                p += 5 ;
            }
            else
            {
                // invalid Matrix Market format
                LG_CHECK (true, -1002, "invalid format") ;
            }

            //------------------------------------------------------------------
            // get the Matrix Market type token
            //------------------------------------------------------------------

            while (*p && isspace (*p)) p++ ;        // skip any leading spaces

            if (strncmp (p, "real", 4) == 0)
            {
                MM_type = MM_real ;
                type = GrB_FP64 ;
                p += 4 ;
            }
            else if (strncmp (p, "integer", 7) == 0)
            {
                MM_type = MM_integer ;
                type = GrB_INT64 ;
                p += 7 ;
            }
            else if (strncmp (p, "complex", 7) == 0)
            {
                MM_type = MM_complex ;
#if 0
                type = GxB_FC64 ;
                p += 7 ;
#endif
                LG_CHECK (true, -1, "complex types not yet supported") ;
            }
            else if (strncmp (p, "pattern", 7) == 0)
            {
                MM_type = MM_pattern ;
                type = GrB_BOOL ;
                p += 7 ;
            }
            else
            {
                // invalid Matrix Market type
                LG_CHECK (true, -1002, "invalid type") ;
            }

            //------------------------------------------------------------------
            // get the storage token
            //------------------------------------------------------------------

            while (*p && isspace (*p)) p++ ;        // skip any leading spaces

            if (strncmp (p, "general", 7) == 0)
            {
                MM_storage = MM_general ;
            }
            else if (strncmp (p, "symmetric", 9) == 0)
            {
                MM_storage = MM_symmetric ;
            }
            else if (strncmp (p, "skew-symmetric", 14) == 0)
            {
                MM_storage = MM_skew_symmetric ;
            }
            else if (strncmp (p, "hermitian", 9) == 0)
            {
                MM_storage = MM_hermitian ;
            }
            else
            {
                // invalid Matrix Market storage
                LG_CHECK (true, -1002, "invalid storage") ;
            }

            //------------------------------------------------------------------
            // ensure the combinations are valid
            //------------------------------------------------------------------

            if (MM_type == MM_pattern)
            {
                // (coodinate) x (pattern) x (general or symmetric)
                LG_CHECK (!
                    (MM_fmt == MM_coordinate &&
                    (MM_storage == MM_general || MM_storage == MM_symmetric)),
                    -1002, "invalid pattern combination") ;
            }

            if (MM_storage == MM_hermitian)
            {
                // (coordinate or array) x (complex) x (Hermitian)
                LG_CHECK (! (MM_type == MM_complex), -1002,
                    "invalid complex combination") ;
            }

        }
        else if (got_mm_header && (line == 2)
                 && (strncmp (buf, "%%graphblas", 11) == 0))
        {

            // -----------------------------------------------------------------
            // %%GraphBLAS <entrytype>
            // -----------------------------------------------------------------

            // This must appear as the 2nd line in the file, after the
            // %%MatrixMarket header (which is required in this case; otherwise
            // the %%GraphBLAS line is treated as a pure comment and the
            // <entrytype> is ignored).

            char *p = buf + 11 ;

            while (*p && isspace (*p)) p++ ;        // skip any leading spaces

            // <entrytype> is one of the 11 real built-in types (GrB_BOOL,
            // GrB_INT8, GrB_INT16, GrB_INT32, GrB_INT64, GrB_UINT8,
            // GrB_UINT16, GrB_UINT32, GrB_UINT64, GrB_FP32, GrB_FP64).  The
            // complex types GxB_FC32, GxB_FC64 are not yet supported.

            if (strncmp (p, "grb_bool", 8) == 0)
            {
                type = GrB_BOOL ;
            }
            else if (strncmp (p, "grb_int8", 8) == 0)
            {
                type = GrB_INT8 ;
            }
            else if (strncmp (p, "grb_int16", 9) == 0)
            {
                type = GrB_INT16 ;
            }
            else if (strncmp (p, "grb_int32", 9) == 0)
            {
                type = GrB_INT32 ;
            }
            else if (strncmp (p, "grb_int64", 9) == 0)
            {
                type = GrB_INT64 ;
            }
            else if (strncmp (p, "grb_uint8", 9) == 0)
            {
                type = GrB_UINT8 ;
            }
            else if (strncmp (p, "grb_uint16", 10) == 0)
            {
                type = GrB_UINT16 ;
            }
            else if (strncmp (p, "grb_uint32", 10) == 0)
            {
                type = GrB_UINT32 ;
            }
            else if (strncmp (p, "grb_uint64", 10) == 0)
            {
                type = GrB_UINT64 ;
            }
            else if (strncmp (p, "grb_fp32", 8) == 0)
            {
                type = GrB_FP32 ;
            }
            else if (strncmp (p, "grb_fp64", 8) == 0)
            {
                type = GrB_FP64 ;
            }
#if 0
            else if (strncmp (p, "grb_fc32", 8) == 0)
            {
                type = GxB_FC32 ;
            }
            else if (strncmp (p, "grb_fc64", 8) == 0)
            {
                type = GxB_FC64 ;
            }
#endif
            else
            {
                // type not supported
                LG_CHECK (true, -1002, "type not supported") ;
            }

            if (MM_storage == MM_skew_symmetric && (type == GrB_BOOL ||
                type == GrB_UINT8  || type == GrB_UINT16 ||
                type == GrB_UINT32 || type == GrB_UINT64))
            {
                // matrices with unsigned types cannot be skew-symmetric
                LG_CHECK (true, -1002, "skew-symmetric matrices cannot have an"
                    " unsigned type") ;
            }

        }
        else if (is_blank_line (buf))
        {

            // -----------------------------------------------------------------
            // blank line or comment line
            // -----------------------------------------------------------------

            continue ;

        }
        else
        {

            // -----------------------------------------------------------------
            // read the first data line
            // -----------------------------------------------------------------

            // format: [nrows ncols nvals] or just [nrows ncols]

            int nitems = sscanf (buf, "%" SCNu64 " %" SCNu64 " %" SCNu64,
                &nrows, &ncols, &nvals) ;

            if (nitems == 2)
            {
                // a dense matrix
                if (!got_mm_header)
                {
                    // if no header, treat it as if it were
                    // %%MatrixMarket matrix array real general
                    MM_fmt = MM_array ;
                    MM_type = MM_real ;
                    MM_storage = MM_general ;
                    type = GrB_FP64 ;
                }
                if (MM_storage == MM_general)
                {
                    // dense general matrix
                    nvals = nrows * ncols ;
                }
                else
                {
                    // dense symmetric, skew-symmetric, or hermitian matrix
                    nvals = nrows + ((nrows * nrows - nrows) / 2) ;
                }
            }
            else if (nitems == 3)
            {
                // a sparse matrix
                if (!got_mm_header)
                {
                    // if no header, treat it as if it were
                    // %%MatrixMarket matrix coordinate real general
                    MM_fmt = MM_coordinate ;
                    MM_type = MM_real ;
                    MM_storage = MM_general ;
                    type = GrB_FP64 ;
                }
            }
            else
            {
                // wrong number of items in first data line
                LG_CHECK (true, -1002, "invalid 1st line") ;
            }

            if (nrows != ncols)
            {
                // a rectangular matrix must be in the general storage
                LG_CHECK (! (MM_storage == MM_general), -1002,
                    "invalid rectangular") ;
            }

            //------------------------------------------------------------------
            // header has been read in
            //------------------------------------------------------------------

            break ;
        }
    }

    //--------------------------------------------------------------------------
    // create the matrix
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Matrix_new (A, type, nrows, ncols)) ;
    *A_type = type;

    //--------------------------------------------------------------------------
    // quick return for empty matrix
    //--------------------------------------------------------------------------

    if (nrows == 0 || ncols == 0 || nvals == 0)
    {
        // success: return an empty matrix.  This is not an error.
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // read the entries
    //--------------------------------------------------------------------------

    GrB_Index i = -1, j = 0 ;
    GrB_Index nvals2 = 0 ;
    for (int64_t k = 0 ; k < nvals ; k++)
    {

        //----------------------------------------------------------------------
        // get the next triplet, skipping blank lines and comment lines
        //----------------------------------------------------------------------

        char x [MAXLINE] ;

        while (true)
        {

            //------------------------------------------------------------------
            // read the file until finding the next triplet
            //------------------------------------------------------------------

            bool ok = get_line (f, buf) ;
            LG_CHECK (!ok, -1002, "premature EOF") ;
            if (is_blank_line (buf))
            {
                // blank line or comment
                continue ;
            }

            //------------------------------------------------------------------
            // get the row and column index
            //------------------------------------------------------------------

            char *p = buf ;
            if (MM_fmt == MM_array)
            {
                // array format, column major order
                i++ ;
                if (i == nrows)
                {
                    j++ ;
                    if (MM_storage == MM_general)
                    {
                        // dense matrix in column major order 
                        i = 0 ;
                    }
                    else
                    {
                        // dense matrix in column major order, only the lower
                        // triangular form is present, including the diagonal
                        i = j ;
                    }
                }
            }
            else
            {
                // coordinate format; read the row and column index
                LG_CHECK (sscanf (p, "%" SCNu64 " %" SCNu64, &i, &j) != 2,
                    -1002, "indices invalid") ;
                // convert from 1-based to 0-based.
                i-- ;
                j-- ;
                // advance p to the 3rd token to get the value of the entry
                while (*p &&  isspace (*p)) p++ ;   // skip any leading spaces
                while (*p && !isspace (*p)) p++ ;   // skip nrows
                while (*p &&  isspace (*p)) p++ ;   // skip any spaces
                while (*p && !isspace (*p)) p++ ;   // skip nrows
            }

            //------------------------------------------------------------------
            // read the value of the entry
            //------------------------------------------------------------------

            while (*p && isspace (*p)) p++ ;        // skip any spaces

            ok = read_entry (p, type, MM_type == MM_pattern, x) ;
            LG_CHECK (!ok, -1002, "entry invalid") ;

            //------------------------------------------------------------------
            // set the value in the matrix
            //------------------------------------------------------------------

            nvals2++ ;
            GrB_TRY (set_value (*A, type, i, j, x)) ;

            //------------------------------------------------------------------
            // also set the A(j,i) entry, if symmetric
            //------------------------------------------------------------------

            if (i != j && MM_storage != MM_general)
            {
                if (MM_storage == MM_symmetric)
                {
                    nvals2++ ;
                    GrB_TRY (set_value (*A, type, j, i, x)) ;
                }
                else if (MM_storage == MM_skew_symmetric)
                {
                    nvals2++ ;
                    negate_scalar (type, x) ;
                    GrB_TRY (set_value (*A, type, j, i, x)) ;
                }
                #if 0
                else if (MM_storage == MM_hermitian)
                {
                    nvals2++ ;
                    double complex *value = (double complex *) x ;
                    (*value) = conj (*value) ;
                    GrB_TRY (set_value (*A, type, j, i, x)) ;
                }
                #endif
            }

            // one more entry has been read in
            break ;
        }
    }

    //--------------------------------------------------------------------------
    // check for duplicates
    //--------------------------------------------------------------------------

    GrB_Index nvals3 = 0 ;
    GrB_TRY (GrB_Matrix_nvals (&nvals3, *A)) ;
    LG_CHECK (nvals2 != nvals3, -1002, "duplicate entries present") ;
    return (0) ;
}

