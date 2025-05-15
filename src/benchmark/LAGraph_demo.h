//------------------------------------------------------------------------------
// LAGraph_demo.h: include file for LAGraph/src/benchmark programs
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2023 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#ifndef LAGRAPH_DEMO_H
#define LAGRAPH_DEMO_H

#include <LAGraph.h>
#include <LG_test.h>
#include "LG_internal.h"

#if defined ( __linux__ )
// for mallopt
#include <malloc.h>
#endif

// set this to 1 to check the results using a slow method
#define LG_CHECK_RESULT 0

#define DEAD_CODE -911
#define CATCH(status)                                                         \
{                                                                             \
    printf ("error: %s line: %d, status: %d\n", __FILE__, __LINE__, status) ; \
    if (msg [0] != '\0') printf ("msg: %s\n", msg) ;                          \
    LG_FREE_ALL ;                                                             \
    return (status) ;                                                         \
}

#undef  LAGRAPH_CATCH
#define LAGRAPH_CATCH(status) CATCH (status)

#undef  GRB_CATCH
#define GRB_CATCH(info) CATCH (info)

#define LAGRAPH_BIN_HEADER 512
#define LEN LAGRAPH_BIN_HEADER

//------------------------------------------------------------------------------
// binwrite: write a matrix to a binary file
//------------------------------------------------------------------------------

#undef  LG_FREE_ALL
#define LG_FREE_ALL                         \
{                                           \
    GrB_free (A) ;                          \
    GrB_free (&hyper_switch) ;              \
    LAGraph_Free ((void **) &Ap, NULL) ;    \
    LAGraph_Free ((void **) &Ab, NULL) ;    \
    LAGraph_Free ((void **) &Ah, NULL) ;    \
    LAGraph_Free ((void **) &Ai, NULL) ;    \
    LAGraph_Free ((void **) &Ax, NULL) ;    \
}

#define FWRITE(p,s,n)                   \
{                                       \
    if (fwrite (p, s, n, f) != n)       \
    {                                   \
        CATCH (LAGRAPH_IO_ERROR) ;      \
    }                                   \
}

static inline int binwrite  // returns 0 if successful, < 0 on error
(
    GrB_Matrix *A,          // matrix to write to the file
    FILE *f,                // file to write it to
    const char *comments    // comments to add to the file, up to 210 characters
                            // in length, not including the terminating null
                            // byte. Ignored if NULL.  Characters past
                            // the 210 limit are silently ignored.
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;
    msg [0] = '\0' ;

#if !LAGRAPH_SUITESPARSE
    printf ("SuiteSparse:GraphBLAS required to write binary *.grb files\n") ;
    return (GrB_NOT_IMPLEMENTED) ;
#else

    GrB_Scalar hyper_switch = NULL ;
    GrB_Index *Ap = NULL, *Ai = NULL, *Ah = NULL ;
    void *Ax = NULL ;
    int8_t *Ab = NULL ;
    if (A == NULL || *A == NULL || f == NULL) CATCH (GrB_NULL_POINTER) ;

    GRB_TRY (GrB_wait (*A, GrB_MATERIALIZE)) ;

    //--------------------------------------------------------------------------
    // determine the basic matrix properties
    //--------------------------------------------------------------------------

    int32_t fmt = -999 ;
    GRB_TRY (GrB_get (*A, &fmt, GrB_STORAGE_ORIENTATION_HINT)) ;

    bool is_hyper = false ;
    bool is_sparse = false ;
    bool is_bitmap = false ;
    bool is_full  = false ;
    int32_t kind ;
    double hyper = -999 ;

    GRB_TRY (GrB_Scalar_new (&hyper_switch, GrB_FP64)) ;
    GRB_TRY (LG_GET_HYPER_SWITCH (*A, hyper_switch)) ;
    GRB_TRY (GrB_Scalar_extractElement (&hyper, hyper_switch)) ;
    GRB_TRY (LG_GET_FORMAT_HINT (*A, &kind)) ;

    switch (kind)
    {
        default :
        case 0 : // for backward compatibility with prior versions
        case 2 : is_sparse = true ; break ;
        case 1 : is_hyper  = true ; break ;
        case 4 : is_bitmap = true ; break ;
        case 8 : is_full   = true ; break ;
    }

    //--------------------------------------------------------------------------
    // export the matrix
    //--------------------------------------------------------------------------

    GrB_Type type ;
    GrB_Index nrows, ncols, nvals, nvec ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, *A)) ;
    size_t typesize ;
    int64_t nonempty = -1 ;
    char *fmt_string ;
    bool jumbled, iso ;
    GrB_Index Ap_size, Ah_size, Ab_size, Ai_size, Ax_size ;

    if (fmt == GrB_COLMAJOR && is_hyper)
    {
        // hypersparse CSC
        GRB_TRY (GxB_Matrix_export_HyperCSC (A, &type, &nrows, &ncols,
            &Ap, &Ah, &Ai, &Ax, &Ap_size, &Ah_size, &Ai_size, &Ax_size,
            &iso, &nvec, &jumbled, NULL)) ;
        fmt_string = "HCSC" ;
    }
    else if (fmt == GrB_ROWMAJOR && is_hyper)
    {
        // hypersparse CSR
        GRB_TRY (GxB_Matrix_export_HyperCSR (A, &type, &nrows, &ncols,
            &Ap, &Ah, &Ai, &Ax, &Ap_size, &Ah_size, &Ai_size, &Ax_size,
            &iso, &nvec, &jumbled, NULL)) ;
        fmt_string = "HCSR" ;
    }
    else if (fmt == GrB_COLMAJOR && is_sparse)
    {
        // standard CSC
        GRB_TRY (GxB_Matrix_export_CSC (A, &type, &nrows, &ncols,
            &Ap, &Ai, &Ax, &Ap_size, &Ai_size, &Ax_size,
            &iso, &jumbled, NULL)) ;
        nvec = ncols ;
        fmt_string = "CSC " ;
    }
    else if (fmt == GrB_ROWMAJOR && is_sparse)
    {
        // standard CSR
        GRB_TRY (GxB_Matrix_export_CSR (A, &type, &nrows, &ncols,
            &Ap, &Ai, &Ax, &Ap_size, &Ai_size, &Ax_size,
            &iso, &jumbled, NULL)) ;
        nvec = nrows ;
        fmt_string = "CSR " ;
    }
    else if (fmt == GrB_COLMAJOR && is_bitmap)
    {
        // bitmap by col
        GRB_TRY (GxB_Matrix_export_BitmapC (A, &type, &nrows, &ncols,
            &Ab, &Ax, &Ab_size, &Ax_size,
            &iso, &nvals, NULL)) ;
        nvec = ncols ;
        fmt_string = "BITMAPC" ;
    }
    else if (fmt == GrB_ROWMAJOR && is_bitmap)
    {
        // bitmap by row
        GRB_TRY (GxB_Matrix_export_BitmapR (A, &type, &nrows, &ncols,
            &Ab, &Ax, &Ab_size, &Ax_size,
            &iso, &nvals, NULL)) ;
        nvec = nrows ;
        fmt_string = "BITMAPR" ;
    }
    else if (fmt == GrB_COLMAJOR && is_full)
    {
        // full by col
        GRB_TRY (GxB_Matrix_export_FullC (A, &type, &nrows, &ncols,
            &Ax, &Ax_size,
            &iso, NULL)) ;
        nvec = ncols ;
        fmt_string = "FULLC" ;
    }
    else if (fmt == GrB_ROWMAJOR && is_full)
    {
        // full by row
        GRB_TRY (GxB_Matrix_export_FullR (A, &type, &nrows, &ncols,
            &Ax, &Ax_size,
            &iso, NULL)) ;
        nvec = nrows ;
        fmt_string = "FULLC" ;
    }
    else
    {
        CATCH (DEAD_CODE) ;    // this "cannot" happen
    }

    //--------------------------------------------------------------------------
    // create the type string
    //--------------------------------------------------------------------------

    LAGRAPH_TRY (LAGraph_SizeOfType (&typesize, type, msg)) ;

    char typename [LEN] ;
    int32_t typecode ;
    if      (type == GrB_BOOL  )
    {
        snprintf (typename, LEN, "GrB_BOOL  ") ;
        typecode = 0 ;
    }
    else if (type == GrB_INT8  )
    {
        snprintf (typename, LEN, "GrB_INT8  ") ;
        typecode = 1 ;
    }
    else if (type == GrB_INT16 )
    {
        snprintf (typename, LEN, "GrB_INT16 ") ;
        typecode = 2 ;
    }
    else if (type == GrB_INT32 )
    {
        snprintf (typename, LEN, "GrB_INT32 ") ;
        typecode = 3 ;
    }
    else if (type == GrB_INT64 )
    {
        snprintf (typename, LEN, "GrB_INT64 ") ;
        typecode = 4 ;
    }
    else if (type == GrB_UINT8 )
    {
        snprintf (typename, LEN, "GrB_UINT8 ") ;
        typecode = 5 ;
    }
    else if (type == GrB_UINT16)
    {
        snprintf (typename, LEN, "GrB_UINT16") ;
        typecode = 6 ;
    }
    else if (type == GrB_UINT32)
    {
        snprintf (typename, LEN, "GrB_UINT32") ;
        typecode = 7 ;
    }
    else if (type == GrB_UINT64)
    {
        snprintf (typename, LEN, "GrB_UINT64") ;
        typecode = 8 ;
    }
    else if (type == GrB_FP32  )
    {
        snprintf (typename, LEN, "GrB_FP32  ") ;
        typecode = 9 ;
    }
    else if (type == GrB_FP64  )
    {
        snprintf (typename, LEN, "GrB_FP64  ") ;
        typecode = 10 ;
    }
    else
    {
        // unsupported type
        CATCH (GrB_NOT_IMPLEMENTED) ;
    }
    typename [72] = '\0' ;

    //--------------------------------------------------------------------------
    // write the header in ascii
    //--------------------------------------------------------------------------

    // The header is informational only, for "head" command, so the file can
    // be visually inspected.

    int32_t ver [3] ;
    GRB_TRY (GrB_get (GrB_GLOBAL, &(ver [0]), GrB_LIBRARY_VER_MAJOR)) ;
    GRB_TRY (GrB_get (GrB_GLOBAL, &(ver [1]), GrB_LIBRARY_VER_MINOR)) ;
    GRB_TRY (GrB_get (GrB_GLOBAL, &(ver [2]), GrB_LIBRARY_VER_PATCH)) ;

    char version [LEN] ;
    snprintf (version, LEN, "%d.%d.%d (LAGraph)", ver [0], ver [1], ver [2]) ;
    version [25] = '\0' ;

    char user [LEN] ;
    for (int k = 0 ; k < LEN ; k++) user [k] = ' ' ;
    user [0] = '\n' ;
    if (comments != NULL)
    {
        strncpy (user, comments, 210) ;
    }
    user [210] = '\0' ;

    char header [LAGRAPH_BIN_HEADER] ;
    int32_t len = snprintf (header, LAGRAPH_BIN_HEADER,
        "SuiteSparse:GraphBLAS matrix\nv%-25s\n"
        "nrows:  %-18" PRIu64 "\n"
        "ncols:  %-18" PRIu64 "\n"
        "nvec:   %-18" PRIu64 "\n"
        "nvals:  %-18" PRIu64 "\n"
        "format: %-8s\n"
        "size:   %-18" PRIu64 "\n"
        "type:   %-72s\n"
        "iso:    %1d\n"
        "%-210s\n\n",
        version, nrows, ncols, nvec, nvals, fmt_string, (uint64_t) typesize,
        typename, iso, user) ;

    // printf ("header len %d\n", len) ;
    for (int32_t k = len ; k < LAGRAPH_BIN_HEADER ; k++) header [k] = ' ' ;
    header [LAGRAPH_BIN_HEADER-1] = '\0' ;
    FWRITE (header, sizeof (char), LAGRAPH_BIN_HEADER) ;

    //--------------------------------------------------------------------------
    // write the scalar content
    //--------------------------------------------------------------------------

    if (iso)
    {
        // kind is 1, 2, 4, or 8: add 100 if the matrix is iso
        kind = kind + 100 ;
    }

    FWRITE (&fmt,      sizeof (int32_t), 1) ;
    FWRITE (&kind,     sizeof (int32_t), 1) ;
    FWRITE (&hyper,    sizeof (double), 1) ;
    FWRITE (&nrows,    sizeof (GrB_Index), 1) ;
    FWRITE (&ncols,    sizeof (GrB_Index), 1) ;
    FWRITE (&nonempty, sizeof (int64_t), 1) ;
    FWRITE (&nvec,     sizeof (GrB_Index), 1) ;
    FWRITE (&nvals,    sizeof (GrB_Index), 1) ;
    FWRITE (&typecode, sizeof (int32_t), 1) ;
    FWRITE (&typesize, sizeof (size_t), 1) ;

    //--------------------------------------------------------------------------
    // write the array content
    //--------------------------------------------------------------------------

    if (is_hyper)
    {
        FWRITE (Ap, sizeof (GrB_Index), nvec+1) ;
        FWRITE (Ah, sizeof (GrB_Index), nvec) ;
        FWRITE (Ai, sizeof (GrB_Index), nvals) ;
        FWRITE (Ax, typesize, (iso ? 1 : nvals)) ;
    }
    else if (is_sparse)
    {
        FWRITE (Ap, sizeof (GrB_Index), nvec+1) ;
        FWRITE (Ai, sizeof (GrB_Index), nvals) ;
        FWRITE (Ax, typesize, (iso ? 1 : nvals)) ;
    }
    else if (is_bitmap)
    {
        FWRITE (Ab, sizeof (int8_t), nrows*ncols) ;
        FWRITE (Ax, typesize, (iso ? 1 : (nrows*ncols))) ;
    }
    else
    {
        FWRITE (Ax, typesize, (iso ? 1 : (nrows*ncols))) ;
    }

    //--------------------------------------------------------------------------
    // re-import the matrix
    //--------------------------------------------------------------------------

    if (fmt == GrB_COLMAJOR && is_hyper)
    {
        // hypersparse CSC
        GRB_TRY (GxB_Matrix_import_HyperCSC (A, type, nrows, ncols,
            &Ap, &Ah, &Ai, &Ax, Ap_size, Ah_size, Ai_size, Ax_size,
            iso, nvec, jumbled, NULL)) ;
    }
    else if (fmt == GrB_ROWMAJOR && is_hyper)
    {
        // hypersparse CSR
        GRB_TRY (GxB_Matrix_import_HyperCSR (A, type, nrows, ncols,
            &Ap, &Ah, &Ai, &Ax, Ap_size, Ah_size, Ai_size, Ax_size,
            iso, nvec, jumbled, NULL)) ;
    }
    else if (fmt == GrB_COLMAJOR && is_sparse)
    {
        // standard CSC
        GRB_TRY (GxB_Matrix_import_CSC (A, type, nrows, ncols,
            &Ap, &Ai, &Ax, Ap_size, Ai_size, Ax_size,
            iso, jumbled, NULL)) ;
    }
    else if (fmt == GrB_ROWMAJOR && is_sparse)
    {
        // standard CSR
        GRB_TRY (GxB_Matrix_import_CSR (A, type, nrows, ncols,
            &Ap, &Ai, &Ax, Ap_size, Ai_size, Ax_size,
            iso, jumbled, NULL)) ;
    }
    else if (fmt == GrB_COLMAJOR && is_bitmap)
    {
        // bitmap by col
        GRB_TRY (GxB_Matrix_import_BitmapC (A, type, nrows, ncols,
            &Ab, &Ax, Ab_size, Ax_size,
            iso, nvals, NULL)) ;
    }
    else if (fmt == GrB_ROWMAJOR && is_bitmap)
    {
        // bitmap by row
        GRB_TRY (GxB_Matrix_import_BitmapR (A, type, nrows, ncols,
            &Ab, &Ax, Ab_size, Ax_size,
            iso, nvals, NULL)) ;
    }
    else if (fmt == GrB_COLMAJOR && is_full)
    {
        // full by col
        GRB_TRY (GxB_Matrix_import_FullC (A, type, nrows, ncols,
            &Ax, Ax_size,
            iso, NULL)) ;
    }
    else if (fmt == GrB_ROWMAJOR && is_full)
    {
        // full by row
        GRB_TRY (GxB_Matrix_import_FullR (A, type, nrows, ncols,
            &Ax, Ax_size,
            iso, NULL)) ;
    }
    else
    {
        CATCH (DEAD_CODE) ;    // this "cannot" happen
    }

    GRB_TRY (GrB_Scalar_setElement (hyper_switch, hyper)) ;
    GRB_TRY (LG_SET_HYPER_SWITCH (*A, hyper_switch)) ;
    GrB_free (&hyper_switch) ;
    return (GrB_SUCCESS) ;
#endif
}

//------------------------------------------------------------------------------
// binread: read a matrix from a binary file
//------------------------------------------------------------------------------

#define FREAD(p,s,n)                    \
{                                       \
    if (fread (p, s, n, f) != n)        \
    {                                   \
        CATCH (-1001) ; /* file I/O error */ \
    }                                   \
}

static inline int binread   // returns 0 if successful, -1 if failure
(
    GrB_Matrix *A,          // matrix to read from the file
    FILE *f                 // file to read it from, already open
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;
    msg [0] = '\0' ;

#if !LAGRAPH_SUITESPARSE
    printf ("SuiteSparse:GraphBLAS required to read binary *.grb files\n") ;
    return (GrB_NOT_IMPLEMENTED) ;
#else

    GrB_Scalar hyper_switch = NULL ;
    GrB_Index *Ap = NULL, *Ai = NULL, *Ah = NULL ;
    int8_t *Ab = NULL ;
    void *Ax = NULL ;
    if (A == NULL || f == NULL) CATCH (GrB_NULL_POINTER) ;
    (*A) = NULL ;

    //--------------------------------------------------------------------------
    // basic matrix properties
    //--------------------------------------------------------------------------

    int32_t fmt = -999 ;
    bool is_hyper, is_sparse, is_bitmap, is_full ;
    int32_t kind, typecode ;
    double hyper = -999 ;
    GrB_Type type ;
    GrB_Index nrows, ncols, nvals, nvec ;
    size_t typesize ;
    int64_t nonempty ;

    //--------------------------------------------------------------------------
    // read the header (and ignore it)
    //--------------------------------------------------------------------------

    // The header is informational only, for "head" command, so the file can
    // be visually inspected.

    char header [LAGRAPH_BIN_HEADER] ;
    FREAD (header, sizeof (char), LAGRAPH_BIN_HEADER) ;
    // printf ("%s\n", header) ;

    //--------------------------------------------------------------------------
    // read the scalar content
    //--------------------------------------------------------------------------

    FREAD (&fmt,      sizeof (int32_t), 1) ;
    FREAD (&kind,     sizeof (int32_t), 1) ;
    FREAD (&hyper,    sizeof (double), 1) ;
    FREAD (&nrows,    sizeof (GrB_Index), 1) ;
    FREAD (&ncols,    sizeof (GrB_Index), 1) ;
    FREAD (&nonempty, sizeof (int64_t), 1) ;
    FREAD (&nvec,     sizeof (GrB_Index), 1) ;
    FREAD (&nvals,    sizeof (GrB_Index), 1) ;
    FREAD (&typecode, sizeof (int32_t), 1) ;
    FREAD (&typesize, sizeof (size_t), 1) ;

    bool iso = false ;
    if (kind > 100)
    {
        iso = true ;
        kind = kind - 100 ;
    }

    is_hyper  = (kind == 1) ;
    is_sparse = (kind == 0 || kind == GxB_SPARSE) ;
    is_bitmap = (kind == GxB_BITMAP) ;
    is_full   = (kind == GxB_FULL) ;

    switch (typecode)
    {
        case 0:  type = GrB_BOOL        ; break ;
        case 1:  type = GrB_INT8        ; break ;
        case 2:  type = GrB_INT16       ; break ;
        case 3:  type = GrB_INT32       ; break ;
        case 4:  type = GrB_INT64       ; break ;
        case 5:  type = GrB_UINT8       ; break ;
        case 6:  type = GrB_UINT16      ; break ;
        case 7:  type = GrB_UINT32      ; break ;
        case 8:  type = GrB_UINT64      ; break ;
        case 9:  type = GrB_FP32        ; break ;
        case 10: type = GrB_FP64        ; break ;
//      case 11: type = GxB_FC32        ; break ;
//      case 12: type = GxB_FC64        ; break ;
        default: CATCH (GrB_NOT_IMPLEMENTED) ;    // unknown or unsupported type
    }

    //--------------------------------------------------------------------------
    // allocate the array content
    //--------------------------------------------------------------------------

    GrB_Index Ap_len = 0, Ap_size = 0 ;
    GrB_Index Ah_len = 0, Ah_size = 0 ;
    GrB_Index Ab_len = 0, Ab_size = 0 ;
    GrB_Index Ai_len = 0, Ai_size = 0 ;
    GrB_Index Ax_len = 0, Ax_size = 0 ;

    bool ok = true ;
    if (is_hyper)
    {
        Ap_len = nvec+1 ;
        Ah_len = nvec ;
        Ai_len = nvals ;
        Ax_len = nvals ;
        LAGraph_Malloc ((void **) &Ap, Ap_len, sizeof (GrB_Index), msg) ;
        LAGraph_Malloc ((void **) &Ah, Ah_len, sizeof (GrB_Index), msg) ;
        LAGraph_Malloc ((void **) &Ai, Ai_len, sizeof (GrB_Index), msg) ;
        Ap_size = Ap_len * sizeof (GrB_Index) ;
        Ah_size = Ah_len * sizeof (GrB_Index) ;
        Ai_size = Ai_len * sizeof (GrB_Index) ;
        ok = (Ap != NULL && Ah != NULL && Ai != NULL) ;
    }
    else if (is_sparse)
    {
        Ap_len = nvec+1 ;
        Ai_len = nvals ;
        Ax_len = nvals ;
        LAGraph_Malloc ((void **) &Ap, Ap_len, sizeof (GrB_Index), msg) ;
        LAGraph_Malloc ((void **) &Ai, Ai_len, sizeof (GrB_Index), msg) ;
        Ap_size = Ap_len * sizeof (GrB_Index) ;
        Ai_size = Ai_len * sizeof (GrB_Index) ;
        ok = (Ap != NULL && Ai != NULL) ;
    }
    else if (is_bitmap)
    {
        Ab_len = nrows*ncols ;
        Ax_len = nrows*ncols ;
        LAGraph_Malloc ((void **) &Ab, nrows*ncols, sizeof (int8_t), msg) ;
        Ab_size = Ab_len * sizeof (GrB_Index) ;
        ok = (Ab != NULL) ;
    }
    else if (is_full)
    {
        Ax_len = nrows*ncols ;
    }
    else
    {
        CATCH (DEAD_CODE) ;    // this "cannot" happen
    }
    LAGraph_Malloc ((void **) &Ax, iso ? 1 : Ax_len, typesize, msg) ;
    Ax_size = (iso ? 1 : Ax_len) * typesize ;
    ok = ok && (Ax != NULL) ;
    if (!ok) CATCH (GrB_OUT_OF_MEMORY) ;        // out of memory

    //--------------------------------------------------------------------------
    // read the array content
    //--------------------------------------------------------------------------

    if (is_hyper)
    {
        FREAD (Ap, sizeof (GrB_Index), Ap_len) ;
        FREAD (Ah, sizeof (GrB_Index), Ah_len) ;
        FREAD (Ai, sizeof (GrB_Index), Ai_len) ;
    }
    else if (is_sparse)
    {
        FREAD (Ap, sizeof (GrB_Index), Ap_len) ;
        FREAD (Ai, sizeof (GrB_Index), Ai_len) ;
    }
    else if (is_bitmap)
    {
        FREAD (Ab, sizeof (int8_t), Ab_len) ;
    }

    FREAD (Ax, typesize, (iso ? 1 : Ax_len)) ;

    //--------------------------------------------------------------------------
    // import the matrix
    //--------------------------------------------------------------------------

    if (fmt == GrB_COLMAJOR && is_hyper)
    {
        // hypersparse CSC
        GRB_TRY (GxB_Matrix_import_HyperCSC (A, type, nrows, ncols,
            &Ap, &Ah, &Ai, &Ax, Ap_size, Ah_size, Ai_size, Ax_size,
            iso, nvec, false, NULL)) ;
    }
    else if (fmt == GrB_ROWMAJOR && is_hyper)
    {
        // hypersparse CSR
        GRB_TRY (GxB_Matrix_import_HyperCSR (A, type, nrows, ncols,
            &Ap, &Ah, &Ai, &Ax, Ap_size, Ah_size, Ai_size, Ax_size,
            iso, nvec, false, NULL)) ;
    }
    else if (fmt == GrB_COLMAJOR && is_sparse)
    {
        // standard CSC
        GRB_TRY (GxB_Matrix_import_CSC (A, type, nrows, ncols,
            &Ap, &Ai, &Ax, Ap_size, Ai_size, Ax_size,
            iso, false, NULL)) ;
    }
    else if (fmt == GrB_ROWMAJOR && is_sparse)
    {
        // standard CSR
        GRB_TRY (GxB_Matrix_import_CSR (A, type, nrows, ncols,
            &Ap, &Ai, &Ax, Ap_size, Ai_size, Ax_size,
            iso, false, NULL)) ;
    }
    else if (fmt == GrB_COLMAJOR && is_bitmap)
    {
        // bitmap by col
        GRB_TRY (GxB_Matrix_import_BitmapC (A, type, nrows, ncols,
            &Ab, &Ax, Ab_size, Ax_size,
            iso, nvals, NULL)) ;
    }
    else if (fmt == GrB_ROWMAJOR && is_bitmap)
    {
        // bitmap by row
        GRB_TRY (GxB_Matrix_import_BitmapR (A, type, nrows, ncols,
            &Ab, &Ax, Ab_size, Ax_size,
            iso, nvals, NULL)) ;
    }
    else if (fmt == GrB_COLMAJOR && is_full)
    {
        // full by col
        GRB_TRY (GxB_Matrix_import_FullC (A, type, nrows, ncols,
            &Ax, Ax_size,
            iso, NULL)) ;
    }
    else if (fmt == GrB_ROWMAJOR && is_full)
    {
        // full by row
        GRB_TRY (GxB_Matrix_import_FullR (A, type, nrows, ncols,
            &Ax, Ax_size,
            iso, NULL)) ;
    }
    else
    {
        CATCH (DEAD_CODE) ;    // this "cannot" happen
    }

    GRB_TRY (GrB_Scalar_new (&hyper_switch, GrB_FP64)) ;
    GRB_TRY (GrB_Scalar_setElement (hyper_switch, hyper)) ;
    GRB_TRY (LG_SET_HYPER_SWITCH (*A, hyper_switch)) ;
    GrB_free (&hyper_switch) ;
    return (GrB_SUCCESS) ;
#endif
}

//------------------------------------------------------------------------------
// readproblem: read a GAP problem from a file
//------------------------------------------------------------------------------

#undef  LG_FREE_WORK
#define LG_FREE_WORK                \
{                                   \
    GrB_free (&A) ;                 \
    GrB_free (&A2) ;                \
    GrB_free (&M) ;                 \
    if (f != NULL) fclose (f) ;     \
    f = NULL ;                      \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    LAGraph_Delete (G, NULL) ;      \
    GrB_free (src_nodes) ;          \
}

// usage:
// test_whatever < matrixfile.mtx
// test_whatever matrixfile.mtx sourcenodes.mtx
// The matrixfile may also have a grb suffix.

static int readproblem          // returns 0 if successful, -1 if failure
(
    // output
    LAGraph_Graph *G,           // graph from the file
    GrB_Matrix *src_nodes,      // source nodes
    // inputs
    bool make_symmetric,        // if true, always return G as undirected
    bool remove_self_edges,     // if true, remove self edges
    bool structural,            // if true, return G->A as bool (all true)
    GrB_Type pref,              // if non-NULL, typecast G->A to this type
    bool ensure_positive,       // if true, ensure all entries are > 0
    int argc,                   // input to main test program
    char **argv                 // input to main test program
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;
    msg [0] = '\0' ;
    GrB_Matrix A = NULL, A2 = NULL, M = NULL ;
    GrB_Type atype = NULL ;
    FILE *f = NULL ;
    if (G == NULL) CATCH (GrB_NULL_POINTER) ;
    (*G) = NULL ;
    if (src_nodes != NULL) (*src_nodes) = NULL ;
    GrB_Type src_type = NULL;

    //--------------------------------------------------------------------------
    // read in a matrix from a file
    //--------------------------------------------------------------------------

    double t_read = LAGraph_WallClockTime ( ) ;

    if (argc > 1)
    {
        // Usage:
        //      ./test_whatever matrixfile.mtx [sources.mtx]
        //      ./test_whatever matrixfile.grb [sources.mtx]

        // read in the file in Matrix Market format from the input file
        char *filename = argv [1] ;
        printf ("matrix: %s\n", filename) ;

        // find the filename extension
        size_t len = strlen (filename) ;
        char *ext = NULL ;
        for (int k = len-1 ; k >= 0 ; k--)
        {
            if (filename [k] == '.')
            {
                ext = filename + k ;
                printf ("[%s]\n", ext) ;
                break ;
            }
        }

        bool is_binary = (ext != NULL && strncmp (ext, ".grb", 4) == 0) ;

        if (is_binary)
        {
            printf ("Reading binary file: %s\n", filename) ;
            f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Binary file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            if (binread (&A, f) < 0) CATCH (-1001) ;    // file I/O error
            fclose (f) ;
            f = NULL ;
        }
        else
        {
            printf ("Reading matrix market file: %s\n", filename) ;
            f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Matrix market file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            int result = LAGraph_MMRead (&A, f, msg) ;
            if (result != GrB_SUCCESS)
            {
                printf ("LAGraph_MMRead failed to read matrix: %s\n",
                    filename) ;
                printf ("result: %d msg: %s\n", result, msg) ;
            }
            LAGRAPH_TRY (result) ;
            fclose (f) ;
            f = NULL ;
        }

        // read in source nodes in Matrix Market format from the input file
        if (argc > 2 && src_nodes != NULL)
        {
            // do not read in the file if the name starts with "-"
            filename = argv [2] ;
            if (filename [0] != '-')
            {
                printf ("sources: %s\n", filename) ;
                f = fopen (filename, "r") ;
                if (f == NULL)
                {
                    printf ("Source node file not found: [%s]\n", filename) ;
                    exit (1) ;
                }
                int result = LAGraph_MMRead (src_nodes, f, msg) ;
                if (result != GrB_SUCCESS)
                {
                    printf ("LAGraph_MMRead failed to read source nodes"
                        " from: %s\n", filename) ;
                    printf ("result: %d msg: %s\n", result, msg) ;
                }
                LAGRAPH_TRY (result) ;
                fclose (f) ;
                f = NULL ;
            }
        }
    }
    else
    {

        // Usage:  ./test_whatever < matrixfile.mtx
        printf ("matrix: from stdin\n") ;

        // read in the file in Matrix Market format from stdin
        int result = LAGraph_MMRead (&A, stdin, msg) ;
        if (result != GrB_SUCCESS)
        {
            printf ("LAGraph_MMRead failed to read: stdin\n") ;
            printf ("result: %d msg: %s\n", result, msg) ;
        }
        LAGRAPH_TRY (result) ;
    }

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols ;
    GRB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_Index n = nrows ;
    if (nrows != ncols) CATCH (GrB_DIMENSION_MISMATCH) ;    // A must be square

    //--------------------------------------------------------------------------
    // convert to 32-bit
    //--------------------------------------------------------------------------

    #if LAGRAPH_SUITESPARSE
    #if GxB_IMPLEMENTATION >= GxB_VERSION (10,0,0)
    GRB_TRY (GrB_Matrix_set_INT32 (A, 32, GxB_ROWINDEX_INTEGER_HINT)) ;
    GRB_TRY (GrB_Matrix_set_INT32 (A, 32, GxB_COLINDEX_INTEGER_HINT)) ;
    GRB_TRY (GrB_Matrix_set_INT32 (A, 32, GxB_OFFSET_INTEGER_HINT)) ;
    printf ("A converted to 32-bit\n") ;
    #endif
    #endif

    //--------------------------------------------------------------------------
    // typecast, if requested
    //--------------------------------------------------------------------------

    #if LAGRAPH_SUITESPARSE
    GRB_TRY (GxB_Matrix_type (&atype, A)) ;
    #else
    char aname [1024] ;
    LG_TRY (LAGraph_Matrix_TypeName (aname, A, msg)) ;
    LG_TRY (LAGraph_TypeFromName (&atype, aname, msg)) ;
    #endif

    if (structural)
    {
        // convert to boolean, with all entries true
        atype = GrB_BOOL ;
        LAGRAPH_TRY (LAGraph_Matrix_Structure (&A2, A, msg)) ;
    }
    else if (pref != NULL && atype != pref)
    {
        // convert to the requested type
        GRB_TRY (GrB_Matrix_new (&A2, pref, n, n)) ;
        atype = pref ;

        GrB_UnaryOp op = NULL ;
        if      (pref == GrB_BOOL  ) op = GrB_IDENTITY_BOOL ;
        else if (pref == GrB_INT8  ) op = GrB_IDENTITY_INT8 ;
        else if (pref == GrB_INT16 ) op = GrB_IDENTITY_INT16 ;
        else if (pref == GrB_INT32 ) op = GrB_IDENTITY_INT32 ;
        else if (pref == GrB_INT64 ) op = GrB_IDENTITY_INT64 ;
        else if (pref == GrB_UINT8 ) op = GrB_IDENTITY_UINT8 ;
        else if (pref == GrB_UINT16) op = GrB_IDENTITY_UINT16 ;
        else if (pref == GrB_UINT32) op = GrB_IDENTITY_UINT32 ;
        else if (pref == GrB_UINT64) op = GrB_IDENTITY_UINT64 ;
        else if (pref == GrB_FP32  ) op = GrB_IDENTITY_FP32 ;
        else if (pref == GrB_FP64  ) op = GrB_IDENTITY_FP64 ;
//      else if (pref == GxB_FC32  ) op = GxB_IDENTITY_FC32 ;
//      else if (pref == GxB_FC64  ) op = GxB_IDENTITY_FC64 ;
        else CATCH (GrB_NOT_IMPLEMENTED) ;    // unsupported type

        GRB_TRY (GrB_apply (A2, NULL, NULL, op, A, NULL)) ;
    }

    if (A2 != NULL)
    {
        GrB_free (&A) ;
        A = A2 ;
        A2 = NULL ;
        GRB_TRY (GrB_wait (A, GrB_MATERIALIZE)) ;
    }

    //--------------------------------------------------------------------------
    // construct the initial graph
    //--------------------------------------------------------------------------

    bool A_is_symmetric =
        (n == 134217726 ||  // HACK for kron
         n == 134217728) ;  // HACK for urand

    LAGraph_Kind G_kind = A_is_symmetric ?  LAGraph_ADJACENCY_UNDIRECTED :
        LAGraph_ADJACENCY_DIRECTED ;
    LAGRAPH_TRY (LAGraph_New (G, &A, G_kind, msg)) ;
    // LAGRAPH_TRY (LAGraph_Graph_Print (*G, 2, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // remove self-edges, if requested
    //--------------------------------------------------------------------------

    if (remove_self_edges)
    {
        LAGRAPH_TRY (LAGraph_DeleteSelfEdges (*G, msg)) ;
    }
    // LAGRAPH_TRY (LAGraph_Graph_Print (*G, 2, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // ensure all entries are > 0, if requested
    //--------------------------------------------------------------------------

    if (!structural && ensure_positive)
    {
        // drop explicit zeros (FUTURE: make this a utility function)
        GrB_IndexUnaryOp idxop = NULL ;
        if      (atype == GrB_BOOL  ) idxop = GrB_VALUENE_BOOL ;
        else if (atype == GrB_INT8  ) idxop = GrB_VALUENE_INT8 ;
        else if (atype == GrB_INT16 ) idxop = GrB_VALUENE_INT16 ;
        else if (atype == GrB_INT32 ) idxop = GrB_VALUENE_INT32 ;
        else if (atype == GrB_INT64 ) idxop = GrB_VALUENE_INT64 ;
        else if (atype == GrB_UINT8 ) idxop = GrB_VALUENE_UINT8 ;
        else if (atype == GrB_UINT16) idxop = GrB_VALUENE_UINT16 ;
        else if (atype == GrB_UINT32) idxop = GrB_VALUENE_UINT32 ;
        else if (atype == GrB_UINT64) idxop = GrB_VALUENE_UINT64 ;
        else if (atype == GrB_FP32  ) idxop = GrB_VALUENE_FP32 ;
        else if (atype == GrB_FP64  ) idxop = GrB_VALUENE_FP64 ;
//      else if (atype == GxB_FC32  ) idxop = GxB_VALUENE_FC32 ;
//      else if (atype == GxB_FC64  ) idxop = GxB_VALUENE_FC64 ;
        if (idxop != NULL)
        {
            GRB_TRY (GrB_select ((*G)->A, NULL, NULL, idxop, (*G)->A, 0, NULL));
        }

        // A = abs (A)
        GrB_UnaryOp op = NULL ;
        if      (atype == GrB_INT8  ) op = GrB_ABS_INT8 ;
        else if (atype == GrB_INT16 ) op = GrB_ABS_INT16 ;
        else if (atype == GrB_INT32 ) op = GrB_ABS_INT32 ;
        else if (atype == GrB_INT64 ) op = GrB_ABS_INT64 ;
        else if (atype == GrB_FP32  ) op = GrB_ABS_FP32 ;
        else if (atype == GrB_FP64  ) op = GrB_ABS_FP64 ;
//      else if (atype == GxB_FC32  ) op = GxB_ABS_FC32 ;
//      else if (atype == GxB_FC64  ) op = GxB_ABS_FC64 ;
        if (op != NULL)
        {
            GRB_TRY (GrB_apply ((*G)->A, NULL, NULL, op, (*G)->A, NULL)) ;
        }
    }

    //--------------------------------------------------------------------------
    // determine the graph properies
    //--------------------------------------------------------------------------

    // LAGRAPH_TRY (LAGraph_Graph_Print (*G, 2, stdout, msg)) ;

    if (!A_is_symmetric)
    {
        // compute G->AT and determine if A has a symmetric structure
        LAGRAPH_TRY (LAGraph_Cached_IsSymmetricStructure (*G, msg)) ;
        if (((*G)->is_symmetric_structure == LAGraph_TRUE) && structural)
        {
            // if G->A has a symmetric structure, declare the graph undirected
            // and free G->AT since it isn't needed.
            (*G)->kind = LAGraph_ADJACENCY_UNDIRECTED ;
            GRB_TRY (GrB_Matrix_free (&((*G)->AT))) ;
        }
        else if (make_symmetric)
        {
            // make sure G->A is symmetric
            bool sym ;
            LAGRAPH_TRY (LAGraph_Matrix_IsEqual (&sym, (*G)->A, (*G)->AT, msg));
            if (!sym)
            {
                printf ("forcing G-> to be symmetric (via A = A+A')\n") ;
                GrB_BinaryOp op = NULL ;
                GrB_Type type ;
                if      (atype == GrB_BOOL  ) op = GrB_LOR ;
                else if (atype == GrB_INT8  ) op = GrB_PLUS_INT8 ;
                else if (atype == GrB_INT16 ) op = GrB_PLUS_INT16 ;
                else if (atype == GrB_INT32 ) op = GrB_PLUS_INT32 ;
                else if (atype == GrB_INT64 ) op = GrB_PLUS_INT64 ;
                else if (atype == GrB_UINT8 ) op = GrB_PLUS_UINT8 ;
                else if (atype == GrB_UINT16) op = GrB_PLUS_UINT16 ;
                else if (atype == GrB_UINT32) op = GrB_PLUS_UINT32 ;
                else if (atype == GrB_UINT64) op = GrB_PLUS_UINT64 ;
                else if (atype == GrB_FP32  ) op = GrB_PLUS_FP32 ;
                else if (atype == GrB_FP64  ) op = GrB_PLUS_FP64 ;
//              else if (type == GxB_FC32  ) op = GxB_PLUS_FC32 ;
//              else if (type == GxB_FC64  ) op = GxB_PLUS_FC64 ;
                else CATCH (GrB_NOT_IMPLEMENTED) ;    // unknown type
                GRB_TRY (GrB_eWiseAdd ((*G)->A, NULL, NULL, op,
                                       (*G)->A, (*G)->AT, NULL)) ;
            }
            // G->AT is not required
            GRB_TRY (GrB_Matrix_free (&((*G)->AT))) ;
            (*G)->kind = LAGraph_ADJACENCY_UNDIRECTED ;
            (*G)->is_symmetric_structure = LAGraph_TRUE ;
        }
    }
    // LAGRAPH_TRY (LAGraph_Graph_Print (*G, 2, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // generate 64 random source nodes, if requested but not provided on input
    //--------------------------------------------------------------------------

    #define NSOURCES 64

    if (src_nodes != NULL && (*src_nodes == NULL))
    {
        src_type = GrB_UINT64;
        GRB_TRY (GrB_Matrix_new (src_nodes, src_type, NSOURCES, 1)) ;
        srand (1) ;
        for (int k = 0 ; k < NSOURCES ; k++)
        {
            uint64_t i = 1 + (rand ( ) % n) ;    // in range 1 to n
            // src_nodes [k] = i
            GRB_TRY (GrB_Matrix_setElement (*src_nodes, i, k, 0)) ;
        }
    }

    if (src_nodes != NULL)
    {
        GRB_TRY (GrB_wait (*src_nodes, GrB_MATERIALIZE)) ;
    }

    //--------------------------------------------------------------------------
    // free workspace, print a summary of the graph, and return result
    //--------------------------------------------------------------------------

    t_read = LAGraph_WallClockTime ( ) - t_read ;
    printf ("read time: %g\n", t_read) ;

    LG_FREE_WORK ;
    // LAGRAPH_TRY (LAGraph_Graph_Print (*G, LAGraph_SHORT, stdout, msg)) ;
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// demo_init: initialize LAGraph for a demo
//------------------------------------------------------------------------------

#undef  LG_FREE_WORK
#undef  LG_FREE_ALL
#define LG_FREE_ALL ;

static inline int demo_init (bool burble)
{
    char msg [LAGRAPH_MSG_LEN] ;
    msg [0] = '\0' ;

    #ifdef __GLIBC__
    // Use mallopt to speedup malloc and free on Linux.  Otherwise, it can take
    // several seconds to free a large block of memory.  For this to be
    // effective, demo_init must be called before calling malloc/free, and
    // before calling LAGraph_Init.
    mallopt (M_MMAP_MAX, 0) ;           // disable mmap; it's too slow
    mallopt (M_TRIM_THRESHOLD, -1) ;    // disable sbrk trimming
    mallopt (M_TOP_PAD, 16*1024*1024) ; // increase padding to speedup malloc
    #endif

    #ifdef _OPENMP
    int bind = omp_get_proc_bind ( ) ;
    int nplaces = omp_get_num_places ( ) ;
    printf ("OPENMP: nplaces %d, bind: %d\n", nplaces, bind) ;
    fprintf (stderr, "OPENMP: nplaces %d, bind: %d\n", nplaces, bind) ;
    fflush (stdout) ;
    fflush (stderr) ;
    #endif

#if defined ( EXPERIMENTAL_GPU )
    // use the GPU
    LAGRAPH_TRY (LAGr_Init (GxB_NONBLOCKING_GPU, NULL, NULL, NULL, NULL, NULL));
#else
    // just use the CPU
    LAGRAPH_TRY (LAGraph_Init (NULL)) ;
#endif

    char library [256], date [256] ;
    GRB_TRY (GrB_get (GrB_GLOBAL, library, GrB_NAME)) ;
    printf ("library: %s ", library) ;

    int32_t version [3] ;
    GRB_TRY (GrB_get (GrB_GLOBAL, &(version [0]), GrB_LIBRARY_VER_MAJOR)) ;
    GRB_TRY (GrB_get (GrB_GLOBAL, &(version [1]), GrB_LIBRARY_VER_MINOR)) ;
    GRB_TRY (GrB_get (GrB_GLOBAL, &(version [2]), GrB_LIBRARY_VER_PATCH)) ;
    printf ("v%d.%d.%d ", version [0], version [1], version [2]) ;

    date [0] = '\0' ;
    GRB_TRY (LG_GET_LIBRARY_DATE (date)) ;
    printf ("[%s]\n", date) ;
    GRB_TRY (LG_SET_BURBLE (burble)) ;

    return (GrB_SUCCESS) ;
}

#undef  LG_FREE_ALL
#endif
