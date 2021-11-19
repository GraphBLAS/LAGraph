//------------------------------------------------------------------------------
// LAGraph_binread.h:
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite

// Contributed by Tim Davis, Texas A&M University

#ifndef LAGRAPH_BINREAD_H
#define LAGRAPH_BINREAD_H

#include <LAGraph.h>

#define CATCH(status)                                                         \
{                                                                             \
    printf ("error: %s line: %d, status: %d\n", __FILE__, __LINE__, status) ; \
    LAGraph_FREE_ALL ;                                                        \
    return (-1) ;                                                             \
}

#undef  LAGraph_CATCH
#define LAGraph_CATCH(status) CATCH (status)

#undef  GrB_CATCH
#define GrB_CATCH(info) CATCH (info)

#define ERROR CATCH(-1)
#define LAGRAPH_BIN_HEADER 512
#define LEN LAGRAPH_BIN_HEADER

// currently, SuiteSparse:GraphBLAS v6 or later is required
#if !defined ( GxB_SUITESPARSE_GRAPHBLAS )
#error "SuiteSparse v6 or later is required"
#endif

#define LAGraph_FREE_ALL                \
{                                       \
    GrB_free (A) ;                      \
    LAGraph_Free ((void **) &Ap) ;      \
    LAGraph_Free ((void **) &Ab) ;      \
    LAGraph_Free ((void **) &Ah) ;      \
    LAGraph_Free ((void **) &Ai) ;      \
    LAGraph_Free ((void **) &Ax) ;      \
}

//------------------------------------------------------------------------------
// binread: read a matrix from a binary file
//------------------------------------------------------------------------------

#define FREAD(p,s,n)                    \
{                                       \
    if (fread (p, s, n, f) != n)        \
    {                                   \
        ERROR ;                         \
    }                                   \
}

int LAGraph_binread   // returns 0 if successful, -1 if failure
(
    GrB_Matrix *A,          // matrix to read from the file
    GrB_Type *A_type,       // Scalar type
    FILE *f                 // file to read it from, already open
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Index *Ap = NULL, *Ai = NULL, *Ah = NULL ;
    int8_t *Ab = NULL ;
    void *Ax = NULL ;
    if (A == NULL || f == NULL) ERROR ;
    (*A) = NULL ;

    //--------------------------------------------------------------------------
    // basic matrix properties
    //--------------------------------------------------------------------------

    GxB_Format_Value fmt = -999 ;
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

    FREAD (&fmt,      sizeof (GxB_Format_Value), 1) ;
    FREAD (&kind,     sizeof (int32_t), 1) ;
    FREAD (&hyper,    sizeof (double), 1) ;
    FREAD (&nrows,    sizeof (GrB_Index), 1) ;
    FREAD (&ncols,    sizeof (GrB_Index), 1) ;
    FREAD (&nonempty, sizeof (int64_t), 1) ;
    FREAD (&nvec,     sizeof (GrB_Index), 1) ;
    FREAD (&nvals,    sizeof (GrB_Index), 1) ;
    FREAD (&typecode, sizeof (int32_t), 1) ;
    FREAD (&typesize, sizeof (size_t), 1) ;

    is_hyper  = (kind == 1) ;
    is_sparse = (kind == 0 || kind == GxB_SPARSE) ;
    is_bitmap = (kind == GxB_BITMAP) ;
    is_full   = (kind == GxB_FULL) ;

    // uniform-valued matrices not yet supported
    bool is_uniform = false ;

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
        #if 0
        case 11: type = GxB_FC32        ; break ;
        case 12: type = GxB_FC64        ; break ;
        #endif
        default: ERROR ;    // unknown or unsupported type
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
        Ap = LAGraph_Malloc (Ap_len, sizeof (GrB_Index)) ;
        Ah = LAGraph_Malloc (Ah_len, sizeof (GrB_Index)) ;
        Ai = LAGraph_Malloc (Ai_len, sizeof (GrB_Index)) ;
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
        Ap = LAGraph_Malloc (Ap_len, sizeof (GrB_Index)) ;
        Ai = LAGraph_Malloc (Ai_len, sizeof (GrB_Index)) ;
        Ap_size = Ap_len * sizeof (GrB_Index) ;
        Ai_size = Ai_len * sizeof (GrB_Index) ;
        ok = (Ap != NULL && Ai != NULL) ;
    }
    else if (is_bitmap)
    {
        Ab_len = nrows*ncols ;
        Ax_len = nrows*ncols ;
        Ab = LAGraph_Malloc (nrows*ncols, sizeof (int8_t)) ;
        Ab_size = Ab_len * sizeof (GrB_Index) ;
        ok = (Ab != NULL) ;
    }
    else if (is_full)
    {
        Ax_len = nrows*ncols ;
    }
    else
    {
        ERROR ;     // unknown matrix format
    }
    Ax = LAGraph_Malloc (Ax_len, typesize) ;
    Ax_size = Ax_len * typesize ;
    if (!ok) ERROR ;        // out of memory

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

    FREAD (Ax, typesize, Ax_len) ;

    //--------------------------------------------------------------------------
    // import the matrix
    //--------------------------------------------------------------------------

    if (fmt == GxB_BY_COL && is_hyper)
    {
        // hypersparse CSC
        GrB_TRY (GxB_Matrix_import_HyperCSC (A, type, nrows, ncols,
            &Ap, &Ah, &Ai, &Ax, Ap_size, Ah_size, Ai_size, Ax_size,
            is_uniform, nvec, false, NULL)) ;
    }
    else if (fmt == GxB_BY_ROW && is_hyper)
    {
        // hypersparse CSR
        GrB_TRY (GxB_Matrix_import_HyperCSR (A, type, nrows, ncols,
            &Ap, &Ah, &Ai, &Ax, Ap_size, Ah_size, Ai_size, Ax_size,
            is_uniform, nvec, false, NULL)) ;
    }
    else if (fmt == GxB_BY_COL && is_sparse)
    {
        // standard CSC
        GrB_TRY (GxB_Matrix_import_CSC (A, type, nrows, ncols,
            &Ap, &Ai, &Ax, Ap_size, Ai_size, Ax_size,
            is_uniform, false, NULL)) ;
    }
    else if (fmt == GxB_BY_ROW && is_sparse)
    {
        // standard CSR
        GrB_TRY (GxB_Matrix_import_CSR (A, type, nrows, ncols,
            &Ap, &Ai, &Ax, Ap_size, Ai_size, Ax_size,
            is_uniform, false, NULL)) ;
    }
    else if (fmt == GxB_BY_COL && is_bitmap)
    {
        // bitmap by col
        GrB_TRY (GxB_Matrix_import_BitmapC (A, type, nrows, ncols,
            &Ab, &Ax, Ab_size, Ax_size,
            is_uniform, nvals, NULL)) ;
    }
    else if (fmt == GxB_BY_ROW && is_bitmap)
    {
        // bitmap by row
        GrB_TRY (GxB_Matrix_import_BitmapR (A, type, nrows, ncols,
            &Ab, &Ax, Ab_size, Ax_size,
            is_uniform, nvals, NULL)) ;
    }
    else if (fmt == GxB_BY_COL && is_full)
    {
        // full by col
        GrB_TRY (GxB_Matrix_import_FullC (A, type, nrows, ncols,
            &Ax, Ax_size,
            is_uniform, NULL)) ;
    }
    else if (fmt == GxB_BY_ROW && is_full)
    {
        // full by row
        GrB_TRY (GxB_Matrix_import_FullR (A, type, nrows, ncols,
            &Ax, Ax_size,
            is_uniform, NULL)) ;
    }
    else
    {
        ERROR ;     // unknown format
    }

    GrB_TRY (GxB_set (*A, GxB_HYPER_SWITCH, hyper)) ;
    (*A_type) = type;
    type = NULL;
    return (0) ;
}


#endif
