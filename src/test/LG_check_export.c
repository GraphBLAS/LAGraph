//------------------------------------------------------------------------------
// LG_check_export: export G->A for testing
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Export G->A in CSR format, for testing only.
// See test_export for a brutal memory test of this method.

#define LAGraph_FREE_ALL                    \
{                                           \
    LAGraph_Free ((void **) Ap_handle) ;    \
    LAGraph_Free ((void **) Aj_handle) ;    \
    LAGraph_Free ((void **) Ax_handle) ;    \
}

#include "LG_internal.h"
#include "LG_test.h"

int LG_check_export
(
    // input
    LAGraph_Graph G,        // export G->A in CSR format
    // output
    GrB_Index **Ap_handle,  // size Ap_len on output
    GrB_Index **Aj_handle,  // size Aj_len on output
    void **Ax_handle,       // size Ax_len * typesize on output
    GrB_Index *Ap_len,
    GrB_Index *Aj_len,
    GrB_Index *Ax_len,
    size_t *typesize,       // size of the type of A
    char *msg
)
{
    LG_CLEAR_MSG ;

    GrB_Index *Ap = NULL, *Aj = NULL ;
    void *Ax = NULL ;
    LG_CHECK (LAGraph_CheckGraph (G, msg), -1002, "graph is invalid") ;
    LG_CHECK (Ap_handle == NULL || Aj_handle == NULL || Ax_handle == NULL
        || Ap_len == NULL || Aj_len == NULL || Ax_len == NULL
        || typesize == NULL, GrB_NULL_POINTER, "inputs are NULL") ;

    size_t s = 0 ;
    if      (G->A_type == GrB_BOOL  ) s = sizeof (bool    ) ;
    else if (G->A_type == GrB_INT8  ) s = sizeof (int8_t  ) ;
    else if (G->A_type == GrB_INT16 ) s = sizeof (int16_t ) ;
    else if (G->A_type == GrB_INT32 ) s = sizeof (int32_t ) ;
    else if (G->A_type == GrB_INT64 ) s = sizeof (int64_t ) ;
    else if (G->A_type == GrB_UINT8 ) s = sizeof (uint8_t ) ;
    else if (G->A_type == GrB_UINT16) s = sizeof (uint16_t) ;
    else if (G->A_type == GrB_UINT32) s = sizeof (uint32_t) ;
    else if (G->A_type == GrB_UINT64) s = sizeof (uint64_t) ;
    else if (G->A_type == GrB_FP32  ) s = sizeof (float   ) ;
    else if (G->A_type == GrB_FP64  ) s = sizeof (double  ) ;
    LG_CHECK (s == 0, -1, "unsupported type") ;
    (*typesize) = s ;

    GrB_TRY (GrB_Matrix_exportSize (Ap_len, Aj_len, Ax_len, GrB_CSR_FORMAT,
        G->A)) ;
    Ap = (GrB_Index *) LAGraph_Malloc (*Ap_len, sizeof (GrB_Index)) ;
    Aj = (GrB_Index *) LAGraph_Malloc (*Aj_len, sizeof (GrB_Index)) ;
    Ax = (void      *) LAGraph_Malloc (*Ax_len, s) ;
    (*Ap_handle) = Ap ;
    (*Aj_handle) = Aj ;
    (*Ax_handle) = Ax ;
    LG_CHECK (Ap == NULL || Aj == NULL || Ax == NULL, GrB_OUT_OF_MEMORY,
        "out of memory") ;

    if      (G->A_type == GrB_BOOL  )
    {
        GrB_TRY (GrB_Matrix_export (Ap, Aj, (bool     *) Ax,
            Ap_len, Aj_len, Ax_len, GrB_CSR_FORMAT, G->A)) ;
    }
    else if (G->A_type == GrB_INT8  )
    {
        GrB_TRY (GrB_Matrix_export (Ap, Aj, (int8_t   *) Ax,
            Ap_len, Aj_len, Ax_len, GrB_CSR_FORMAT, G->A)) ;
    }
    else if (G->A_type == GrB_INT16 )
    {
        GrB_TRY (GrB_Matrix_export (Ap, Aj, (int16_t  *) Ax, 
            Ap_len, Aj_len, Ax_len, GrB_CSR_FORMAT, G->A)) ;
    }
    else if (G->A_type == GrB_INT32 )
    {
        GrB_TRY (GrB_Matrix_export (Ap, Aj, (int32_t  *) Ax, 
            Ap_len, Aj_len, Ax_len, GrB_CSR_FORMAT, G->A)) ;
    }
    else if (G->A_type == GrB_INT64 )
    {
        GrB_TRY (GrB_Matrix_export (Ap, Aj, (int64_t  *) Ax,
            Ap_len, Aj_len, Ax_len, GrB_CSR_FORMAT, G->A)) ;
    }
    else if (G->A_type == GrB_UINT8 )
    {
        GrB_TRY (GrB_Matrix_export (Ap, Aj, (uint8_t  *) Ax,
            Ap_len, Aj_len, Ax_len, GrB_CSR_FORMAT, G->A)) ;
    }
    else if (G->A_type == GrB_UINT16)
    {
        GrB_TRY (GrB_Matrix_export (Ap, Aj, (uint16_t *) Ax,
            Ap_len, Aj_len, Ax_len, GrB_CSR_FORMAT, G->A)) ;
    }
    else if (G->A_type == GrB_UINT32)
    {
        GrB_TRY (GrB_Matrix_export (Ap, Aj, (uint32_t *) Ax,
            Ap_len, Aj_len, Ax_len, GrB_CSR_FORMAT, G->A)) ;
    }
    else if (G->A_type == GrB_UINT64)
    {
        GrB_TRY (GrB_Matrix_export (Ap, Aj, (uint64_t *) Ax,
            Ap_len, Aj_len, Ax_len, GrB_CSR_FORMAT, G->A)) ;
    }
    else if (G->A_type == GrB_FP32  )
    {
        GrB_TRY (GrB_Matrix_export (Ap, Aj, (float    *) Ax,
            Ap_len, Aj_len, Ax_len, GrB_CSR_FORMAT, G->A)) ;
    }
    else if (G->A_type == GrB_FP64  )
    {
        GrB_TRY (GrB_Matrix_export (Ap, Aj, (double   *) Ax,
            Ap_len, Aj_len, Ax_len, GrB_CSR_FORMAT, G->A)) ;
    }

    return (GrB_SUCCESS) ;
}
