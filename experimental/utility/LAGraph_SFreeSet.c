//------------------------------------------------------------------------------
// LAGraph_SFreeSet: free a set of matrices
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"
#include "LAGraphX.h"

void LAGraph_SFreeSet           // free a set of matrices
(
    // input/output
    GrB_Matrix **Set_handle,    // array of GrB_Matrix of size nmatrices
    GrB_Index nmatrices         // # of matrices in the set
)
{
    if (Set_handle != NULL)
    {
        GrB_Matrix *Set = (*Set_handle) ;
        if (Set != NULL)
        {
            for (GrB_Index i = 0 ; i < nmatrices ; i++)
            {
                GrB_free (&(Set [i])) ;
            }
        }
        LAGraph_Free ((void **) Set_handle, NULL) ;
    }
}

