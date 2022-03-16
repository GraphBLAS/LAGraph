//------------------------------------------------------------------------------
// LAGraph_SFreeContents: free the Contents returned by LAGraph_SRead.
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"
#include "LAGraphX.h"

void LAGraph_SFreeContents  // free the Contents returned by LAGraph_SRead
(
    // input/output
    LAGraph_Contents **Contents_handle,     // array of size ncontents
    GrB_Index ncontents
)
{
    if (Contents_handle != NULL)
    {
        LAGraph_Contents *Contents = (*Contents_handle) ;
        if (Contents != NULL)
        {
            for (GrB_Index i = 0 ; i < ncontents ; i++)
            {
                LAGraph_Free ((void **) &(Contents [i].blob), NULL) ;
            }
        }
        LAGraph_Free ((void **) Contents_handle, NULL) ;
    }
}

