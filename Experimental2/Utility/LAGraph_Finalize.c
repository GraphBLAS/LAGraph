//------------------------------------------------------------------------------
// LAGraph_Finalize: finish LAGraph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

#include "LAGraph_Internal.h"

int LAGraph_Finalize (char *msg)    // returns 0 if successful, -1 if failure
{
    // finalize GraphBLAS
    GrB_TRY (GrB_finalize ( )) ;
}

