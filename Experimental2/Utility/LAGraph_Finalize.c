//------------------------------------------------------------------------------
// LAGraph_Finalize: finish LAGraph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_Finalize (char *msg)    // returns 0 if successful, -1 if failure
{
    // finalize GraphBLAS
    GrB_TRY (GrB_finalize ( )) ;
}

