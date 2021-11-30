//------------------------------------------------------------------------------
// LG_setup.c: setup an LAGraph test
// -----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#include "LG_internal.h"
#include "LG_test.h"

int LG_setup (bool brutal_test, char *msg)
{
    LG_brutal = -1 ;        // disable brutal testing for now
    LG_nmalloc = 0 ;        // assuming nothing is malloc'd
    if (brutal_test)
    {
        return (LAGraph_Xinit (LG_check_malloc, LG_check_calloc,
            LG_check_realloc, LG_check_free, msg)) ;
    }
    else
    {
        return (LAGraph_Init (msg)) ;
    }
}

