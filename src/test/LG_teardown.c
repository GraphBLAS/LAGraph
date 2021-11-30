//------------------------------------------------------------------------------
// LG_teardown.c: teardown an LAGraph test
// -----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#include "LG_internal.h"
#include "LG_test.h"

int LG_teardown (bool brutal_test, char *msg)
{
    LG_CHECK (LAGraph_Finalize (msg), -1, "finalize failed") ;
    if (brutal_test)
    {
        // nothing must be left allocated
        LG_CHECK (LG_nmalloc != 0, -1, "memory leak") ;
    }
}

