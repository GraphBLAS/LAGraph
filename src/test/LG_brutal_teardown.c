//------------------------------------------------------------------------------
// LG_brutal_teardown.c: teardown an LAGraph test with brutal memory testing
// -----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"
#include "LG_test.h"

int LG_brutal_teardown (char *msg)
{
    LG_TRY (LAGraph_Finalize (msg)) ;
    // nothing must be left allocated
    if (LG_nmalloc != 0) printf ("Leak! %g\n", (double) LG_nmalloc) ;
    return ((LG_nmalloc == 0) ? 0 : -911) ;
}

