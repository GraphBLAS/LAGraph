//------------------------------------------------------------------------------
// LAGraph/src/test/LG_get_vector: extract contents of a vector
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

#include "LAGraph_test.h"

bool LG_get_vector
(
    int64_t *x,
    GrB_Vector X,
    int64_t n,
    int64_t missing
)
{

    for (int64_t i = 0 ; i < n ; i++)
    {
        int64_t t ;
        int info = GrB_Vector_extractElement_INT64 (&t, X, i) ;
        TEST_CHECK (info == GrB_SUCCESS || info == GrB_NO_VALUE) ;
        if (info == GrB_SUCCESS)
        {
            x [i] = t ;
        }
        else
        {
            x [i] = missing ;
        }
    }
    return (true) ;             // success
}

