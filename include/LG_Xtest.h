//------------------------------------------------------------------------------
// LG_Xtest.h: include file for LAGraphX test library
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#ifndef LG_XTEST_H
#define LG_XTEST_H

#include <LAGraphX.h>

int LG_check_mis        // check if iset is a valid MIS of A
(
    GrB_Matrix A,
    GrB_Vector iset,
    char *msg
) ;

#endif
