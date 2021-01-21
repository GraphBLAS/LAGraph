//------------------------------------------------------------------------------
// LAGraph_Test.h: include file for LAGraph/Test2
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M University

#ifndef LAGRAPH_TEST_H
#define LAGRAPH_TEST_H

#include "LAGraph2.h"

#ifndef LAGraph_CATCH
#define LAGraph_CATCH(status)                                               \
{                                                                           \
    printf ("LAGraph error: %s line: %d, status: %d: %s\n", __FILE__,       \
        __LINE__, status, msg) ;                                            \
    LAGRAPH_FREE_ALL ;                                                      \
    return (-1) ;                                                           \
}
#endif

#ifndef GrB_CATCH
#define GrB_CATCH(info)                                                     \
{                                                                           \
    printf ("GraphBLAS error: %s line: %d, info: %d: %s\n", __FILE__,       \
        __LINE__, info, msg) ;                                              \
    LAGRAPH_FREE_ALL ;                                                      \
    return (-1) ;                                                           \
}
#endif

#endif

