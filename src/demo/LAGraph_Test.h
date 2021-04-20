//------------------------------------------------------------------------------
// LAGraph_Test.h: include file for LAGraph/Test2
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M University

#ifndef LAGRAPH_TEST_H
#define LAGRAPH_TEST_H

#include <LAGraph.h>

#ifndef LAGraph_CATCH
#define LAGraph_CATCH(status)                                               \
{                                                                           \
    printf ("LAGraph error: %s line: %d, status: %d: [%s]\n", __FILE__,     \
        __LINE__, status, msg) ;                                            \
    LAGRAPH_FREE_ALL ;                                                      \
    return (-1) ;                                                           \
}
#endif

#ifndef GrB_CATCH
#define GrB_CATCH(info)                                                     \
{                                                                           \
    printf ("GraphBLAS error: %s line: %d, info: %d: [%s]\n", __FILE__,     \
        __LINE__, info, msg) ;                                              \
    LAGRAPH_FREE_ALL ;                                                      \
    return (-1) ;                                                           \
}
#endif

// LAGraph_Test_ReadProblem: read in a graph from a file
int LAGraph_Test_ReadProblem    // returns 0 if successful, -1 if failure
(
    // output
    LAGraph_Graph *G,           // graph from the file
    GrB_Matrix *SourceNodes,    // source nodes
    // inputs
    bool make_symmetric,        // if true, always return G as undirected
    bool remove_self_edges,     // if true, remove self edges
    bool pattern,               // if true, return G->A as bool (all true)
    GrB_Type pref,              // if non-NULL, typecast G->A to this type
    bool ensure_positive,       // if true, ensure all entries are > 0
    int argc,                   // input to main test program
    char **argv,                // input to main test program
    char *msg
) ;

#endif
