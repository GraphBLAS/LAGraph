//------------------------------------------------------------------------------
// GB_msort_2.h: definitions for GB_msort_2.c
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2019, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// A parallel mergesort of an array of 2-by-n integers.  Each key consists
// of two integers.

#ifndef GB_MSORT_2_H
#define GB_MSORT_2_H

//------------------------------------------------------------------------------
// prototypes only needed for GB_msort_2
//------------------------------------------------------------------------------

void GB_merge_parallel_2                   // parallel merge
(
    int64_t *LAGRAPH_RESTRICT S_0,              // output of length nbigger + nsmaller
    int64_t *LAGRAPH_RESTRICT S_1,
    const int64_t *LAGRAPH_RESTRICT Bigger_0,   // Bigger [0..nbigger-1]
    const int64_t *LAGRAPH_RESTRICT Bigger_1,
    const int64_t nbigger,
    const int64_t *LAGRAPH_RESTRICT Smaller_0,  // Smaller [0..nsmaller-1]
    const int64_t *LAGRAPH_RESTRICT Smaller_1,
    const int64_t nsmaller
) ;

void GB_merge_select_2      // parallel or sequential merge of 2-by-n arrays
(
    int64_t *LAGRAPH_RESTRICT S_0,              // output of length nleft+nright
    int64_t *LAGRAPH_RESTRICT S_1,
    const int64_t *LAGRAPH_RESTRICT Left_0,     // Left [0..nleft-1]
    const int64_t *LAGRAPH_RESTRICT Left_1,
    const int64_t nleft,
    const int64_t *LAGRAPH_RESTRICT Right_0,    // Right [0..nright-1]
    const int64_t *LAGRAPH_RESTRICT Right_1,
    const int64_t nright
) ;

void GB_mergesort_2 // sort array A of size 2-by-n, using 2 keys (A [0:1][])
(
    int64_t *LAGRAPH_RESTRICT A_0,      // size n array
    int64_t *LAGRAPH_RESTRICT A_1,      // size n array
    int64_t *LAGRAPH_RESTRICT W_0,      // size n array, workspace
    int64_t *LAGRAPH_RESTRICT W_1,      // size n array, workspace
    const int64_t n
) ;


#define GB_BASECASE (64 * 1024)

//------------------------------------------------------------------------------
// GB_lt_2: sorting comparator function, two keys
//------------------------------------------------------------------------------

// A [a] and B [b] are keys of two integers.

// GB_lt_2 returns true if A [a] < B [b], for GB_qsort_2 and GB_msort_2

#define GB_lt_2(A_0, A_1, a, B_0, B_1, b)                                   \
(                                                                           \
    (A_0 [a] < B_0 [b]) ?                                                   \
    (                                                                       \
        true                                                                \
    )                                                                       \
    :                                                                       \
    (                                                                       \
        (A_0 [a] == B_0 [b]) ?                                              \
        (                                                                   \
            /* primary key is the same; tie-break on the 2nd key */         \
            (A_1 [a] < B_1 [b])                                             \
        )                                                                   \
        :                                                                   \
        (                                                                   \
            false                                                           \
        )                                                                   \
    )                                                                       \
)

#endif
