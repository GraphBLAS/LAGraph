//------------------------------------------------------------------------------
// LG_qsort_1b: sort a 2-by-n list, using A [0][ ] as the sort key
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// These functions are not user-callable.  TODO: These functions do not yet
// have a user-callable interface, and they might not need one since they are
// currently unused.

#include "LG_internal.h"

// returns true if A [a] < B [b]
#define LG_lt(A,a,B,b) LG_lt_1 (A ## _0, a, B ## _0, b)

// each entry has a single key
#define LG_K 1

//------------------------------------------------------------------------------
// LG_qsort_1b: generic method for any data type
//------------------------------------------------------------------------------

// argument list for calling a function
#define LG_arg(A)                       \
    A ## _0, A ## _1, xsize

// argument list for calling a function, with offset
#define LG_arg_offset(A,x)              \
    A ## _0 + (x), A ## _1 + (x)*xsize, xsize

// argument list for defining a function
#define LG_args(A)                      \
    int64_t *LG_RESTRICT A ## _0,       \
    LG_void *LG_RESTRICT A ## _1,       \
    size_t xsize

// swap A [a] and A [b]
#define LG_swap(A,a,b)                                                        \
{                                                                             \
    int64_t t0 = A ## _0 [a] ; A ## _0 [a] = A ## _0 [b] ; A ## _0 [b] = t0 ; \
    memcpy (tx, A ## _1 + (a)*xsize, xsize) ;                                 \
    memcpy (A ## _1 + (a)*xsize, A ## _1 + (b)*xsize, xsize) ;                \
    memcpy (A ## _1 + (b)*xsize, tx, xsize) ;                                 \
}

#define LG_partition LG_partition_1b
#define LG_quicksort LG_quicksort_1b

#include "LG_qsort_template.h"

int LG_qsort_1b    // sort array A of size 2-by-n, using 1 key (A [0][])
(
    int64_t *LG_RESTRICT A_0,       // size n array
    LG_void *LG_RESTRICT A_1,       // size n array
    const size_t xsize,             // size of entries in A_1
    const int64_t n
)
{
    uint64_t seed = n ;
    LG_void *tx = LAGraph_Malloc (1, xsize) ;
    if (tx == NULL) return (-1) ;
    LG_quicksort (LG_arg (A), n, &seed, tx) ;
    LAGraph_Free ((void **) &tx) ;
    return (0) ;
}

//------------------------------------------------------------------------------
// LG_qsort_1b_size1:  quicksort with A_1 of type that has sizeof 1
//------------------------------------------------------------------------------

// for GrB_BOOL, GrB_INT8, GrB_UINT8, and user-defined types with sizeof(...)=1

#define A1_type uint8_t

// argument list for calling a function
#undef  LG_arg
#define LG_arg(A)                       \
    A ## _0, A ## _1

// argument list for calling a function, with offset
#undef  LG_arg_offset
#define LG_arg_offset(A,x)              \
    A ## _0 + (x), A ## _1 + (x)

// argument list for defining a function
#undef  LG_args
#define LG_args(A)                      \
    int64_t *LG_RESTRICT A ## _0,       \
    A1_type *LG_RESTRICT A ## _1        \

// swap A [a] and A [b]
#undef  LG_swap
#define LG_swap(A,a,b)                  \
{                                       \
    int64_t t0 = A ## _0 [a] ; A ## _0 [a] = A ## _0 [b] ; A ## _0 [b] = t0 ; \
    A1_type t1 = A ## _1 [a] ; A ## _1 [a] = A ## _1 [b] ; A ## _1 [b] = t1 ; \
}

#undef  LG_partition
#define LG_partition LG_partition_1b_size1
#undef  LG_quicksort
#define LG_quicksort LG_quicksort_1b_size1

#include "LG_qsort_template.h"

void LG_qsort_1b_size1  // LG_qsort_1b with A_1 with sizeof = 1
(
    int64_t *LG_RESTRICT A_0,       // size n array
    uint8_t *LG_RESTRICT A_1,       // size n array
    const int64_t n
)
{
    uint64_t seed = n ;
    LG_quicksort (LG_arg (A), n, &seed, NULL) ;
}

//------------------------------------------------------------------------------
// LG_qsort_1b_size2:  quicksort with A_1 of type that has sizeof 2
//------------------------------------------------------------------------------

// for GrB_INT16, GrB_UINT16, and user-defined types of sizeof(...) = 2

#undef  A1_type
#define A1_type uint16_t
#undef  LG_partition
#define LG_partition LG_partition_1b_size2
#undef  LG_quicksort
#define LG_quicksort LG_quicksort_1b_size2

#include "LG_qsort_template.h"

void LG_qsort_1b_size2  // LG_qsort_1b with A_1 with sizeof = 2
(
    int64_t *LG_RESTRICT A_0,       // size n array
    uint16_t *LG_RESTRICT A_1,      // size n array
    const int64_t n
)
{
    uint64_t seed = n ;
    LG_quicksort (LG_arg (A), n, &seed, NULL) ;
}

//------------------------------------------------------------------------------
// LG_qsort_1b_size4:  quicksort with A_1 of type that has sizeof 4
//------------------------------------------------------------------------------

// for GrB_INT32, GrB_UINT32, GrB_FP32, and user-defined types with
// sizeof(...) = 4.

#undef  A1_type
#define A1_type uint32_t
#undef  LG_partition
#define LG_partition LG_partition_1b_size4
#undef  LG_quicksort
#define LG_quicksort LG_quicksort_1b_size4

#include "LG_qsort_template.h"

void LG_qsort_1b_size4  // LG_qsort_1b with A_1 with sizeof = 4
(
    int64_t *LG_RESTRICT A_0,       // size n array
    uint32_t *LG_RESTRICT A_1,      // size n array
    const int64_t n
)
{
    uint64_t seed = n ;
    LG_quicksort (LG_arg (A), n, &seed, NULL) ;
}

//------------------------------------------------------------------------------
// LG_qsort_1b_size8:  quicksort with A_1 of type that has sizeof 8
//------------------------------------------------------------------------------

// for GrB_INT64, GrB_UINT64, GrB_FP64, GxB_FC32, and user-defined types
// with sizeof(...) = 8.

#undef  A1_type
#define A1_type uint64_t
#undef  LG_partition
#define LG_partition LG_partition_1b_size8
#undef  LG_quicksort
#define LG_quicksort LG_quicksort_1b_size8

#include "LG_qsort_template.h"

void LG_qsort_1b_size8  // LG_qsort_1b with A_1 with sizeof = 8
(
    int64_t *LG_RESTRICT A_0,       // size n array
    uint64_t *LG_RESTRICT A_1,      // size n array
    const int64_t n
)
{
    uint64_t seed = n ;
    LG_quicksort (LG_arg (A), n, &seed, NULL) ;
}

//------------------------------------------------------------------------------
// LG_qsort_1b_size16:  quicksort with A_1 of type that has sizeof 16
//------------------------------------------------------------------------------

// for GxB_FC64 and user-defined types with sizeof(...) = 16.

#undef  A1_type
#define A1_type LG_blob16
#undef  LG_partition
#define LG_partition LG_partition_1b_size16
#undef  LG_quicksort
#define LG_quicksort LG_quicksort_1b_size16

#include "LG_qsort_template.h"

void LG_qsort_1b_size16 // LG_qsort_1b with A_1 with sizeof = 16
(
    int64_t *LG_RESTRICT A_0,       // size n array
    LG_blob16 *LG_RESTRICT A_1,     // size n array
    const int64_t n
)
{
    ASSERT (sizeof (LG_blob16) == 16) ;
    uint64_t seed = n ;
    LG_quicksort (LG_arg (A), n, &seed, NULL) ;
}
