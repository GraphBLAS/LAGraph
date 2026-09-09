//------------------------------------------------------------------------------
// LAGraph_Finalize: finish LAGraph
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

#if LAGRAPH_SUITESPARSE
#if GxB_IMPLEMENTATION <= GxB_VERSION (10,4,0)
void GB_Global_GrB_init_called_set (bool GrB_init_called) ;
#endif
#endif

int LAGraph_Finalize (char *msg)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    // ensure LAGraph_Finalized has not already been called
    LG_ASSERT_MSG (LG_get_LAGr_Init_has_been_called ( ), GrB_INVALID_VALUE,
        "LAGraph is already finalized") ;

    //--------------------------------------------------------------------------
    // free global objects
    //--------------------------------------------------------------------------

    LG_TRY (LG_Random_Finalize (msg)) ;

    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_first_int8  )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_first_int16 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_first_int32 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_first_int64 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_first_uint8 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_first_uint16)) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_first_uint32)) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_first_uint64)) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_first_fp32  )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_first_fp64  )) ;

    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_second_int8  )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_second_int16 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_second_int32 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_second_int64 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_second_uint8 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_second_uint16)) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_second_uint32)) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_second_uint64)) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_second_fp32  )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_second_fp64  )) ;

    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_one_int8  )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_one_int16 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_one_int32 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_one_int64 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_one_uint8 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_one_uint16)) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_one_uint32)) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_one_uint64)) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_one_fp32  )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_plus_one_fp64  )) ;

    GRB_TRY (GrB_Semiring_free (&LAGraph_any_one_bool  )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_any_one_int8  )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_any_one_int16 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_any_one_int32 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_any_one_int64 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_any_one_uint8 )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_any_one_uint16)) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_any_one_uint32)) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_any_one_uint64)) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_any_one_fp32  )) ;
    GRB_TRY (GrB_Semiring_free (&LAGraph_any_one_fp64  )) ;

    //--------------------------------------------------------------------------
    // finalize GraphBLAS
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_finalize ( )) ;
    #if LAGRAPH_SUITESPARSE
    #if GxB_IMPLEMENTATION <= GxB_VERSION (10,4,0)
    GB_Global_GrB_init_called_set (false) ;
    #endif
    #endif
    LG_set_LAGr_Init_has_been_called (false) ;
    return (GrB_SUCCESS) ;
}

