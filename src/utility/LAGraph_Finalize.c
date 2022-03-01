//------------------------------------------------------------------------------
// LAGraph_Finalize: finish LAGraph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_Finalize (char *msg)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    //--------------------------------------------------------------------------
    // free global objects
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_first_int8  )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_first_int16 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_first_int32 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_first_int64 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_first_uint8 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_first_uint16)) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_first_uint32)) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_first_uint64)) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_first_fp32  )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_first_fp64  )) ;

    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_second_int8  )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_second_int16 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_second_int32 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_second_int64 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_second_uint8 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_second_uint16)) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_second_uint32)) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_second_uint64)) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_second_fp32  )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_second_fp64  )) ;

    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_one_int8  )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_one_int16 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_one_int32 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_one_int64 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_one_uint8 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_one_uint16)) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_one_uint32)) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_one_uint64)) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_one_fp32  )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_plus_one_fp64  )) ;

    GrB_TRY (GrB_Semiring_free (&LAGraph_structural_bool  )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_structural_int8  )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_structural_int16 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_structural_int32 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_structural_int64 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_structural_uint8 )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_structural_uint16)) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_structural_uint32)) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_structural_uint64)) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_structural_fp32  )) ;
    GrB_TRY (GrB_Semiring_free (&LAGraph_structural_fp64  )) ;

    //--------------------------------------------------------------------------
    // finalize GraphBLAS
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_finalize ( )) ;
    return (GrB_SUCCESS) ;
}

