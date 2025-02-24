
//------------------------------------------------------------------------------
// LAGraph_Fast_Build: Uses saxpy methods for faster builds, especially powerful
// when output is bitmap.
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2024 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Gabriel Gomez, Texas A&M University

//------------------------------------------------------------------------------


#include "LG_internal.h"
#include "LAGraphX.h"
#if GxB_IMPLEMENTATION >= GxB_VERSION (10,0,0)
#include <omp.h>
#undef LG_FREE_ALL
#define LG_FREE_ALL                                           \
{                                                             \
    GrB_free(&P);                                             \
    GrB_free(&ramp_loc);                                      \
    GrB_free(&con);                                           \
    LAGraph_Free(&ramp_a, msg);                               \
}                                                     

int LAGraph_FastAssign
(
    // output
    // Vector to be built (or assigned): initialized with correct dimensions.
    GrB_Vector c, 
    // inputs
    const GrB_Vector mask,
    const GrB_BinaryOp accum, 
    const GrB_Vector i, // Indecies  (duplicates allowed)
    const GrB_Vector x, // Values
    // Optional (Give me a ramp with size > x.size for faster calculations) 
    const GrB_Vector ramp, 
    const GrB_Monoid dup, // Applied to duplicates
    char *msg
)
{
    // TODO: put data from ALL input vectors into a GxB_IS_READONLY vector
    // to be sure it remains completely unchanged? 
    GrB_Vector ramp_loc = NULL;
    GrB_Matrix P = NULL;
    int64_t n, nrows;
    GxB_Container con = NULL;
    void *ramp_a = NULL;
    int ramp_h = 0;
    int64_t ramp_n = 0, ramp_size = 0;

    bool iso = false;
    //TODO: assert inputs are full etc.
    LG_ASSERT (c != NULL, GrB_NULL_POINTER);
    LG_ASSERT (i != NULL, GrB_NULL_POINTER);
    LG_ASSERT (x != NULL, GrB_NULL_POINTER);
    LG_ASSERT_MSG (c != x, GrB_NOT_IMPLEMENTED, "c cannot be aliased with x.");   

    GRB_TRY (GrB_Vector_size(&n, i));
    GRB_TRY (GrB_Vector_size(&nrows, c));
    GRB_TRY (GrB_Vector_get_INT32(x, (int32_t *) &iso, GxB_ISO));

    GrB_Type x_type = NULL;
    char typename[LAGRAPH_MAX_NAME_LEN];
    LG_TRY (LAGraph_Vector_TypeName(typename, x, msg));
    LG_TRY (LAGraph_TypeFromName (&x_type, typename, msg)) ;

    GrB_Type ramp_type = (n + 1 <= INT32_MAX)? GrB_UINT32: GrB_UINT64;
    GrB_IndexUnaryOp idxnum = (n <= INT32_MAX)? 
            GrB_ROWINDEX_INT32: GrB_ROWINDEX_INT64;
    GRB_TRY (GrB_Vector_new(&ramp_loc, ramp_type, n + 1));
    if(ramp == NULL)
    {
        
        GRB_TRY (GrB_assign (ramp_loc, NULL, NULL, 0, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_apply (ramp_loc, NULL, NULL, idxnum, ramp_loc, 0, NULL)) ;
    }
    else
    {
        GRB_TRY (GxB_Vector_unload(
            ramp, &ramp_a, &ramp_type, &ramp_n, &ramp_size, &ramp_h, NULL)) ;
        LG_ASSERT (ramp_n > n, GrB_DIMENSION_MISMATCH);
        GRB_TRY (GxB_Vector_load(
            ramp_loc, &ramp_a, ramp_type, n + 1, (n + 1) * (ramp_size / ramp_n),
            GxB_IS_READONLY, NULL)) ;
    }
    GRB_TRY (GrB_Matrix_new(&P, x_type, nrows, n));
    // GxB_fprint(ramp_loc, GxB_COMPLETE, stdout);
    GRB_TRY (GxB_Container_new(&con));
    GRB_TRY (GrB_free(&con->p)) ;
    GRB_TRY (GrB_free(&con->i)) ;
    GRB_TRY (GrB_free(&con->x)) ;
    con->p = ramp_loc;
    if (c == i)
    {
        GRB_TRY (GrB_Vector_dup(&con->i, i)) ;
    }
    else
        con->i = i;
    con->x = x;
    con->format = GxB_SPARSE;
    con->orientation = GrB_COLMAJOR;
    con->nrows = nrows;
    con->ncols = n;
    con->iso = iso;
    con->nvals = n;
    con->jumbled = false;
    GRB_TRY (GxB_load_Matrix_from_Container(P, con, NULL));
    GRB_TRY (GrB_reduce(
        c, mask, accum, dup, P, NULL)) ;
    GRB_TRY (GxB_unload_Matrix_into_Container(P, con, NULL));
    // Don't let inputs get freed
    con->p = NULL;
    if (c != i)
        con->i = NULL;
    con->x = NULL;
    if (ramp)
    {
        GRB_TRY (GxB_Vector_load(
            ramp, &ramp_a, ramp_type, ramp_n, ramp_size, ramp_h, NULL)) ;
    }
    LG_FREE_ALL;
}
#endif
