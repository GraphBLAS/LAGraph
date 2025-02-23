
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
    GrB_free(&ramp);                                          \
    GrB_free(&con);                                           \
    GrB_free(&temp);                                          \
}                                                     

int LAGraph_FastAssign
(
    GrB_Vector c, // Vector to be built (or assigned): initialized with correct dimensions.
    GrB_Vector mask,
    GrB_BinaryOp accum, 
    GrB_Vector i, // Indecies  (duplicates allowed)
    GrB_Vector x, // Values
    // GrB_Vector ramp, // Optional (makes P load O(1))
    GrB_Monoid dup, // Applied to duplicates
    char *msg
)
{
    GrB_Vector ramp = NULL;
    GrB_Matrix P = NULL;
    int64_t n, nrows;
    GxB_Container con = NULL;
    GrB_Vector temp = NULL;

    bool iso;
    //TODO: allow user to input a ramp for faster times
    //TODO: assert inputs are full etc
    LG_ASSERT (c != NULL, GrB_NULL_POINTER);
    LG_ASSERT (i != NULL, GrB_NULL_POINTER);
    LG_ASSERT (x != NULL, GrB_NULL_POINTER);

    GRB_TRY (GrB_Vector_size(&n, i));
    GRB_TRY (GrB_Vector_size(&nrows, c));
    GRB_TRY (GrB_Vector_get_INT32(x, (int32_t *) &iso, GxB_ISO));

    GrB_Type ramp_type = (n + 1 <= INT32_MAX)? GrB_UINT32: GrB_UINT64;
    GrB_Type x_type = NULL;
    char typename[LAGRAPH_MAX_NAME_LEN];
    LG_TRY (LAGraph_Vector_TypeName(typename, x, msg));
    LG_TRY (LAGraph_TypeFromName (&x_type, typename, msg)) ;

    GrB_IndexUnaryOp idxnum = (n <= INT32_MAX)? 
                GrB_ROWINDEX_INT32: GrB_ROWINDEX_INT64;
    GRB_TRY (GrB_Vector_new(&ramp, ramp_type, n + 1));
    GRB_TRY (GrB_Matrix_new(&P, x_type, nrows, n));
    GRB_TRY (GrB_assign (ramp, NULL, NULL, 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_apply (ramp, NULL, NULL, idxnum, ramp, 0, NULL)) ;
    // GxB_fprint(ramp, GxB_COMPLETE, stdout);
    GRB_TRY (GxB_Container_new(&con));
    temp = con->p;
    con->p = ramp;
    ramp = temp;
    temp = con->i;
    con->i = i;
    i = temp;
    temp = con->x;
    con->x = x;
    x = temp;
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
    temp = con->p;
    con->p = ramp;
    ramp = temp;
    temp = con->i;
    con->i = i;
    i = temp;
    temp = con->x;
    con->x = x;
    x = temp;
    LG_FREE_ALL;
}
#endif
