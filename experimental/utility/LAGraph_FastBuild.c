//------------------------------------------------------------------------------
// LAGraph_Fast_Build: Uses saxpy methods for faster builds, especially powerful
// when output is bitmap.
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.


#include "LG_internal.h"
#include "LAGraphX.h"

#include <omp.h>
#undef LG_FREE_ALL
#define LG_FREE_ALL                                           \
{                                                             \
    GrB_free(&P);                                             \
    GrB_free(&ramp);                                          \
}                                                     

int LAGraph_Fast_Build
(
    GrB_Vector c, // Vector to be built: initialized with correct dimensions.
    GrB_Vector i, // Indecies 
    GrB_Vector x, // Values
    // GrB_Vector ramp, // Optional (makes P load O(1))
    GrB_Monoid dup, // Applied to duplicates
    char *msg
)
{
    GrB_Vector ramp;
    GrB_Matrix P;
    int64_t n, nrows;
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
    GrB_Type x_type;
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
    GxB_Container con;
    GRB_TRY (GxB_Container_new(&con));
    GrB_Vector temp;
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
        c, NULL, NULL, dup, P, NULL)) ;
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
    GrB_free(&con);
}
