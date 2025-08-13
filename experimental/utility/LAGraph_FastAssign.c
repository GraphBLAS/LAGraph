
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

// This method is allows for assigns or build into a vector. It is created to 
// target some specific mxv GraphBLAS methods that occasionally tend to be 
// faster than both build and assign. 

// This method is used for cases where you want to build a vector equivilent to 
// the following for loop:
//      for (j = 0 ; j < n ; j++)
//      {
//          uint64_t i = I_vec [j] ;
//          c [i] += X_vec [j] ;
//      }

// It is fastest when the accum biop is equivalent to the dup moniod and c is a 
// full vector or when there is no accumulator and I_vec.nvals is 
// sufficiently large when compared to c.nrows. 

// It builds a matix P which is coustructed in O(1) time if a 
// filled ramp vector is passed in. 

// P is built by column and will contain one entry per column placed at row 
// I_vec [j]. So, P*X_vec will find for every row the intersections with X_vec,
// apply the dup monoid and put the result with the accumulator back into c at 
// that row. 

// X_vec can be sparse if dup is a semiring (dup_second), in which case the jth
// column correspoding to an empty X_vec[j] have no effect on the outcome.

// desc can also be modified to affect the mask or to pass in GRB_TRAN for
// GrB_INP0, this has some interesting usecases, specifically when X_vec is a 
// "map". That is:
// LAGraph_FastAssign(c, NULL, NULL, I_vec, map, NULL, any_second, NULL, msg)
// is equivalent to:
//      for (j = 0 ; j < n ; j++)
//      {
//          c [j] = map[I_vec[j]] ;
//      }
// Currently, this is not implemented for LAGraph_FastAssign_Monoid so a 
// semiring must be passed in.


#include "LG_internal.h"
#include "LAGraphX.h"
#if LG_SUITESPARSE_GRAPHBLAS_V10

// Uncomment if you would like to use the monoid version of FastAssign.
// Passing in a semiring is faster but this may be more convienient.
#if 0
#undef LG_FREE_ALL
#define LG_FREE_ALL                                           \
{                                                             \
    GrB_free(&sem);                                           \
}                                                     

int LAGraph_FastAssign_Monoid
(
    // output
    // Vector to be built (or assigned): initialized with correct dimensions.
    GrB_Vector c, 
    // inputs
    const GrB_Vector mask,
    const GrB_BinaryOp accum, 
    const GrB_Vector I_vec, // Indecies  (duplicates allowed)
    const GrB_Vector X_vec, // Values
    // Optional (Give me a ramp with size > X_vec.size for faster calculations) 
    const GrB_Vector ramp, 
    const GrB_Monoid dup, // Applied to duplicates
    const GrB_Descriptor desc,
    char *msg
)
{
    GrB_BinaryOp op = NULL;
    GrB_Semiring sem = NULL;
    int code = 0;
    GrB_get(X_vec, &code, GrB_EL_TYPE_CODE);
    switch (code)
    {
        case GrB_BOOL_CODE   : op = GrB_SECOND_BOOL   ; break ;
        case GrB_INT8_CODE   : op = GrB_SECOND_INT8   ; break ;
        case GrB_INT16_CODE  : op = GrB_SECOND_INT16  ; break ;
        case GrB_INT32_CODE  : op = GrB_SECOND_INT32  ; break ;
        case GrB_INT64_CODE  : op = GrB_SECOND_INT64  ; break ;
        case GrB_UINT8_CODE  : op = GrB_SECOND_UINT8  ; break ;
        case GrB_UINT16_CODE : op = GrB_SECOND_UINT16 ; break ;
        case GrB_UINT32_CODE : op = GrB_SECOND_UINT32 ; break ;
        case GrB_UINT64_CODE : op = GrB_SECOND_UINT64 ; break ;
        case GrB_FP32_CODE   : op = GrB_SECOND_FP32   ; break ;
        case GrB_FP64_CODE   : op = GrB_SECOND_FP64   ; break ;
        case GxB_FC32_CODE   : op = GxB_SECOND_FC32   ; break ;
        case GxB_FC64_CODE   : op = GxB_SECOND_FC64   ; break ;
        default : 
            LG_ERROR_MSG("LAGraph failed (file %s, line %d):" \
            " LAGraph_FastAssign_Monoid not implemented for UDTs", 
            __FILE__, __LINE__);
            break ;
    }
    GRB_TRY (GrB_Semiring_new(&sem, dup, op)) ;
    // GxB_print(sem, stdout) ;
    LG_TRY (LAGraph_FastAssign_Semiring 
        (c, mask, accum, I_vec, X_vec, ramp, sem, desc, msg)) ;
    LG_FREE_ALL ;
    return (GrB_SUCCESS);
}
#endif

#undef LG_FREE_ALL
#define LG_FREE_ALL                                           \
{                                                             \
    GrB_free(&P);                                             \
    GrB_free(&con);                                           \
    LAGraph_Free(&ramp_a, msg);                               \
    LAGraph_Free(&i_a, msg);                                  \
}                                                     

// This method can be faster if given a builtin semiring. 
int LAGraph_FastAssign_Semiring
(
    // output
    // Vector to be built (or assigned): initialized with correct dimensions.
    GrB_Vector c, 
    // inputs
    const GrB_Vector mask,
    const GrB_BinaryOp accum, 
    const GrB_Vector I_vec, // Indecies  (duplicates allowed)
    const GrB_Vector X_vec, // Values
    // Optional (Give me a ramp with size > X_vec.size for faster calculations) 
    const GrB_Vector ramp, 
    // monoid is applied to duplicates. Binary op should be SECOND.
    const GrB_Semiring semiring, 
    const GrB_Descriptor desc,
    char *msg
)
{
    // TODO: Let indicies be specified by value of by index in I_vec via 
    // descriptor (although build or assign would be better if I_vec is by 
    // index since it is sorted and has no dups). By value could be useful if 
    // I_vec is not full.
    // TODO: Ditto for X_vec

    GrB_Matrix P = NULL;
    GrB_Index n, nrows;
    GxB_Container con = NULL;
    void *ramp_a = NULL, *i_a =NULL;
    int ramp_h = GrB_DEFAULT, trsp = GrB_DEFAULT, i_h = GrB_DEFAULT;
    uint64_t ramp_n = 0, ramp_size = 0, i_n = 0, i_size= 0;
    GrB_Type i_type = NULL, ramp_type = NULL;
    //----------------------------------------------------------------------
    // Check inputs
    //----------------------------------------------------------------------
    //TODO: assert inputs are full or desc says to use by value etc.
    LG_ASSERT (c != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (I_vec != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (X_vec != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT_MSG (c != X_vec && c != I_vec && c != ramp && c != mask, 
        GrB_NOT_IMPLEMENTED, "c cannot be aliased with any input.") ; 

    //----------------------------------------------------------------------
    // Find dimensions and type
    //----------------------------------------------------------------------
    GRB_TRY (GrB_Vector_size(&n, I_vec)) ;

    if(desc != NULL)
    {
        GRB_TRY (GrB_get(desc, &trsp, GrB_INP0)) ;
    }
    if(trsp == GrB_TRAN)
    {
        GRB_TRY (GrB_Vector_size(&nrows, X_vec)) ;
    }
    else 
    {
        GRB_TRY (GrB_Vector_size(&nrows, c)) ;
    }

    //----------------------------------------------------------------------
    // Load up containers
    //----------------------------------------------------------------------
    GRB_TRY (GrB_Matrix_new(&P, GrB_BOOL, nrows, n));
    GRB_TRY (GxB_Container_new(&con));

    if(ramp == NULL)
    {
        //FUTURE: maybe let user input a size 0 ramp and build it for them?
        GRB_TRY (GrB_free(&(con->p))) ;
        ramp_type = (n + 1 <= INT32_MAX)? GrB_UINT32: GrB_UINT64;
        GrB_IndexUnaryOp idxnum = (n + 1 <= INT32_MAX)? 
                GrB_ROWINDEX_INT32: GrB_ROWINDEX_INT64;
        GRB_TRY (GrB_Vector_new(&(con->p), ramp_type, n + 1));
        GRB_TRY (GrB_assign (con->p, NULL, NULL, 0, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_apply (con->p, NULL, NULL, idxnum, con->p, 0, NULL)) ;
    }
    else
    {
        GRB_TRY (GxB_Vector_unload(
            ramp, &ramp_a, &ramp_type, &ramp_n, &ramp_size, &ramp_h, NULL)) ;
        LG_ASSERT_MSG (ramp_n > n, GrB_DIMENSION_MISMATCH, "Ramp too small!");
        GRB_TRY (GxB_Vector_load(
            con->p, &ramp_a, ramp_type, n + 1, ramp_size,
            GxB_IS_READONLY, NULL)) ;
        // Since con->p won't free this array I should be safe to load it back 
        // into ramp.
        GRB_TRY (GxB_Vector_load(
            ramp, &ramp_a, ramp_type, ramp_n, ramp_size, ramp_h, NULL)) ;
        ramp_a = NULL;
    }
    // con->i = I_vec;
    GRB_TRY (GxB_Vector_unload(
        I_vec, &i_a, &i_type, &i_n, &i_size, &i_h, NULL)) ;
    GRB_TRY (GxB_Vector_load(
        con->i, &i_a, i_type, i_n, i_size, GxB_IS_READONLY, NULL)) ;
    // Since con->i won't free this array I should be safe to load it back 
    // into I_vec.
    GRB_TRY (GxB_Vector_load(
        I_vec, &i_a, i_type, i_n, i_size, i_h, NULL)) ;
    i_a = NULL;

    // con->x [0] = false, of length 1
    GRB_TRY (GrB_free(&(con->x))) ;
    GRB_TRY (GrB_Vector_new (&(con->x), GrB_BOOL, 1)) ;
    GRB_TRY (GrB_assign (con->x, NULL, NULL, 0, GrB_ALL, 1, NULL)) ;
    con->format = GxB_SPARSE;
    con->orientation = GrB_COLMAJOR;
    con->nrows = nrows;
    con->ncols = n ;
    con->nvals = n ;
    con->nrows_nonempty = -1 ;
    con->ncols_nonempty = n ;
    con->iso = true ;
    con->jumbled = false ;
    con->format = GxB_SPARSE ;
    con->orientation = GrB_COLMAJOR ;
    con->Y = NULL ;
    //----------------------------------------------------------------------
    // Load P and do the mxv
    //----------------------------------------------------------------------
    GRB_TRY (GxB_load_Matrix_from_Container(P, con, NULL));
    // GRB_TRY (GxB_fprint(P, GxB_SHORT, stdout));
    GRB_TRY (GrB_mxv(c, mask, accum, semiring, P, X_vec, desc));
    //----------------------------------------------------------------------
    // Free work. 
    // Note: this does not free inputs since they are marked GxB_IS_READONLY
    //----------------------------------------------------------------------
    GrB_free(&P) ;
    GrB_free(&con) ;
    return (GrB_SUCCESS) ;
}
#endif
