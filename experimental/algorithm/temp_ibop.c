#include "LG_internal.h"
#include <LAGraphX.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

typedef struct
{
    int64_t k;
    double v;
} tuple_fp64;
#define FP64_K "typedef struct { int64_t k ; double v ; } tuple_fp64 ;"
void make_fp64(tuple_fp64 *z,
               const double *x, GrB_Index ix, GrB_Index jx,
               const void *y, GrB_Index iy, GrB_Index jy,
               const void *theta)
{
    z->k = (int64_t)jx;
    z->v = (*x);
}
void max_fp64(tuple_fp64 *z, const tuple_fp64 *x, const tuple_fp64 *y)
{
    if (x->v > y->v || (x->v == y->v && x->k < y->k))
    {
        z->k = x->k;
        z->v = x->v;
    }
    else
    {
        z->k = y->k;
        z->v = y->v;
    }
}

#define MAX_FP64 (a string containing the max_fp64 function above)
// create the types and operators:
GrB_Scalar Theta;
// unused, but cannot be NULL
GrB_Scalar_new(&Theta, GrB_BOOL);
GrB_Scalar_setElement_BOOL(Theta, 0);
GzB_IndexBinaryOp Iop;
GrB_BinaryOp Bop;
GrB_Type Tuple;
GxB_Type_new(&Tuple, sizeof(tuple_fp64), "tuple_fp64", FP64_K);
 GzB_IndexBinaryOp_new (&Iop, make_fp64, Tuple, GrB_FP64, GrB_BOOL, GrB_BOOL,
 "make_fp64", MAKE_FP64)) ;
 OK(GzB_BinaryOp_new_IndexOp(&Bop, Iop, Theta));
 tuple_fp64 id;
 memset(&id, 0, sizeof(tuple_fp64));
 id.k = INT64_MAX;
 id.v = (double)(-INFINITY);
 OK(GxB_BinaryOp_new(&MonOp, max_fp64, Tuple, Tuple, Tuple, "max_fp64", MAX_FP64));
 GrB_Monoid MonOp;
 GrB_Semiring Semiring;
 GrB_Monoid_new_UDT(&Monoid, MonOp, &id);
 GrB_Semiring_new(&Semiring, Monoid, Bop);
 // compute the argmax of each row of a GrB_FP64 matrix A:
 // y = zeros (ncols,1) ;
 GrB_Vector y;
 GrB_Matrix_new (&y, GrB_BOOL, ncols, 1)) ;
 GrB_Matrix_assign_BOOL (y, NULL, NULL, 0, GrB_ALL, ncols, GrB_ALL, 1, NULL)) ;
 // c = A*y using the argmax semiring
 GrB_Vector_new (&c, Tuple, nrows, 1)) ;
 GrB_mxv(c, NULL, NULL, Semiring, A, y, NULL);