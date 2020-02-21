//------------------------------------------------------------------------------
// LAGraph_complex:  complex number support for LAGraph
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// Contributed by Michel Pelletier. Adapted from 'usercomplex.c' code
// in SuiteSparse Demo by Dr. Tim Davis.

#include "LAGraph_internal.h"

// a global value for returning the complex type in a Matrix Market file:
GrB_Type LAGraph_Complex = NULL ;

#if defined __INTEL_COMPILER
#pragma warning (disable: 58 167 144 161 177 181 186 188 589 593 869 981 1418 1419 1572 1599 2259 2282 2557 2547 3280 )
#elif defined __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wincompatible-pointer-types"
#endif

#define C double complex
#define X *x
#define Y *y
#define Z *z

#define ONE  CMPLX(1,0)
#define ZERO CMPLX(0,0)
#define T ONE
#define F ZERO
#define BOOL(X) (X != ZERO)

//------------------------------------------------------------------------------
// 8 binary functions, z=f(x,y), where CxC -> C
//------------------------------------------------------------------------------

void complex_first  (C Z, const C X, const C Y) { Z = X ; }
void complex_second (C Z, const C X, const C Y) { Z = Y ; }
void complex_plus   (C Z, const C X, const C Y) { Z = X + Y ; }
void complex_minus  (C Z, const C X, const C Y) { Z = X - Y ; }
void complex_rminus (C Z, const C X, const C Y) { Z = Y - X ; }
void complex_times  (C Z, const C X, const C Y) { Z = X * Y ; }
void complex_div    (C Z, const C X, const C Y) { Z = X / Y ; }
void complex_rdiv   (C Z, const C X, const C Y) { Z = Y / X ; }

void complex_min (C Z, const C X, const C Y)
{
    // min (x,y): complex number with smallest magnitude.  If tied, select the
    // one with the smallest phase angle (same as MATLAB definition).
    // No special cases for NaNs.
    double absx = cabs (X) ;
    double absy = cabs (Y) ;
    if (absx < absy)
    {
        Z = X ;
    }
    else if (absx > absy)
    {
        Z = Y ;
    }
    else
    {
        if (carg (X) < carg (Y))
        {
            Z = X ;
        }
        else
        {
            Z = Y ;
        }
    }
}

void complex_max (C Z, const C X, const C Y)
{
    // max (x,y): complex number with largest magnitude.  If tied, select the
    // one with the largest phase angle (same as MATLAB definition).
    // No special cases for NaNs.
    double absx = cabs (X) ;
    double absy = cabs (Y) ;
    if (absx > absy)
    {
        Z = X ;
    }
    else if (absx < absy)
    {
        Z = Y ;
    }
    else
    {
        if (carg (X) > carg (Y))
        {
            Z = X ;
        }
        else
        {
            Z = Y ;
        }
    }
}

void complex_skew
(
    bool *z,
    const double complex *x,
    const double complex *y
)
{
    (*z) = (*x) == -(*y) ;
}

void complex_pair
(
    double complex *z,
    const double complex *x,
    const double complex *y
)
{
    (*z) = ONE ;
}

void complex_any
(
    double complex *z,
    const double complex *x,
    const double complex *y
)
{
    (*z) = (*y) ;
}

void complex_hermitian
(
    bool *z,
    const double complex *x,
    const double complex *y
)
{
    (*z) = (*x) == conj (*y) ;
}



GrB_BinaryOp
    LAGraph_FIRST_Complex = NULL          ,
    LAGraph_SECOND_Complex = NULL         ,
    LAGraph_MIN_Complex = NULL            ,
    LAGraph_MAX_Complex = NULL            ,
    LAGraph_PLUS_Complex = NULL           ,
    LAGraph_MINUS_Complex = NULL          ,
    LAGraph_TIMES_Complex = NULL          ,
    LAGraph_DIV_Complex = NULL            ,
    LAGraph_RMINUS_Complex = NULL         ,
    LAGraph_RDIV_Complex = NULL           ,
    LAGraph_SKEW_Complex = NULL           ,
    LAGraph_PAIR_Complex = NULL           ,
    LAGraph_ANY_Complex = NULL            ,
    LAGraph_HERMITIAN_Complex = NULL      ;

//------------------------------------------------------------------------------
// 6 binary functions, z=f(x,y), where CxC -> C ; (1,0) = true, (0,0) = false
//------------------------------------------------------------------------------

// inequality operators follow the MATLAB convention

#define R(x) creal(x)

void complex_iseq (C Z, const C X, const C Y) { Z = (X == Y) ? T : F ; }
void complex_isne (C Z, const C X, const C Y) { Z = (X != Y) ? T : F ; }
void complex_isgt (C Z, const C X, const C Y) { Z = (R(X) >  R(Y)) ? T : F ; }
void complex_islt (C Z, const C X, const C Y) { Z = (R(X) <  R(Y)) ? T : F ; }
void complex_isge (C Z, const C X, const C Y) { Z = (R(X) >= R(Y)) ? T : F ; }
void complex_isle (C Z, const C X, const C Y) { Z = (R(X) <= R(Y)) ? T : F ; }

GrB_BinaryOp
    LAGraph_ISEQ_Complex = NULL       ,
    LAGraph_ISNE_Complex = NULL       ,
    LAGraph_ISGT_Complex = NULL       ,
    LAGraph_ISLT_Complex = NULL       ,
    LAGraph_ISGE_Complex = NULL       ,
    LAGraph_ISLE_Complex = NULL       ;

//------------------------------------------------------------------------------
// binary boolean functions, z=f(x,y), where CxC -> C
//------------------------------------------------------------------------------

void complex_or (C Z, const C X, const C Y)
{
    Z = (BOOL (X) || BOOL (Y)) ? T : F ;
}

void complex_and (C Z, const C X, const C Y)
{
    Z = (BOOL (X) && BOOL (Y)) ? T : F ;
}

void complex_xor (C Z, const C X, const C Y)
{
    Z = (BOOL (X) != BOOL (Y)) ? T : F ;
}

GrB_BinaryOp
    LAGraph_OR_Complex = NULL        ,
    LAGraph_AND_Complex = NULL       ,
    LAGraph_XOR_Complex = NULL       ;

//------------------------------------------------------------------------------
// 6 binary functions, z=f(x,y), where CxC -> bool
//------------------------------------------------------------------------------

// inequality operators follow the MATLAB convention

void complex_eq (bool Z, const C X, const C Y) { Z = (X == Y) ; }
void complex_ne (bool Z, const C X, const C Y) { Z = (X != Y) ; }
void complex_gt (bool Z, const C X, const C Y) { Z = (R (X) >  R (Y)) ;}
void complex_lt (bool Z, const C X, const C Y) { Z = (R (X) <  R (Y)) ;}
void complex_ge (bool Z, const C X, const C Y) { Z = (R (X) >= R (Y)) ;}
void complex_le (bool Z, const C X, const C Y) { Z = (R (X) <= R (Y)) ;}

GrB_BinaryOp
    LAGraph_EQ_Complex = NULL        ,
    LAGraph_NE_Complex = NULL        ,
    LAGraph_GT_Complex = NULL        ,
    LAGraph_LT_Complex = NULL        ,
    LAGraph_GE_Complex = NULL        ,
    LAGraph_LE_Complex = NULL        ;

//------------------------------------------------------------------------------
// binary functions, z=f(x,y), where double x double -> complex
//------------------------------------------------------------------------------

void complex_complex (C Z, const double X, const double Y) { Z = CMPLX (X,Y) ; }

GrB_BinaryOp LAGraph_COMPLEX_Complex = NULL ;

//------------------------------------------------------------------------------
// unary functions, z=f(x) where C -> C
//------------------------------------------------------------------------------

void complex_one      (C Z, const C X) { Z =       1. ; }
void complex_identity (C Z, const C X) { Z =       X  ; }
void complex_ainv     (C Z, const C X) { Z =      -X  ; }
void complex_abs      (C Z, const C X) { Z = CMPLX (cabs (X), 0) ; }
void complex_minv     (C Z, const C X) { Z =  1. / X  ; } 
void complex_not      (C Z, const C X) { Z = BOOL (X) ? F : T ; }
void complex_conj     (C Z, const C X) { Z = conj (X) ; }

void complex_isone
(
    bool *z,
    const double complex *x
)
{
    (*z) = ((*x) == 1) ;
}


void complex_true_bool
(
    bool *z,
    const double complex *x     // ignored
)
{
    (*z) = true ;
}



GrB_UnaryOp
    LAGraph_IDENTITY_Complex = NULL        ,
    LAGraph_AINV_Complex = NULL            ,
    LAGraph_MINV_Complex = NULL            ,
    LAGraph_NOT_Complex = NULL             ,
    LAGraph_CONJ_Complex = NULL            ,
    LAGraph_ONE_Complex = NULL             ,
    LAGraph_ABS_Complex  = NULL            ,
    LAGraph_TRUE_BOOL_Complex = NULL       ,     
    LAGraph_ISONE_Complex  = NULL          ;

//------------------------------------------------------------------------------
// unary functions, z=f(x) where C -> double
//------------------------------------------------------------------------------

void complex_real  (double Z, const C X) { Z = creal (X) ; }
void complex_imag  (double Z, const C X) { Z = cimag (X) ; }
void complex_cabs  (double Z, const C X) { Z = cabs  (X) ; }
void complex_angle (double Z, const C X) { Z = carg  (X) ; }

GrB_UnaryOp
    LAGraph_REAL_Complex = NULL            ,
    LAGraph_IMAG_Complex = NULL            ,
    LAGraph_CABS_Complex = NULL            ,
    LAGraph_ANGLE_Complex = NULL           ;

//------------------------------------------------------------------------------
// unary functions, z=f(x) where double -> C
//------------------------------------------------------------------------------

void complex_complex_real (C Z, const double X) { Z = CMPLX (X, 0) ; }
void complex_complex_imag (C Z, const double X) { Z = CMPLX (0, X) ; }

GrB_UnaryOp
    LAGraph_COMPLEX_REAL_Complex = NULL    ,
    LAGraph_COMPLEX_IMAG_Complex = NULL    ;

//------------------------------------------------------------------------------
// Complex type, scalars, monoids, and semiring
//------------------------------------------------------------------------------

GrB_Monoid
    LAGraph_PLUS_Complex_MONOID = NULL     ,
    LAGraph_TIMES_Complex_MONOID = NULL    ;
    
GrB_Semiring LAGraph_PLUS_TIMES_Complex = NULL ;
C LAGraph_Complex_1  = ONE ;
C LAGraph_Complex_0 = ZERO ;

#define OK(method)                      \
    info = method ;                     \
    if (info != GrB_SUCCESS)            \
    {                                   \
        LAGraph_Complex_finalize ( ) ;  \
        return (info) ;                 \
    }

//------------------------------------------------------------------------------
// LAGraph_Complex_init: create the complex type, operators, monoids, and semiring
//------------------------------------------------------------------------------

GrB_Info LAGraph_Complex_init ( )
{

    GrB_Info info ;

    //--------------------------------------------------------------------------
    // create the Complex type
    //--------------------------------------------------------------------------

    OK (GrB_Type_new (&LAGraph_Complex, sizeof (C))) ;

    #undef C
    #undef D
    #define C LAGraph_Complex
    #define D GrB_FP64

    //--------------------------------------------------------------------------
    // create the Complex binary operators, CxC->C
    //--------------------------------------------------------------------------

    OK (GrB_BinaryOp_new (&LAGraph_FIRST_Complex       , complex_first       , C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_SECOND_Complex      , complex_second      , C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_MIN_Complex         , complex_min         , C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_MAX_Complex         , complex_max         , C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_PLUS_Complex        , complex_plus        , C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_MINUS_Complex       , complex_minus       , C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_RMINUS_Complex      , complex_rminus      , C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_TIMES_Complex       , complex_times       , C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_DIV_Complex         , complex_div         , C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_RDIV_Complex        , complex_rdiv        , C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_PAIR_Complex        , complex_pair        , C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_ANY_Complex         , complex_any         , C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_SKEW_Complex        , complex_skew        , GrB_BOOL, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_HERMITIAN_Complex   , complex_hermitian   , GrB_BOOL, C, C)) ;

    //--------------------------------------------------------------------------
    // create the Complex binary comparison operators, CxC -> C
    //--------------------------------------------------------------------------

    OK (GrB_BinaryOp_new (&LAGraph_ISEQ_Complex , complex_iseq ,  C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_ISNE_Complex , complex_isne ,  C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_ISGT_Complex , complex_isgt ,  C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_ISLT_Complex , complex_islt ,  C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_ISGE_Complex , complex_isge ,  C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_ISLE_Complex , complex_isle ,  C, C, C)) ;

    //--------------------------------------------------------------------------
    // create the Complex boolean operators, CxC -> C
    //--------------------------------------------------------------------------

    OK (GrB_BinaryOp_new (&LAGraph_OR_Complex  , complex_or  ,  C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_AND_Complex , complex_and ,  C, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_XOR_Complex , complex_xor ,  C, C, C)) ;

    //--------------------------------------------------------------------------
    // create the Complex binary operators, CxC -> bool
    //--------------------------------------------------------------------------

    OK (GrB_BinaryOp_new (&LAGraph_EQ_Complex , complex_eq ,  GrB_BOOL, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_NE_Complex , complex_ne ,  GrB_BOOL, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_GT_Complex , complex_gt ,  GrB_BOOL, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_LT_Complex , complex_lt ,  GrB_BOOL, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_GE_Complex , complex_ge ,  GrB_BOOL, C, C)) ;
    OK (GrB_BinaryOp_new (&LAGraph_LE_Complex , complex_le ,  GrB_BOOL, C, C)) ;

    //--------------------------------------------------------------------------
    // create the Complex binary operator, double x double -> C
    //--------------------------------------------------------------------------

    OK (GrB_BinaryOp_new (&LAGraph_COMPLEX_Complex, complex_complex, C, D, D)) ;

    //--------------------------------------------------------------------------
    // create the Complex unary operators, C->C
    //--------------------------------------------------------------------------

    OK (GrB_UnaryOp_new (&LAGraph_ONE_Complex       , complex_one       , C, C)) ;
    OK (GrB_UnaryOp_new (&LAGraph_IDENTITY_Complex  , complex_identity  , C, C)) ;
    OK (GrB_UnaryOp_new (&LAGraph_AINV_Complex      , complex_ainv      , C, C)) ;
    OK (GrB_UnaryOp_new (&LAGraph_ABS_Complex       , complex_abs       , C, C)) ;
    OK (GrB_UnaryOp_new (&LAGraph_MINV_Complex      , complex_minv      , C, C)) ;
    OK (GrB_UnaryOp_new (&LAGraph_NOT_Complex       , complex_not       , C, C)) ;
    OK (GrB_UnaryOp_new (&LAGraph_CONJ_Complex      , complex_conj      , C, C)) ;
    OK (GrB_UnaryOp_new (&LAGraph_ISONE_Complex     , complex_isone     , GrB_BOOL, C)) ;
    OK (GrB_UnaryOp_new (&LAGraph_TRUE_BOOL_Complex , complex_true_bool , GrB_BOOL, C)) ;

    //--------------------------------------------------------------------------
    // create the unary functions, C -> double
    //--------------------------------------------------------------------------

    OK (GrB_UnaryOp_new (&LAGraph_REAL_Complex  , complex_real  , D, C)) ;
    OK (GrB_UnaryOp_new (&LAGraph_IMAG_Complex  , complex_imag  , D, C)) ;
    OK (GrB_UnaryOp_new (&LAGraph_CABS_Complex  , complex_cabs  , D, C)) ;
    OK (GrB_UnaryOp_new (&LAGraph_ANGLE_Complex , complex_angle , D, C)) ;

    //--------------------------------------------------------------------------
    // create the unary functions, double -> C
    //--------------------------------------------------------------------------

    OK (GrB_UnaryOp_new (&LAGraph_COMPLEX_REAL_Complex , complex_complex_real , C, D)) ;
    OK (GrB_UnaryOp_new (&LAGraph_COMPLEX_IMAG_Complex , complex_complex_imag , C, D)) ;

    //--------------------------------------------------------------------------
    // create the Complex monoids
    //--------------------------------------------------------------------------

    OK (GrB_Monoid_new_UDT (&LAGraph_PLUS_Complex_MONOID,  LAGraph_PLUS_Complex,  &LAGraph_Complex_0)) ;
    OK (GrB_Monoid_new_UDT (&LAGraph_TIMES_Complex_MONOID, LAGraph_TIMES_Complex, &LAGraph_Complex_1)) ;

    //--------------------------------------------------------------------------
    // create the Complex plus-times semiring
    //--------------------------------------------------------------------------

    // more could be created, but this suffices for testing GraphBLAS
    OK (GrB_Semiring_new
        (&LAGraph_PLUS_TIMES_Complex, LAGraph_PLUS_Complex_MONOID, LAGraph_TIMES_Complex)) ;

    return (GrB_SUCCESS) ;
}


//------------------------------------------------------------------------------
// LAGraph_Complex_finalize: free all complex types, operators, monoids, and semiring
//------------------------------------------------------------------------------

GrB_Info LAGraph_Complex_finalize ( )
{

    //--------------------------------------------------------------------------
    // free the Complex plus-times semiring
    //--------------------------------------------------------------------------

    GrB_Semiring_free (&LAGraph_PLUS_TIMES_Complex) ;

    //--------------------------------------------------------------------------
    // free the Complex monoids
    //--------------------------------------------------------------------------

    GrB_Monoid_free (&LAGraph_PLUS_Complex_MONOID ) ;
    GrB_Monoid_free (&LAGraph_TIMES_Complex_MONOID) ;

    //--------------------------------------------------------------------------
    // free the Complex binary operators, CxC->C
    //--------------------------------------------------------------------------

    GrB_BinaryOp_free (&LAGraph_FIRST_Complex ) ;
    GrB_BinaryOp_free (&LAGraph_SECOND_Complex) ;
    GrB_BinaryOp_free (&LAGraph_MIN_Complex   ) ;
    GrB_BinaryOp_free (&LAGraph_MAX_Complex   ) ;
    GrB_BinaryOp_free (&LAGraph_PLUS_Complex  ) ;
    GrB_BinaryOp_free (&LAGraph_MINUS_Complex ) ;
    GrB_BinaryOp_free (&LAGraph_RMINUS_Complex) ;
    GrB_BinaryOp_free (&LAGraph_TIMES_Complex ) ;
    GrB_BinaryOp_free (&LAGraph_DIV_Complex   ) ;
    GrB_BinaryOp_free (&LAGraph_RDIV_Complex  ) ;
    GrB_BinaryOp_free (&LAGraph_PAIR_Complex  ) ;
    GrB_BinaryOp_free (&LAGraph_ANY_Complex  ) ;
    GrB_BinaryOp_free (&LAGraph_SKEW_Complex  ) ;
    GrB_BinaryOp_free (&LAGraph_HERMITIAN_Complex  ) ;

    GrB_BinaryOp_free (&LAGraph_ISEQ_Complex) ;
    GrB_BinaryOp_free (&LAGraph_ISNE_Complex) ;
    GrB_BinaryOp_free (&LAGraph_ISGT_Complex) ;
    GrB_BinaryOp_free (&LAGraph_ISLT_Complex) ;
    GrB_BinaryOp_free (&LAGraph_ISGE_Complex) ;
    GrB_BinaryOp_free (&LAGraph_ISLE_Complex) ;

    GrB_BinaryOp_free (&LAGraph_OR_Complex) ;
    GrB_BinaryOp_free (&LAGraph_AND_Complex) ;
    GrB_BinaryOp_free (&LAGraph_XOR_Complex) ;

    //--------------------------------------------------------------------------
    // free the Complex binary operators, CxC -> bool
    //--------------------------------------------------------------------------

    GrB_BinaryOp_free (&LAGraph_EQ_Complex) ;
    GrB_BinaryOp_free (&LAGraph_NE_Complex) ;
    GrB_BinaryOp_free (&LAGraph_GT_Complex) ;
    GrB_BinaryOp_free (&LAGraph_LT_Complex) ;
    GrB_BinaryOp_free (&LAGraph_GE_Complex) ;
    GrB_BinaryOp_free (&LAGraph_LE_Complex) ;

    //--------------------------------------------------------------------------
    // free the Complex binary operator, double x double -> complex
    //--------------------------------------------------------------------------

    GrB_BinaryOp_free (&LAGraph_COMPLEX_Complex) ;

    //--------------------------------------------------------------------------
    // free the Complex unary operators, C->C
    //--------------------------------------------------------------------------

    GrB_UnaryOp_free (&LAGraph_ONE_Complex       ) ;
    GrB_UnaryOp_free (&LAGraph_IDENTITY_Complex  ) ;
    GrB_UnaryOp_free (&LAGraph_AINV_Complex      ) ;
    GrB_UnaryOp_free (&LAGraph_ABS_Complex       ) ;
    GrB_UnaryOp_free (&LAGraph_MINV_Complex      ) ;
    GrB_UnaryOp_free (&LAGraph_NOT_Complex       ) ;
    GrB_UnaryOp_free (&LAGraph_CONJ_Complex      ) ;
    GrB_UnaryOp_free (&LAGraph_ISONE_Complex     ) ;
    GrB_UnaryOp_free (&LAGraph_TRUE_BOOL_Complex ) ;

    //--------------------------------------------------------------------------
    // free the unary functions, C -> double
    //--------------------------------------------------------------------------

    GrB_UnaryOp_free (&LAGraph_REAL_Complex ) ;
    GrB_UnaryOp_free (&LAGraph_IMAG_Complex ) ;
    GrB_UnaryOp_free (&LAGraph_CABS_Complex ) ;
    GrB_UnaryOp_free (&LAGraph_ANGLE_Complex) ;

    //--------------------------------------------------------------------------
    // free the unary functions, double -> C
    //--------------------------------------------------------------------------

    GrB_UnaryOp_free (&LAGraph_COMPLEX_REAL_Complex) ;
    GrB_UnaryOp_free (&LAGraph_COMPLEX_IMAG_Complex) ;

    //--------------------------------------------------------------------------
    // free the Complex type
    //--------------------------------------------------------------------------

    GrB_Type_free (&LAGraph_Complex) ;

    return (GrB_SUCCESS) ;
}

