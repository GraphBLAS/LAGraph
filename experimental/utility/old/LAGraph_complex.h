
#if 0

//------------------------------------------------------------------------------
// Complex type, scalars, monoids, and semiring
//------------------------------------------------------------------------------

// See:
// https://www.drdobbs.com/complex-arithmetic-in-the-intersection-o/184401628#

#if defined ( __cplusplus )

    extern "C++"
    {
        // C++ complex types
        #include <cmath>
        #include <complex>
        #undef I
        typedef std::complex<float>  LAGraph_FC32_t ;
        typedef std::complex<double> LAGraph_FC64_t ;
    }

    #define LAGraph_CMPLXF(r,i) LAGraph_FC32_t(r,i)
    #define LAGraph_CMPLX(r,i)  LAGraph_FC64_t(r,i)

#elif ( _MSC_VER && !__INTEL_COMPILER )

    // Microsoft Windows complex types
    #include <complex.h>
    #undef I
    typedef _Fcomplex LAGraph_FC32_t ;
    typedef _Dcomplex LAGraph_FC64_t ;

    #define LAGraph_CMPLXF(r,i) (_FCbuild (r,i))
    #define LAGraph_CMPLX(r,i)  ( _Cbuild (r,i))

#else

    // ANSI C11 complex types
    #include <complex.h>
    #undef I
    typedef float  complex LAGraph_FC32_t ;
    typedef double complex LAGraph_FC64_t ;

    #ifndef CMPLX
        // gcc 6.2 on the the Mac doesn't #define CMPLX
        #define LAGraph_CMPLX(r,i) \
        ((LAGraph_FC64_t)((double)(r)) + \
         (LAGraph_FC64_t)((double)(i) * _Complex_I))
    #else
        // use the ANSI C11 CMPLX macro
        #define LAGraph_CMPLX(r,i) CMPLX (r,i)
    #endif

    #ifndef CMPLXF
        // gcc 6.2 on the the Mac doesn't #define CMPLXF
        #define LAGraph_CMPLXF(r,i) \
        ((LAGraph_FC32_t)((float)(r)) + \
         (LAGraph_FC32_t)((float)(i) * _Complex_I))
    #else
        // use the ANSI C11 CMPLXF macro
        #define LAGraph_CMPLXF(r,i) CMPLXF (r,i)
    #endif

#endif

extern GrB_Type LAGraph_ComplexFP64 ;

extern GrB_Monoid
    LAGraph_PLUS_ComplexFP64_MONOID       ,
    LAGraph_TIMES_ComplexFP64_MONOID      ;

extern GrB_Semiring LAGraph_PLUS_TIMES_ComplexFP64 ;

extern LAGraph_FC64_t LAGraph_ComplexFP64_1 ;
extern LAGraph_FC64_t LAGraph_ComplexFP64_0 ;

GrB_Info LAGraph_Complex_init ( ) ;
GrB_Info LAGraph_Complex_finalize ( ) ;

#endif
