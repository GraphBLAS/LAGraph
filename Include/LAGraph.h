//------------------------------------------------------------------------------
// LAGraph.h:  include file for user applications that use LAGraph
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// TODO: add more comments to this file.

//------------------------------------------------------------------------------
// include files and global #defines
//------------------------------------------------------------------------------

#include "GraphBLAS.h"
#include <complex.h>

#define LAGRAPH_RAND_MAX 32767

//------------------------------------------------------------------------------
// global objects
//------------------------------------------------------------------------------

// LAGraph_Complex is a GrB_Type containing the ANSI C11 double complex
// type.  This is required so that any arbitrary Matrix Market format
// can be read into GraphBLAS.
extern GrB_Type LAGraph_Complex ;

// binary operators to test for symmetry, skew-symmetry and Hermitian property
extern GrB_BinaryOp LAGraph_EQ_Complex ;
extern GrB_BinaryOp LAGraph_SKEW_INT8 ;
extern GrB_BinaryOp LAGraph_SKEW_INT16 ;
extern GrB_BinaryOp LAGraph_SKEW_INT32 ;
extern GrB_BinaryOp LAGraph_SKEW_INT64 ;
extern GrB_BinaryOp LAGraph_SKEW_FP32 ;
extern GrB_BinaryOp LAGraph_SKEW_FP64 ;
extern GrB_BinaryOp LAGraph_SKEW_Complex ;
extern GrB_BinaryOp LAGraph_Hermitian ;

// unary operators to check if the entry is equal to 1
extern GrB_UnaryOp LAGraph_ISONE_INT8 ;
extern GrB_UnaryOp LAGraph_ISONE_INT16 ;
extern GrB_UnaryOp LAGraph_ISONE_INT32 ;
extern GrB_UnaryOp LAGraph_ISONE_INT64 ;
extern GrB_UnaryOp LAGraph_ISONE_UINT8 ;
extern GrB_UnaryOp LAGraph_ISONE_UINT16 ;
extern GrB_UnaryOp LAGraph_ISONE_UINT32 ;
extern GrB_UnaryOp LAGraph_ISONE_UINT64 ;
extern GrB_UnaryOp LAGraph_ISONE_FP32 ;
extern GrB_UnaryOp LAGraph_ISONE_FP64 ;
extern GrB_UnaryOp LAGraph_ISONE_Complex ;

// boolean monoid
extern GrB_Monoid LAGraph_LAND_MONOID ;

//------------------------------------------------------------------------------
// user-callable functions
//------------------------------------------------------------------------------

GrB_Info LAGraph_init ( ) ;         // start LAGraph

GrB_Info LAGraph_finalize ( ) ;     // end LAGraph

GrB_Info LAGraph_mmread
(
    GrB_Matrix *A,      // handle of matrix to create
    FILE *f             // file to read from, already open
) ;

GrB_info LAGraph_mmwrite
(
    GrB_Matrix A,           // matrix to write to the file
    FILE *f                 // file to write it to
    // TODO , FILE *fcomments         // optional file with extra comments
) ;

GrB_Info LAGraph_ispattern  // return GrB_SUCCESS if successful
(
    bool *result,           // true if A is all one, false otherwise
    GrB_Matrix A,
    GrB_UnaryOp userop      // for A with user-defined type
) ;

GrB_Info LAGraph_isequal    // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp userop     // for A and B with user-defined types.  ignored
                            // if A and B are of built-in types
) ;

GrB_Info LAGraph_isall      // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp op         // GrB_EQ_<type>, for the type of A and B,
                            // to check for equality.  Or use any desired
                            // operator.  The operator should return GrB_BOOL.
) ;

uint64_t LAGraph_rand (uint64_t *seed) ;

uint64_t LAGraph_rand64 (uint64_t *seed) ;

double LAGraph_randx (uint64_t *seed) ;

GrB_Info LAGraph_random         // create a random matrix
(
    GrB_Matrix *A,              // handle of matrix to create
    GrB_Type type,              // built-in type, or LAGraph_Complex
    GrB_Index nrows,            // number of rows
    GrB_Index ncols,            // number of columns
    GrB_Index nvals,            // number of values
    bool make_pattern,          // if true, A is a pattern
    bool make_symmetric,        // if true, A is symmetric
    bool make_skew_symmetric,   // if true, A is skew-symmetric
    bool make_hermitian,        // if trur, A is hermitian
    bool no_diagonal,           // if true, A has no entries on the diagonal
    uint64_t *seed              // random number seed; modified on return
) ;

GrB_Info LAGraph_alloc_global ( ) ;

GrB_Info LAGraph_free_global ( ) ;

GrB_Info LAGraph_malloc
(
    void **p,               // pointer to allocated block of memory
    size_t nitems,          // number of items
    size_t size_of_item     // size of each item
) ;

void LAGraph_free
(
    void **p                // *p is freed and set to NULL
) ;

