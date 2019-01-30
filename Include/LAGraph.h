//------------------------------------------------------------------------------
// LAGraph.h: 
//------------------------------------------------------------------------------

#include "GraphBLAS.h"
#include <complex.h>

// LAGraph_Complex is a GrB_Type containing the ANSI C11 double complex
// type.  This is required so that any arbitrary Matrix Market format
// can be read into GraphBLAS.
extern GrB_Type LAGraph_Complex ;

GrB_Info LAGraph_init ( ) ;         // start LAGraph
GrB_Info LAGraph_finalize ( ) ;     // end LAGraph

GrB_Info LAGraph_mmread
(
    GrB_Matrix *A,      // handle of matrix to create
    FILE *f             // file to read from, already open
) ;

GrB_Info LAGraph_isequal_type   // return GrB_SUCCESS if successful
(
    bool *result,               // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp op             // GrB_EQ_<type>, for the type of A and B
) ;

GrB_Info LAGraph_isequal    // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp userop     // for A and B with user-defined types.  ignored
                            // if A and B are of built-in types
) ;

