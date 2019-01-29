//------------------------------------------------------------------------------
// LAGraph.h: 
//------------------------------------------------------------------------------

#include "GraphBLAS.h"
#include <complex.h>

// LAGraph_Complex is a GrB_Type containing the ANSI C11 double complex
// type.  This is required so that any arbitrary Matrix Market format
// can be read into GraphBLAS.
extern GrB_Type LAGraph_Complex ;

GrB_Info LAGraph_mmread
(
    GrB_Matrix *A,      // handle of matrix to create
    FILE *f             // file to read from, already open
) ;

GrB_Info LAGraph_init ( ) ;         // start LAGraph
GrB_Info LAGraph_finalize ( ) ;     // end LAGraph
