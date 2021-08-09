//------------------------------------------------------------------------------
// LAGraph_tsvread: read a tsv file
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// LAGraph_tsvread: read a tsv file.  Contributed by Tim Davis, Texas A&M
// University.

// Reads a tsv file.  Each line in the file specifies a single entry: i, j, x.
// The indices i and j are assumed to be one-based.  The dimensions of the
// matrix must be provided by the caller.  This format is used for matrices at
// http://graphchallenge.org.  The Matrix Market format is recommended instead;
// it is more flexible and easier to use, since that format includes the matrix
// type and size in the file itself.  See LAGraph_mmread and LAGraph_mmwrite.

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING

#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL GrB_free (Chandle) ;

GrB_Info LAGraph_tsvread        // returns GrB_SUCCESS if successful
(
    GrB_Matrix *Chandle,        // C, created on output
    FILE *f,                    // file to read from (already open)
    GrB_Type type,              // the type of C to create
    GrB_Index nrows,            // C is nrows-by-ncols
    GrB_Index ncols
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (Chandle == NULL || f == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    //--------------------------------------------------------------------------
    // create the output matrix
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Matrix C = NULL ;
    (*Chandle) = NULL ;
    LAGRAPH_OK (GrB_Matrix_new (&C, type, nrows, ncols)) ;

    //--------------------------------------------------------------------------
    // read the entries
    //--------------------------------------------------------------------------

    GrB_Index i, j ;

    if (type == GrB_INT64)
    {

        //----------------------------------------------------------------------
        // read the entries as int64
        //----------------------------------------------------------------------

        int64_t x ;
        while (fscanf (f, "%"PRIu64"%"PRIu64"%"PRId64"\n", &i, &j, &x) != EOF)
        {
            LAGRAPH_OK (GrB_Matrix_setElement (C, x, i-1, j-1)) ;
        }

    }
    else if (type == GrB_UINT64)
    {

        //----------------------------------------------------------------------
        // read the entries as uint64
        //----------------------------------------------------------------------

        uint64_t x ;
        while (fscanf (f, "%"PRIu64"%"PRIu64"%"PRIu64"\n", &i, &j, &x) != EOF)
        {
            LAGRAPH_OK (GrB_Matrix_setElement (C, x, i-1, j-1)) ;
        }

    }
    else
    {

        //----------------------------------------------------------------------
        // read the entries as double, and typecast to the matrix type
        //----------------------------------------------------------------------

        double x ;
        while (fscanf (f, "%"PRIu64"%"PRIu64"%lg\n", &i, &j, &x) != EOF)
        {
            LAGRAPH_OK (GrB_Matrix_setElement (C, x, i-1, j-1)) ;
        }
    }

    //--------------------------------------------------------------------------
    // finalize the matrix and return the result
    //--------------------------------------------------------------------------

    GrB_Index ignore ;
    LAGRAPH_OK (GrB_Matrix_nvals (&ignore, C)) ;
    (*Chandle) = C ;
    return (GrB_SUCCESS) ;
}
