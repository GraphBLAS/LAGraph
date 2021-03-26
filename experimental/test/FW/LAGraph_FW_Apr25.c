//------------------------------------------------------------------------------
// LAGraph_FW: Floyd-Warshall method: all pairs shortest paths
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

// The input is a square unsymmetric matrix G, for a directed graph.  G can be
// of any type.  If it is real (float or double) or a 64-bit integer, then the
// output is of type GrB_FP64.  Otherwise, the output is of type GrB_INT32.

// G(i,j) is the edge weight for edge (i,j).  D(i,j) on output is the length of
// the shortest path from node i to j, if the entry is present.  If D(i,j) is
// not present then there is no path from i to j.  The shortest path itself
// is not returned.

// Negative weights are OK, unless there is a negative weight cycle.  In
// that case, the output is undefined.

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include "LAGraph.h"

#define LAGRAPH_FREE_ALL                                                   \
    GrB_free (D);                                                          \
    GrB_free (&A);                                                         \
    GrB_free (&B);

GrB_Info LAGraph_FW
(
    const GrB_Matrix G,     // input graph, with edge weights
    GrB_Matrix *D           // output graph, created on output
)
{
    GrB_Info info;
    GrB_Matrix A = NULL, B = NULL ;

    GxB_print (G, 1) ;

    // make sure D is valid
    if (D == NULL)
    {
        return (GrB_INVALID_VALUE) ;
    }
    (*D) = NULL ;

    // determine the type of the input and output graphs
    GrB_Type gtype, otype ;
    GrB_BinaryOp op ;
    GrB_Semiring semiring ;
    LAGRAPH_OK (GxB_Matrix_type (&gtype, G)) ;
    if (gtype == GrB_FP64 || gtype == GrB_FP32 ||
        gtype == GrB_INT64 || gtype == GrB_UINT64)
    {
        otype = GrB_FP64 ;
        semiring = GxB_MIN_PLUS_FP64 ;
        op = GrB_MIN_FP64 ;
    }
    else
    {
        otype = GrB_INT32 ;
        semiring = GxB_MIN_PLUS_INT32 ;
        op = GrB_MIN_INT32 ;
    }

    GrB_Index n, ncols;
    LAGRAPH_OK (GrB_Matrix_nrows(&n, G));
    LAGRAPH_OK (GrB_Matrix_ncols(&ncols, G));
    if (n != ncols)
    {
        return (GrB_INVALID_VALUE) ;
    }

    // printf ("here\n") ;

    LAGRAPH_OK (GrB_Matrix_new (D,  otype, n, n)) ;
    LAGRAPH_OK (GrB_Matrix_new (&A, otype, n, 1));
    LAGRAPH_OK (GrB_Matrix_new (&B, otype, 1, n));

    // D = G, with possible typecasting
    LAGRAPH_OK (GrB_assign (*D, GrB_NULL, GrB_NULL, G, GrB_ALL, n, GrB_ALL, n, GrB_NULL)) ;

    GxB_print (*D, 1) ;

    for (GrB_Index i = 0; i < n; i++)
    {
        // printf ("i %g\n", (double) i) ;
        // A = D (:,i), the ith column
        LAGRAPH_OK (GrB_extract(A, GrB_NULL, GrB_NULL, *D, GrB_ALL, n, &i, 1, GrB_NULL));
        // B = D (i,:), the ith row
        LAGRAPH_OK (GrB_extract(B, GrB_NULL, GrB_NULL, *D, &i, 1, GrB_ALL, n, GrB_NULL));
        // D = min (D,A*B) with "*" being the min-plus semiring
        LAGRAPH_OK (GrB_mxm(*D, GrB_NULL, op, semiring, A, B, GrB_NULL));
    }

    return GrB_SUCCESS;
}

