//------------------------------------------------------------------------------
// binread: read a SuiteSparse:GraphBLAS binary file
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

// Contributed by Tim Davis, Texas A&M University

// usage:
// binread infile.grb

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&A) ;                 \
}

int main (int argc, char **argv)
{
    GrB_Info info ;
    GrB_Matrix A = NULL ;

    if (argc < 2)
    {
        LAGRAPH_ERROR ("Usage: binread infile.grb",
            GrB_INVALID_VALUE) ;
    }

    printf ("infile:  %s\n", argv [1]) ;

    LAGRAPH_OK (LAGraph_init ( )) ;

    //--------------------------------------------------------------------------
    // read matrix from input file
    //--------------------------------------------------------------------------

    double tic [2] ;
    LAGraph_tic (tic) ;

    // read in the file
    LAGRAPH_OK (LAGraph_binread (&A, argv [1])) ;

    double t_read = LAGraph_toc (tic) ;
    printf ("read time: %g sec\n", t_read) ;

    LAGRAPH_OK (GxB_fprint (A, 2, stdout)) ;

    //--------------------------------------------------------------------------
    // free everthing
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    GrB_finalize ( ) ;
    return (GrB_SUCCESS) ;
}
