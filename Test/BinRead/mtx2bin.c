//------------------------------------------------------------------------------
// mtx2bin: convert Matrix Market file to SuiteSparse:GraphBLAS binary file
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
// mtx2bin infile.mtx outfile.grb

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&A) ;                 \
}

int main (int argc, char **argv)
{
    GrB_Info info ;
    GrB_Matrix A = NULL ;

    if (argc < 3)
    {
        LAGRAPH_ERROR ("Usage: mxt2bin infile.mtx outfile.grb",
            GrB_INVALID_VALUE) ;
    }

    printf ("infile:  %s\n", argv [1]) ;
    printf ("outfile: %s\n", argv [2]) ;

    LAGRAPH_OK (LAGraph_init ( )) ;

    //--------------------------------------------------------------------------
    // read matrix from input file
    //--------------------------------------------------------------------------

    double tic [2] ;
    LAGraph_tic (tic) ;

    // read in the file in Matrix Market format from the input file
    FILE *f = fopen (argv [1], "r") ;
    if (f == NULL)
    {
        printf ("Matrix file not found: [%s]\n", argv [1]) ;
        exit (1) ;
    }
    LAGRAPH_OK (LAGraph_mmread(&A, f)) ;
    fclose (f) ;

    GrB_Index nvals ;
    LAGRAPH_OK (GrB_Matrix_nvals (&nvals, A)) ;
    LAGRAPH_OK (GxB_fprint (A, 2, stdout)) ;

    double t_read = LAGraph_toc (tic) ;
    printf ("read time: %g sec\n", t_read) ;

    //--------------------------------------------------------------------------
    // write to output file
    //--------------------------------------------------------------------------

    LAGraph_tic (tic) ;
    LAGRAPH_OK (LAGraph_binwrite (&A, argv [2], argv [1])) ;

    double t_binwrite = LAGraph_toc (tic) ;
    printf ("binary write time: %g sec\n", t_binwrite) ;

    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS) ;
}
