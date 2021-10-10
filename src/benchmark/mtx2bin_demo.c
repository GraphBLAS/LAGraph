//------------------------------------------------------------------------------
// mtx2bin: convert Matrix Market file to SuiteSparse:GraphBLAS binary file
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M University

// usage:
// mtx2bin infile.mtx outfile.grb

#include "LAGraph_demo.h"

#define LAGraph_FREE_ALL            \
{                                   \
    GrB_free (&A) ;                 \
}

int main (int argc, char **argv)
{
    GrB_Info info ;
    GrB_Matrix A = NULL ;
    GrB_Type atype = NULL ;
    char msg [LAGRAPH_MSG_LEN] ;

    if (argc < 3)
    {
        printf ("Usage: mxt2bin infile.mtx outfile.grb\n") ;
        exit (1) ;
    }

    printf ("infile:  %s\n", argv [1]) ;
    printf ("outfile: %s\n", argv [2]) ;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read matrix from input file
    //--------------------------------------------------------------------------

    double tic [2] ;
    LAGraph_TRY (LAGraph_Tic (tic, msg)) ;

    // read in the file in Matrix Market format from the input file
    FILE *f = fopen (argv [1], "r") ;
    if (f == NULL)
    {
        printf ("Matrix file not found: [%s]\n", argv [1]) ;
        exit (1) ;
    }
    LAGraph_TRY (LAGraph_MMRead (&A, &atype, f, msg)) ;
    fclose (f) ;

    LAGraph_TRY (LAGraph_Matrix_wait (A, msg)) ;

    double t_read ;
    LAGraph_TRY (LAGraph_Toc (&t_read, tic, msg)) ;
    printf ("read time: %g sec\n", t_read) ;

    //--------------------------------------------------------------------------
    // write to output file
    //--------------------------------------------------------------------------

    LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
    f = fopen (argv [2], "w") ;
    if (f == NULL)
    {
        printf ("Unable to open binary output file: [%s]\n", argv [2]) ;
        exit (1) ;
    }
    if (binwrite (&A, f, argv [1]) != 0)
    {
        printf ("Unable to create binary file\n") ;
        exit (1) ;
    }
    double t_binwrite ;
    LAGraph_TRY (LAGraph_Toc (&t_binwrite, tic, msg)) ;
    printf ("binary write time: %g sec\n", t_binwrite) ;

    LAGraph_FREE_ALL ;
    return (0) ;
}

