//------------------------------------------------------------------------------
// gtest: read in a graph from a binary file
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

    
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M.

// Usage:

//      ./build/gtest filename edgetype

// See LAGraph/Source/Utility/LAGraph_grread.c for a description of the
// file format.

// The filename is required.  See t1.gr and t2.gr as examples.  The edgetype
// is optional.  If not present, then the expected edge weight size is zero
// (the graph has no edge weights), and the matrix G is read in as GrB_BOOL
// with all edges with weight of 1.

// Otherwise, the following edgetypes may be used:

//      edgetype    G type      edgesize in the file must be:
//      --------    ------      -----------------------------
//      bool        GrB_BOOL    1
//      int8        GrB_INT8    1
//      int16       GrB_INT16   2
//      int32       GrB_INT32   4
//      int64       GrB_INT64   8
//      uint8       GrB_UINT8   1
//      uint16      GrB_UINT16  2
//      uint32      GrB_UINT32  4
//      uint64      GrB_UINT64  8
//      float       GrB_FP32    4
//      double      GrB_FP64    8

//------------------------------------------------------------------------------

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&G) ;                 \
}

//------------------------------------------------------------------------------
// gtest main:
//------------------------------------------------------------------------------

int main (int argc, char **argv)
{

    GrB_Info info ;
    GrB_Matrix G = NULL ;
    uint64_t G_version = 0 ;

    printf ("\n\n================= Gr_Read/gtest: test LAGraph_grread\n") ;
    LAGRAPH_OK (LAGraph_init ( )) ;

    if (argc < 2 || argc > 3)
    {
        LAGRAPH_ERROR ("usage: gtest grbinaryfilename.gr edgetype",
            GrB_INVALID_VALUE) ;
    }

    printf ("filename: %s\n", argv [1]) ;

    // TODO: this is ugly; can't the edge type be determined from the file
    // itself?  Note that an esize of 4 in the header of the file is ambigous.
    // Is it int32? uint32? float?  All of them have sizeof(...) = 4.
    // So the *caller* of LAGraph_grread must state the type of the graph.

    // determine the edgetype
    GrB_Type gtype ;
    if (argc < 3)
    {
        // the file has no edge weights
        gtype = NULL ;
    }
    else if (strcmp (argv [2], "bool") == 0)
    {
        gtype = GrB_BOOL ;
    }
    else if (strcmp (argv [2], "int8") == 0)
    {
        gtype = GrB_INT8 ;
    }
    else if (strcmp (argv [2], "int16") == 0)
    {
        gtype = GrB_INT16 ;
    }
    else if (strcmp (argv [2], "int32") == 0)
    {
        gtype = GrB_INT32 ;
    }
    else if (strcmp (argv [2], "int64") == 0)
    {
        gtype = GrB_INT64 ;
    }
    else if (strcmp (argv [2], "uint8") == 0)
    {
        gtype = GrB_UINT8 ;
    }
    else if (strcmp (argv [2], "uint16") == 0)
    {
        gtype = GrB_UINT16 ;
    }
    else if (strcmp (argv [2], "uint32") == 0)
    {
        gtype = GrB_UINT32 ;
    }
    else if (strcmp (argv [2], "uint64") == 0)
    {
        gtype = GrB_UINT64 ;
    }
    else if (strcmp (argv [2], "float") == 0)
    {
        gtype = GrB_FP32 ;
    }
    else if (strcmp (argv [2], "double") == 0)
    {
        gtype = GrB_FP64 ;
    }
    else
    {
        LAGRAPH_ERROR ("unknown type", GrB_INVALID_VALUE) ;
    }

    if (gtype == NULL)
    {
        printf ("Graph is unweighted; G will be GrB_BOOL\n") ;
    }
    else
    {
        printf ("Graph is weighted, with the type: %s\n", argv [2]) ;
        LAGRAPH_OK (GxB_print (gtype, GxB_COMPLETE)) ;
    }

    // read the graph
    LAGRAPH_OK (LAGraph_grread (&G, &G_version, argv [1], gtype)) ;

    // TODO what is this?
    printf ("G_version: %g\n", (double) G_version) ;

    // print and check the graph
    LAGRAPH_OK (GxB_print (G, GxB_SHORT)) ;

    printf ("gtest: all tests passed\n\n") ;
    GrB_free (&G) ;
    LAGraph_finalize ( ) ;
    return (0) ;
}

