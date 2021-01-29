//------------------------------------------------------------------------------
// LAGraph_Test_ReadProblem: read in a graph from a file
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Scott Kolodziej and Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// usage:
// test_whatever < matrixfile.mtx
// test_whatever matrixfile.mtx sourcenodes.mtx

// The matrixfile may also have a grb suffix.

#define LAGRAPH_FREE_WORK           \
{                                   \
    GrB_free (&thunk) ;             \
    GrB_free (&A) ;                 \
    GrB_free (&A2) ;                \
    if (f != NULL) fclose (f) ;     \
    f = NULL ;                      \
}

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE_WORK ;             \
    LAGraph_Delete (G, NULL) ;      \
    GrB_free (SourceNodes) ;        \
}

#include "LG_internal.h"

int LAGraph_Test_ReadProblem    // returns 0 if successful, -1 if failure
(
    // output
    LAGraph_Graph *G,           // graph from the file
    GrB_Matrix *SourceNodes,    // source nodes
    // inputs
    bool make_symmetric,        // if true, always return G as undirected
    bool remove_self_edges,     // if true, remove self edges
    bool pattern,               // if true, return G->A as bool (all true)
    GrB_Type pref,              // if non-NULL, typecast G->A to this type
    bool ensure_positive,       // if true, ensure all entries are > 0
    int argc,                   // input to main test program
    char **argv,                // input to main test program
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Matrix A = NULL, A2 = NULL ;
    GxB_Scalar thunk = NULL ;
    FILE *f = NULL ;
    LG_CHECK (G == NULL, -1, "G is missing") ;
    (*G) = NULL ;
    if (SourceNodes != NULL) (*SourceNodes) = NULL ;

    //--------------------------------------------------------------------------
    // read in a matrix from a file
    //--------------------------------------------------------------------------

    double tic [2] ;
    LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;

    if (argc > 1)
    {
        // Usage:
        //      ./test_whatever matrixfile.mtx sources.mtx
        //      ./test_whatever matrixfile.grb sources.mtx

        // read in the file in Matrix Market format from the input file
        char *filename = argv [1] ;
        printf ("matrix: %s\n", filename) ;

        // find the filename extension
        size_t len = strlen (filename) ;
        char *ext = NULL ;
        for (int k = len-1 ; k >= 0 ; k--)
        {
            if (filename [k] == '.')
            {
                ext = filename + k ;
                printf ("[%s]\n", ext) ;
                break ;
            }
        }
        bool is_binary = (ext != NULL && strncmp (ext, ".grb", 4) == 0) ;

        if (is_binary)
        {
            printf ("Reading binary file: %s\n", filename) ;
            LAGraph_TRY (LAGraph_BinRead (&A, filename, msg)) ;
        }
        else
        {
            printf ("Reading Matrix Market file: %s\n", filename) ;
            f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Matrix file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            LAGraph_TRY (LAGraph_MMRead (&A, f, msg)) ;
            fclose (f) ;
            f = NULL ;
        }

        // read in source nodes in Matrix Market format from the input file
        if (argc > 2 && SourceNodes != NULL)
        {
            // do not read in the file if the name starts with "-"
            filename = argv [2] ;
            if (filename [0] != '-')
            {
                printf ("sources: %s\n", filename) ;
                f = fopen (filename, "r") ;
                if (f == NULL)
                {
                    printf ("Source node file not found: [%s]\n", filename) ;
                    exit (1) ;
                }
                LAGraph_TRY (LAGraph_MMRead (SourceNodes, f, msg)) ;
                fclose (f) ;
                f = NULL ;
            }
        }
    }
    else
    {

        // Usage:  ./test_whatever < matrixfile.mtx
        printf ("matrix: from stdin\n") ;

        // read in the file in Matrix Market format from stdin
        LAGraph_TRY (LAGraph_MMRead (&A, stdin, msg)) ;
    }

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols ;
    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_Index n = nrows ;
    LG_CHECK (nrows != ncols, -1, "A must be square") ;

    //--------------------------------------------------------------------------
    // typecast, if requested
    //--------------------------------------------------------------------------

    GrB_Type type ;
    GrB_TRY (GxB_Matrix_type (&type, A)) ;
    if (pattern)
    {
        // convert to boolean, pattern-only, with all entries true
        LAGraph_TRY (LAGraph_Pattern (&A2, A, msg)) ;
    }
    else if (pref != NULL && type != pref)
    {
        // convert to the requested type
        GrB_TRY (GrB_Matrix_new (&A2, pref, n, n)) ;

        GrB_UnaryOp op = NULL ;
        if      (pref == GrB_BOOL  ) op = GrB_IDENTITY_BOOL ;
        else if (pref == GrB_INT8  ) op = GrB_IDENTITY_INT8 ;
        else if (pref == GrB_INT16 ) op = GrB_IDENTITY_INT16 ;
        else if (pref == GrB_INT32 ) op = GrB_IDENTITY_INT32 ;
        else if (pref == GrB_INT64 ) op = GrB_IDENTITY_INT64 ;
        else if (pref == GrB_UINT8 ) op = GrB_IDENTITY_UINT8 ;
        else if (pref == GrB_UINT16) op = GrB_IDENTITY_UINT16 ;
        else if (pref == GrB_UINT32) op = GrB_IDENTITY_UINT32 ;
        else if (pref == GrB_UINT64) op = GrB_IDENTITY_UINT64 ;
        else if (pref == GrB_FP32  ) op = GrB_IDENTITY_FP32 ;
        else if (pref == GrB_FP64  ) op = GrB_IDENTITY_FP64 ;
        else if (pref == GxB_FC32  ) op = GxB_IDENTITY_FC32 ;
        else if (pref == GxB_FC64  ) op = GxB_IDENTITY_FC64 ;

        GrB_TRY (GrB_apply (A2, NULL, NULL, op, A, NULL)) ;
    }

    if (A2 != NULL)
    {
        GrB_free (&A) ;
        A = A2 ;
        A2 = NULL ;
        GrB_TRY (GrB_wait (&A)) ;
    }

    //--------------------------------------------------------------------------
    // remove self-edges, if requested
    //--------------------------------------------------------------------------

    if (remove_self_edges)
    {
        GrB_TRY (GxB_Scalar_new (&thunk, GrB_INT64)) ;
        GrB_TRY (GxB_Scalar_setElement (thunk, 0)) ;
        GrB_TRY (GxB_select (A, NULL, NULL, GxB_OFFDIAG, A, thunk, NULL)) ;
        GrB_free (&thunk) ;
    }

    //--------------------------------------------------------------------------
    // ensure all entries are > 0, if requested
    //--------------------------------------------------------------------------

    if (!pattern && ensure_positive)
    {
        // drop explicit zeros
        GrB_TRY (GxB_select (A, NULL, NULL, GxB_NONZERO, A, NULL, NULL)) ;

        // A = abs (A)
        GrB_UnaryOp op = NULL ;
             if (type == GrB_INT8  ) op = GrB_ABS_INT8 ;
        else if (type == GrB_INT16 ) op = GrB_ABS_INT16 ;
        else if (type == GrB_INT32 ) op = GrB_ABS_INT32 ;
        else if (type == GrB_INT64 ) op = GrB_ABS_INT64 ;
        else if (type == GrB_FP32  ) op = GrB_ABS_FP32 ;
        else if (type == GrB_FP64  ) op = GrB_ABS_FP64 ;
        else if (type == GxB_FC32  ) op = GxB_ABS_FC32 ;
        else if (type == GxB_FC64  ) op = GxB_ABS_FC64 ;
        if (op != NULL)
        {
            GrB_TRY (GrB_apply (A, NULL, NULL, op, A, NULL)) ;
        }
    }

    //--------------------------------------------------------------------------
    // construct the graph
    //--------------------------------------------------------------------------

    bool A_is_symmetric =
        (n == 134217726 ||  // HACK for kron
         n == 134217728) ;  // HACK for urand

    if (A_is_symmetric)
    {
        // A is known to be symmetric
        LAGraph_TRY (LAGraph_New (G, &A, LAGRAPH_ADJACENCY_UNDIRECTED, msg)) ;
        ASSERT ((*G)->A_pattern_is_symmetric == true) ;
    }
    else
    {
        // compute G->AT and determine if A has a symmetric pattern
        LAGraph_TRY (LAGraph_New (G, &A, LAGRAPH_ADJACENCY_DIRECTED, msg)) ;
        LAGraph_TRY (LAGraph_Property_ASymmetricPattern (*G, msg)) ;
        if ((*G)->A_pattern_is_symmetric && pattern)
        {
            // if G->A has a symmetric pattern, declare the graph undirected
            // and free G->AT since it isn't needed.
            (*G)->kind = LAGRAPH_ADJACENCY_UNDIRECTED ;
            GrB_TRY (GrB_Matrix_free (&((*G)->AT))) ;
        }
        else if (make_symmetric)
        {
            // make sure G->A is symmetric
            bool sym ;
            LAGraph_TRY (LAGraph_IsEqual (&sym, (*G)->A, (*G)->AT, NULL, msg)) ;
            if (!sym)
            {
                GrB_BinaryOp op = NULL ;
                GrB_Type type ;
                GrB_TRY (GxB_Matrix_type (&type, (*G)->A)) ;
                if      (type == GrB_BOOL  ) op = GrB_LOR ;
                else if (type == GrB_INT8  ) op = GrB_PLUS_INT8 ;
                else if (type == GrB_INT16 ) op = GrB_PLUS_INT16 ;
                else if (type == GrB_INT32 ) op = GrB_PLUS_INT32 ;
                else if (type == GrB_INT64 ) op = GrB_PLUS_INT64 ;
                else if (type == GrB_UINT8 ) op = GrB_PLUS_UINT8 ;
                else if (type == GrB_UINT16) op = GrB_PLUS_UINT16 ;
                else if (type == GrB_UINT32) op = GrB_PLUS_UINT32 ;
                else if (type == GrB_UINT64) op = GrB_PLUS_UINT64 ;
                else if (type == GrB_FP32  ) op = GrB_PLUS_FP32 ;
                else if (type == GrB_FP64  ) op = GrB_PLUS_FP64 ;
                else if (type == GxB_FC32  ) op = GxB_PLUS_FC32 ;
                else if (type == GxB_FC64  ) op = GxB_PLUS_FC64 ;
                GrB_TRY (GrB_eWiseAdd ((*G)->A, NULL, NULL, op,
                    (*G)->A, (*G)->AT, NULL)) ;
                GrB_TRY (GrB_Matrix_free (&((*G)->AT))) ;
            }
            (*G)->kind = LAGRAPH_ADJACENCY_UNDIRECTED ;
            (*G)->A_pattern_is_symmetric = true ;
        }
    }

    (*G)->ndiag = (remove_self_edges) ? 0 : LAGRAPH_UNKNOWN ;

    //--------------------------------------------------------------------------
    // generate 64 random source nodes, if requested but not provided on input
    //--------------------------------------------------------------------------

    #define NSOURCES 64

    if (SourceNodes != NULL && (*SourceNodes == NULL))
    {
        GrB_TRY (GrB_Matrix_new (SourceNodes, GrB_INT64, NSOURCES, 1)) ;
        srand (1) ;
        for (int k = 0 ; k < NSOURCES ; k++)
        {
            int64_t i = 1 + (rand ( ) % n) ;    // in range 1 to n
            // SourceNodes [k] = i
            GrB_TRY (GrB_Matrix_setElement (*SourceNodes, i, k, 0)) ;
        }
    }

    if (SourceNodes != NULL)
    {
        GrB_TRY (GrB_wait (SourceNodes)) ;
    }

    //--------------------------------------------------------------------------
    // free workspace, print a summary of the graph, and return result
    //--------------------------------------------------------------------------

    double t_read ;
    LAGraph_TRY (LAGraph_Toc (&t_read, tic, msg)) ;
    printf ("read time: %g\n", t_read) ;

    LAGRAPH_FREE_WORK ;
    LAGraph_TRY (LAGraph_DisplayGraph (*G, 0, msg)) ;
    return (0) ;
}

