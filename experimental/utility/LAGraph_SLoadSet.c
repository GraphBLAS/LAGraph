//------------------------------------------------------------------------------
// LAGraph_SLoadSet: load a set of matrices from a *.lagraph file
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// Contributed by Tim Davis, Texas A&M University

//------------------------------------------------------------------------------

// LAgraph_SLoadSet loads a set of GrB_Matrix objects from a *.lagraph file.
// It returns a GrB_Matrix array of size nmatrices.  In the future, it will
// also return a set of GrB_Vectors, and a an array of uncompressed ascii
// texts.  The caller is responsible for freeing the output of this method,
// via:

//      LAGraph_Free ((void **) &collection) ;
//      LAGraph_SFreeSet (&Set, nmatrices) ;

// See also LAGraph_SRead, which just reads in the serialized objects and
// does not convert them to their corresponding GrB_Matrix, GrB_Vector, or
// uncompressed texts.

//------------------------------------------------------------------------------

#define LAGraph_FREE_WORK                                           \
{                                                                   \
    if (f != NULL) fclose (f) ;                                     \
    f = NULL ;                                                      \
    LAGraph_SFreeContents (&Contents, ncontents) ;                  \
}

#define LAGraph_FREE_ALL                                            \
{                                                                   \
    LAGraph_FREE_WORK ;                                             \
    LAGraph_SFreeSet (&Set, nmatrices) ;                            \
    LAGraph_Free ((void **) &collection) ;                          \
}

#include "LG_internal.h"
#include "LAGraphX.h"

//------------------------------------------------------------------------------
// LAGraph_SLoadSet
//------------------------------------------------------------------------------

int LAGraph_SLoadSet            // load a set of matrices from a *.lagraph file
(
    // input:
    char *filename,                 // name of file to read; NULL for stdin
    // outputs:
    GrB_Matrix **Set_handle,        // array of GrB_Matrix of size nmatrices
    GrB_Index *nmatrices_handle,    // # of matrices loaded from *.lagraph file
//  TODO:
//  GrB_Vector **Set_handle,        // array of GrB_Vector of size nvector
//  GrB_Index **nvectors_handle,    // # of vectors loaded from *.lagraph file
//  char **Text_handle,             // array of pointers to (char *) strings
//  GrB_Index **ntext_handle,       // # of texts loaded from *.lagraph file
    char **collection_handle,       // name of this collection of matrices
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    FILE *f = NULL ;
    char *collection = NULL ;
    GrB_Matrix *Set = NULL ;
    LAGraph_Contents *Contents = NULL ;
    GrB_Index ncontents = 0 ;
    GrB_Index nmatrices = 0 ;
//  GrB_Index nvectors = 0 ;
//  GrB_Index ntexts = 0 ;

    LG_CHECK (Set_handle == NULL || nmatrices_handle == NULL
        || collection_handle == NULL, GrB_NULL_POINTER, "inputs are NULL") ;

    //--------------------------------------------------------------------------
    // read the file
    //--------------------------------------------------------------------------

    if (filename == NULL)
    {
        LAGraph_TRY (LAGraph_SRead (stdin, &collection, &Contents, &ncontents,
            msg)) ;
    }
    else
    {
        f = fopen (filename, "r") ;
        LG_CHECK (f == NULL, -1002, "unable to open input file") ;
        LAGraph_TRY (LAGraph_SRead (f, &collection, &Contents, &ncontents,
            msg)) ;
        fclose (f) ;
        f = NULL ;
    }

    //--------------------------------------------------------------------------
    // count the matrices/vectors/texts in the Contents
    //--------------------------------------------------------------------------

    // TODO: for now, all Contents are matrices
    nmatrices = ncontents ;

#if 0
    for (GrB_Index i = 0 ; i < ncontents ; i++)
    {
        switch (Contents [i].kind)
        {
            case LAGraph_matrix_kind : nmatrices++ ; break ;
            case LAGraph_vector_kind : nvectors++  ; break ;
            case LAGraph_text_kind   : ntexts++    ; break ;
            default : LG_CHECK (true, GrB_INVALID_VALUE, "unknown kind") ;
        }
    }
    if (nvectors > 0 || ntexts > 0)
    {
        // TODO
        printf ("Warning: %lu vectors and %lu texts ignored\n",
            nvectors, ntexts) ;
    }
#endif

    //--------------------------------------------------------------------------
    // convert all the matrices (skip vectors and text content for now)
    //--------------------------------------------------------------------------

    Set = LAGraph_Calloc (nmatrices, sizeof (GrB_Matrix)) ;
    LG_CHECK (Set == NULL, GrB_OUT_OF_MEMORY, "out of memory") ;

    GrB_Index kmatrices = 0 ;
    for (GrB_Index i = 0 ; i < ncontents ; i++)
    {
        // convert Contents [i]
        void *blob = Contents [i].blob ;
        size_t blob_size = Contents [i].blob_size ;

        if (Contents [i].kind == LAGraph_matrix_kind)
        {
            // convert Contents [i].typename to a GrB_Type ctype.
            // SuiteSparse:GraphBLAS allows this to be NULL for built-in types.
            GrB_Type ctype = NULL ;
            #if LG_SUITESPARSE
            // TODO For user-defined types, LAGraph would need to pass in an
            // array of registered user-defined types.  If GxB_Type_from_name
            // returns NULL, it would then look through that list of types.
            GrB_TRY (GxB_Type_from_name (&ctype, Contents [i].type_name)) ;
            #endif
            GrB_TRY (GrB_Matrix_deserialize (&(Set [kmatrices]), ctype, blob,
                blob_size)) ;
            kmatrices++ ;
        }
        // TODO: handle vectors and texts
        // else if (Content [i].kind == LAGraph_vector_kind) ...
        // else if (Content [i].kind == LAGraph_text_kind) ...

        // free the ith blob
        LAGraph_Free ((void **) &(Contents [i].blob)) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGraph_FREE_WORK ;
    (*Set_handle) = Set ;
    (*collection_handle) = collection ;
    (*nmatrices_handle) = nmatrices ;
    return (0) ;
}

