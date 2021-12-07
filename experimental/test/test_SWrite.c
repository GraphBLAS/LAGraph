//----------------------------------------------------------------------------
// LAGraph/src/test/test_SWrite.c: test cases for LAGraph_SWrite and SRead
// ----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//-----------------------------------------------------------------------------

#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;
GrB_Matrix A = NULL ;
GrB_Matrix B = NULL ;
GrB_Matrix *S = NULL ;

GrB_Type atype = NULL ;
#define LEN 512
char filename [LEN+1] ;

#define NFILES 51
const char *files [ ] =
{
    "A.mtx",
    "cover.mtx",
    "cover_structure.mtx",
    "jagmesh7.mtx",
    "ldbc-cdlp-directed-example.mtx",
    "ldbc-cdlp-undirected-example.mtx",
    "ldbc-directed-example-bool.mtx",
    "ldbc-directed-example.mtx",
    "ldbc-directed-example-unweighted.mtx",
    "ldbc-undirected-example-bool.mtx",
    "ldbc-undirected-example.mtx",
    "ldbc-undirected-example-unweighted.mtx",
    "ldbc-wcc-example.mtx",
    "LFAT5.mtx",
    "msf1.mtx",
    "msf2.mtx",
    "msf3.mtx",
    "sample2.mtx",
    "sample.mtx",
    "sources_7.mtx",
    "olm1000.mtx",
    "bcsstk13.mtx",
    "cryg2500.mtx",
    "tree-example.mtx",
    "west0067.mtx",
    "lp_afiro.mtx",
    "lp_afiro_structure.mtx",
    "karate.mtx",
    "matrix_bool.mtx",
    "matrix_int8.mtx",
    "matrix_int16.mtx",
    "matrix_int32.mtx",
    "matrix_int64.mtx",
    "matrix_uint8.mtx",
    "matrix_uint16.mtx",
    "matrix_uint32.mtx",
    "matrix_uint64.mtx",
    "matrix_fp32.mtx",
    "matrix_fp32_structure.mtx",
    "matrix_fp64.mtx",
    "west0067_jumbled.mtx",
    "skew_fp32.mtx",
    "skew_fp64.mtx",
    "skew_int8.mtx",
    "skew_int16.mtx",
    "skew_int32.mtx",
    "skew_int64.mtx",
    "structure.mtx",
    "full.mtx",
    "full_symmetric.mtx",
    "empty.mtx",
    "",
} ;

//****************************************************************************
void test_SWrite (void)
{
    LAGraph_Init (msg) ;
    GrB_Descriptor desc = NULL ;
    #if LG_SUITESPARSE
    OK (GrB_Descriptor_new (&desc)) ;
    OK (GxB_set (desc, GxB_COMPRESSION, GxB_COMPRESSION_LZ4HC + 9)) ;
    #endif

    for (int k = 0 ; k < NFILES ; k++)
    {

        // load the matrix as A
        const char *aname = files [k] ;
        if (strlen (aname) == 0) break;
        printf ("\n================================== %d %s:\n", k, aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        fclose (f) ;
        // GxB_print (A, 3) ;

        #if LG_SUITESPARSE
        for (int scon = 1 ; scon <= 8 ; scon = 2*scon)
        #endif
        {
            #if LG_SUITESPARSE
            // for SuiteSparse only: test all sparsity formats
            OK (GxB_set (A, GxB_SPARSITY_CONTROL, scon)) ;
            #endif

            // workaround for bug in v6.0.0 to v6.0.2:
            // ensure the matrix is not iso
            #if LG_SUITESPARSE
            #if GxB_IMPLEMENTATION < GxB_VERSION (6,0,3)
            printf ("workaround for bug in SS:GrB v6.0.2 (fixed in v6.0.3)\n") ;
            OK (GrB_Matrix_setElement (A, 0, 0, 0)) ;
            OK (GrB_wait (A, GrB_MATERIALIZE)) ;
            #endif
            #endif

            // open a temporary *.lagraph file to hold the matrix
            f = tmpfile ( ) ;
            // snprintf (filename, LEN, "%s.lagraph", aname) ;
            // f = fopen (filename, "w") ;
            TEST_CHECK (f != NULL) ;

            // serialize the matrix
            // GxB_set (GxB_BURBLE, true) ;
            bool ok ;
            void *blob = NULL ;
            GrB_Index blob_size = 0 ;
            #if LG_SUITESPARSE
            if (k % 2 == 0)
            {
                // for SuiteSparse: try GxB for every other matrix
                OK (GxB_Matrix_serialize (&blob, &blob_size, A, desc)) ;
            }
            else
            #endif
            {
                // try GrB version
                OK (GrB_Matrix_serializeSize (&blob_size, A)) ;
                GrB_Index blob_size_old = blob_size ;
                blob = LAGraph_Malloc (blob_size, sizeof (uint8_t)) ;
                TEST_CHECK (blob != NULL) ;
                OK (GrB_Matrix_serialize (blob, &blob_size, A)) ;
                blob = LAGraph_Realloc (blob_size, blob_size_old,
                    sizeof (uint8_t), blob, &ok) ;
                TEST_CHECK (ok) ;
            }

            // deserialize the matrix
            int rr = (GrB_Matrix_deserialize (&B, atype, blob, blob_size)) ;
            // printf ("A:\n") ; GxB_print (A, 2) ;
            // printf ("B:\n") ; GxB_print (B, 2) ;
            // printf ("rr: %d\n", rr) ;
            OK (rr) ;
            // GxB_set (GxB_BURBLE, false) ;

            // ensure the matrices A and B are the same
            // GxB_print (A,3) ;
            // GxB_print (B,3) ;
            OK (LAGraph_IsEqual_type (&ok, A, B, atype, msg)) ;
            TEST_CHECK (ok) ;
            OK (GrB_free (&B)) ;

            // get the name of the C typedef for the matrix
            char *typename ;
            OK (LAGraph_TypeName (&typename, atype, msg)) ;

            // write the header for a single matrix
            OK (LAGraph_SWrite_HeaderStart (f, "lagraph_test", msg)) ;
            OK (LAGraph_SWrite_HeaderItem (f, LAGraph_matrix_kind, "A",
                typename, 0, blob_size, msg)) ;
            OK (LAGraph_SWrite_HeaderEnd (f, msg)) ;

            // write the binary blob to the file then free the blob
            OK (LAGraph_SWrite_Item (f, blob, blob_size, msg)) ;
            LAGraph_Free (&blob) ;

            // open the file and load back the contents
            rewind (f) ;
            // fclose (f) ;
            // f = fopen (filename, "r") ;

            char *collection = NULL ;
            LAGraph_Contents *Contents = NULL ;
            GrB_Index ncontents ;
            OK (LAGraph_SRead (f, &collection, &Contents, &ncontents, msg)) ;
            TEST_CHECK (collection != NULL) ;
            if (collection == NULL) abort ( ) ;
            // printf ("collection %s\n", collection) ;
            TEST_CHECK (strcmp (collection, "lagraph_test") == 0) ;
            TEST_CHECK (ncontents == 1) ;
            fclose (f) ;

            // convert the contents to a matrix B
            void *blob2 = Contents [0].blob ;
            size_t blob_size2 = Contents [0].blob_size ;
            // printf ("blob_size2 %lu\n", blob_size2) ;
            TEST_CHECK (blob_size == blob_size2) ;

            OK (GrB_Matrix_deserialize (&B, atype, blob2, blob_size2)) ;
            // GxB_set (GxB_BURBLE, false) ;

            // ensure the matrices A and B are the same
            // GxB_print (A,3) ;
            // GxB_print (B,3) ;
            OK (LAGraph_IsEqual_type (&ok, A, B, atype, msg)) ;
            TEST_CHECK (ok) ;
            OK (GrB_free (&B)) ;

            // free the contents: TODO make this a utility function
            LAGraph_Free ((void **) &collection) ;
            for (int i = 0 ; i < ncontents ; i++)
            {
                LAGraph_Contents *Item = &(Contents [i]) ;
                LAGraph_Free ((void **) &(Item->blob)) ;
            }
            LAGraph_Free ((void **) &Contents) ;
        }

        OK (GrB_free (&A)) ;
    }

    OK (GrB_free (&desc)) ;
    LAGraph_Finalize (msg) ;
}

//****************************************************************************

TEST_LIST = {
    {"SWrite", test_SWrite},
    {NULL, NULL}
};
