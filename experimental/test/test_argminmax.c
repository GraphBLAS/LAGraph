// FIXME: this test needs to check its results.
// Use a brute force method (extractTuples and do it in plain C, perhaps).

#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <LAGraph.h>
#include <LG_internal.h>

char msg [LAGRAPH_MSG_LEN] ;

#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] =
{
    { "structure.mtx" },
    { "karate.mtx" },
    { "west0067.mtx" },
    { "bcsstk13.mtx" },
    { "" },
} ;

void test_argminmax (void)
{

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------

    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL, C = NULL ;
    GrB_Matrix x = NULL, p = NULL ;
    GrB_Index nrows, ncols ;

    for (int k = 0 ; ; k++)
    {
        // load the matrix as A
        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break ;
        printf ("\n %s: ==================================\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;
        fclose (f) ;
        OK (GrB_Matrix_nrows (&nrows, A)) ;
        OK (GrB_Matrix_ncols (&ncols, A)) ;

        for (int ktype = 0 ; ktype < 11 ; ktype++)
        {
            GrB_Type type ;
            switch (ktype)
            {
                case  0: type = GrB_BOOL    ; break ;
                case  1: type = GrB_INT8    ; break ;
                case  2: type = GrB_INT16   ; break ;
                case  3: type = GrB_INT32   ; break ;
                case  4: type = GrB_INT64   ; break ;
                case  5: type = GrB_UINT8   ; break ;
                case  6: type = GrB_UINT16  ; break ;
                case  7: type = GrB_UINT32  ; break ;
                case  8: type = GrB_UINT64  ; break ;
                case  9: type = GrB_FP32    ; break ;
                default:
                case 10: type = GrB_FP64    ; break ;
            }

            // typecast A into a different type
            OK (GrB_Matrix_new (&C, type, nrows, ncols)) ;
            OK (GrB_assign (C, NULL, NULL, A,
                GrB_ALL, nrows, GrB_ALL, ncols, NULL)) ;

            printf ("\nA:\n") ;
            OK (LAGraph_Matrix_Print (A, 2, stdout, msg)) ;

            printf ("\nC:\n") ;
            OK (LAGraph_Matrix_Print (C, 2, stdout, msg)) ;

            for (int is_min = 0 ; is_min <= 1 ; is_min++)
            {
                for (int dim = 0 ; dim <= 2 ; dim++)
                {
                    printf ("\nis_min: %d dim: %d\n", is_min, dim) ;
                    // test the algorithm
                    OK (LAGraph_argminmax (&x, &p, C, dim, is_min, msg)) ;
                    // print the result
                    printf ("\nx:\n") ;
                    OK (LAGraph_Matrix_Print (x, 2, stdout, msg)) ;
                    printf ("\np:\n") ;
                    OK (LAGraph_Matrix_Print (p, 2, stdout, msg)) ;
                    OK (GrB_free (&x)) ;
                    OK (GrB_free (&p)) ;
                }
            }
            OK (GrB_free (&C)) ;
        }
        OK (GrB_free (&A)) ;
    }

    //--------------------------------------------------------------------------
    // finalize LAGraph
    //--------------------------------------------------------------------------

    LAGraph_Finalize (msg) ;
}

//----------------------------------------------------------------------------
// the main program is created by acutest, and it runs a list of tests:
//----------------------------------------------------------------------------

TEST_LIST =
{
    {"argminmax", test_argminmax},
    {NULL, NULL}
} ;

