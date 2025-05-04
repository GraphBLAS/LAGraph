#include <stdio.h>
#include <acutest.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>

#define LEN 512
#define MAX_LABELS 3
#define MAX_RESULTS 2000000

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G[MAX_LABELS] ;
LAGraph_Graph R[MAX_LABELS] ;
GrB_Matrix A ;

char testcase_name [LEN+1] ;
char filename [LEN+1] ;

typedef struct
{
    const char* name ;
    const char* graphs[MAX_LABELS] ;
    const char* fas[MAX_LABELS] ;
    const char* fa_meta ;
    const char* sources ;
    const GrB_Index expected[MAX_RESULTS] ;
    const size_t expected_count ;
}
matrix_info ;

const matrix_info files [ ] =
{
    {"simple 1 or more",
     {"rpq_data/a.mtx",   "rpq_data/b.mtx",   NULL},
     {"rpq_data/1_a.mtx", NULL },                    // Regex: a+
     "rpq_data/1_meta.txt",
     "rpq_data/1_sources.txt",
     {2, 4, 6, 7}, 4},
    {"simple kleene star",
     {"rpq_data/a.mtx",   "rpq_data/b.mtx",   NULL},
     {"rpq_data/2_a.mtx", "rpq_data/2_b.mtx", NULL}, // Regex: (a b)*
     "rpq_data/2_meta.txt",
     "rpq_data/2_sources.txt",
     {2, 6, 8}, 3},
    {"kleene star of the conjunction",
     {"rpq_data/a.mtx",   "rpq_data/b.mtx",   NULL},
     {"rpq_data/3_a.mtx", "rpq_data/3_b.mtx", NULL}, // Regex: (a | b)*
     "rpq_data/3_meta.txt",
     "rpq_data/3_sources.txt",
     {3, 6}, 2},
    {"simple repeat from n to m times",
     {"rpq_data/a.mtx",   "rpq_data/b.mtx",   NULL},
     {"",                 "rpq_data/4_b.mtx", NULL}, // Regex: b b b (b b)?
     "rpq_data/4_meta.txt",
     "rpq_data/4_sources.txt",
     {3, 4, 6}, 3},
    {NULL, NULL, NULL, NULL},
} ;

//****************************************************************************
void test_RegularPathQueryBasic (void)
{
    LAGraph_Init (msg) ;

    for (int k = 0 ; ; k++)
    {
        if (files[k].sources == NULL) break ;

        snprintf (testcase_name, LEN, "basic regular path query %s", files[k].name) ;
        TEST_CASE (testcase_name) ;

        for (int check_symmetry = 0 ; check_symmetry <= 1 ; check_symmetry++)
        {
            // Load graph from MTX files representing its adjacency matrix
            // decomposition
            for (int i = 0 ; ; i++)
            {
                const char *name = files[k].graphs[i] ;

                if (name == NULL) break ;
                if (strlen(name) == 0) continue ;

                snprintf (filename, LEN, LG_DATA_DIR "%s", name) ;
                FILE *f = fopen (filename, "r") ;
                TEST_CHECK (f != NULL) ;
                OK (LAGraph_MMRead (&A, f, msg)) ;
                OK (fclose (f));

                OK (LAGraph_New (&(G[i]), &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;

                TEST_CHECK (A == NULL) ;
            }

            // Load NFA from MTX files representing its adjacency matrix
            // decomposition
            for (int i = 0 ; ; i++)
            {
                const char *name = files[k].fas[i] ;

                if (name == NULL) break ;
                if (strlen(name) == 0) continue ;

                snprintf (filename, LEN, LG_DATA_DIR "%s", name) ;
                FILE *f = fopen (filename, "r") ;
                TEST_CHECK (f != NULL) ;
                OK (LAGraph_MMRead (&A, f, msg)) ;
                OK (fclose (f)) ;

                OK (LAGraph_New (&(R[i]), &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;

                if (check_symmetry)
                {
                    // Check if the pattern is symmetric - if it isn't make it.
                    // Note this also computes R[i]->AT
                    OK (LAGraph_Cached_IsSymmetricStructure (R[i], msg)) ;
                }

                TEST_CHECK (A == NULL) ;
            }

            // Note the matrix rows/cols are enumerated from 0 to n-1.
            // Meanwhile, in MTX format they are enumerated from 1 to n. Thus,
            // when loading/comparing the results these values should be
            // decremented/incremented correspondingly.

            // Load graph source nodes from the sources file
            GrB_Index s ;
            GrB_Index S[16] ;
            size_t ns = 0 ;

            const char *name = files[k].sources ;
            snprintf (filename, LEN, LG_DATA_DIR "%s", name) ;
            FILE *f = fopen (filename, "r") ;
            TEST_CHECK (f != NULL) ;

            while (fscanf(f, "%" PRIu64, &s) != EOF)
            {
                S[ns++] = s - 1 ;
            }

            OK (fclose(f)) ;

            // Load NFA starting states from the meta file
            GrB_Index qs ;
            GrB_Index QS[16] ;
            size_t nqs = 0 ;

            name = files[k].fa_meta ;
            snprintf (filename, LEN, LG_DATA_DIR "%s", name) ;
            f = fopen (filename, "r") ;
            TEST_CHECK (f != NULL) ;

            uint64_t nqs64 = 0 ;
            TEST_CHECK (fscanf(f, "%" PRIu64, &nqs64) != EOF) ;
            nqs = (size_t) nqs64 ;

            for (uint64_t i = 0; i < nqs; i++) {
                TEST_CHECK (fscanf(f, "%" PRIu64, &qs) != EOF) ;
                QS[i] = qs - 1 ;
            }

            // Load NFA final states from the same file
            uint64_t qf ;
            uint64_t QF[16] ;
            size_t nqf = 0 ;
            uint64_t  nqf64 = 0 ;

            TEST_CHECK (fscanf(f, "%" PRIu64, &nqf64) != EOF) ;
            nqf = (size_t) nqf64 ;

            for (uint64_t i = 0; i < nqf; i++) {
                TEST_CHECK (fscanf(f, "%" PRIu64, &qf) != EOF) ;
                QF[i] = qf - 1 ;
            }

            OK (fclose(f)) ;

            // Evaluate the algorithm
            GrB_Vector r = NULL ;

            OK (LAGraph_RegularPathQuery (&r, R, MAX_LABELS, QS, nqs,
                                             QF, nqf, G, S, ns, msg)) ;

            // Extract results from the output vector
            GrB_Index *reachable ;
            bool *values ;

            GrB_Index nvals ;
            GrB_Vector_nvals (&nvals, r) ;

            OK (LAGraph_Malloc ((void **) &reachable, MAX_RESULTS, sizeof (GrB_Index), msg)) ;
            OK (LAGraph_Malloc ((void **) &values, MAX_RESULTS, sizeof (GrB_Index), msg)) ;

            GrB_Vector_extractTuples (reachable, values, &nvals, r) ;

            // Compare the results with expected values
            TEST_CHECK (nvals == files[k].expected_count) ;
            for (uint64_t i = 0 ; i < nvals ; i++)
                TEST_CHECK (reachable[i] + 1 == files[k].expected[i]) ;

            // Cleanup
            OK (LAGraph_Free ((void **) &values, NULL)) ;
            OK (LAGraph_Free ((void **) &reachable, NULL)) ;

            OK (GrB_free (&r)) ;

            for (uint64_t i = 0 ; i < MAX_LABELS ; i++)
            {
                if (G[i] == NULL) continue ;
                OK (LAGraph_Delete (&(G[i]), msg)) ;
            }

            for (uint64_t i = 0 ; i < MAX_LABELS ; i++ )
            {
                if (R[i] == NULL) continue ;
                OK (LAGraph_Delete (&(R[i]), msg)) ;
            }
        }
    }

    LAGraph_Finalize (msg) ;
}


TEST_LIST = {
    {"RegularPathQueryBasic", test_RegularPathQueryBasic},
    {NULL, NULL}
};
