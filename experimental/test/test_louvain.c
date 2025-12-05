#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>
#include "LG_Xtest.h"

#define dbg(x) GxB_print(x, 5)
#define err(x, info)                                    \
    if (!(info == GrB_SUCCESS || info == GrB_NO_VALUE)) \
    {                                                   \
        char **err;                                     \
        GrB_error(err, x);                              \
        printf("\ninfo: %d error: %s\n", info, err);    \
    }
char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;
GrB_Matrix A = NULL;
#define LEN 512
char filename[LEN + 1];
typedef struct
{
    const char *matrix_file; // Adjeancy matrix or graph
    const double mod;
} matrix_info;

const matrix_info files[] = {

    {"comm0.mtx", 0.357142857142857},
    // {"karate.mtx", .42},
    // {"50node.mtx", .0000},
    // {"20000node.mtx", .42},
    // {"50000node.mtx", .42},
    // {"ca-GrQc.mtx", .42},
    // {"email-Enron.mtx", .42},
    // {"com-Amazon.mtx", .42},
    // {"com-Youtube.mtx", .42},
    {"", -1}};


void test_LouvainSeq(void)
{
     LAGraph_Init(msg);

    for (int k = 0;; k++)
    {
        uint64_t seed = 1224;

        const char *aname = files[k].matrix_file;
        if (strlen(aname) == 0)
            break;
        printf("\n================================== %s:\n", aname);
        snprintf(filename, LEN, LG_DATA_DIR "%s", files[k].matrix_file);
        FILE *f = fopen(filename, "r");
        TEST_CHECK(f != NULL);
        OK(LAGraph_MMRead(&A, f, msg));
        fclose(f);

        OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
        TEST_CHECK(A == NULL);


        OK(LAGraph_Cached_AT(G, msg));
        // check if the pattern is symmetric - if it isn't make it.
        OK(LAGraph_Cached_IsSymmetricStructure(G, msg));
        GrB_Matrix S = NULL;
        double tsimple = LAGraph_WallClockTime();
        OK(LAGraph_LouvainSeq(&S, G,seed, msg));

        // OK(LAGraph_Louvain_res(&S,G,.3,msg));
        tsimple = LAGraph_WallClockTime() - tsimple;
        double Q = 0.0;
        OK(LAGr_Modularity2(&Q, 1.0, G->A, S, msg));
        printf("Q:%f\n", Q);
        // printf("Number of Communities: %d",comms);
        printf(" time: %f\n", tsimple);
        OK(LAGraph_Delete(&G, msg));
    }
    LAGraph_Finalize(msg);
}
void test_LouvainIS(void)
{
     LAGraph_Init(msg);

    for (int k = 0;; k++)
    {
        uint64_t seed = 1249141465;
        const char *aname = files[k].matrix_file;
        if (strlen(aname) == 0)
            break;
        printf("\n================================== %s:\n", aname);
        snprintf(filename, LEN, LG_DATA_DIR "%s", files[k].matrix_file);
        FILE *f = fopen(filename, "r");
        TEST_CHECK(f != NULL);
        OK(LAGraph_MMRead(&A, f, msg));
        fclose(f);

        OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
        TEST_CHECK(A == NULL);


        OK(LAGraph_Cached_AT(G, msg));
        // check if the pattern is symmetric - if it isn't make it.
        OK(LAGraph_Cached_IsSymmetricStructure(G, msg));
        GrB_Matrix S = NULL;
        double tsimple = LAGraph_WallClockTime();
        OK(LAGraph_LouvainIS(&S,seed, G, msg));

        // OK(LAGraph_Louvain_res(&S,G,.3,msg));
        tsimple = LAGraph_WallClockTime() - tsimple;
        double Q = 0.0; 
        double tsimple2 = LAGraph_WallClockTime();
        OK(LAGr_Modularity2(&Q, 1.0, G->A, S, msg));
        tsimple2 = LAGraph_WallClockTime() - tsimple2;

        printf("Q:%f time to calc Q: %f\n", Q,tsimple2);
        // printf("Number of Communities: %d",comms);
        printf(" time: %f\n", tsimple);
        GrB_free(&S);
        OK(LAGraph_Delete(&G, msg));

    }

    LAGraph_Finalize(msg);
}

TEST_LIST = {

    {"LouvainSeq", test_LouvainSeq},
    {"LouvainIS", test_LouvainIS},

    {NULL, NULL}};



    