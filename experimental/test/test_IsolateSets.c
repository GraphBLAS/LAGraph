#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>
#include "LG_Xtest.h"

char msg [LAGRAPH_MSG_LEN];
LAGraph_Graph G;
GrB_Matrix A;
#define LEN 512
char filename [LEN+1];
typedef struct
{
    const char * matrix_file; //Adjeancy matrix or graph
} matrix_info;
const matrix_info files[] = {

    {"comm1.mtx"},
    {"comm0.mtx"},
    {"res1.mtx"},
    {"karate.mtx"},
    {""} 
};

void test_IsolateSets(void){
    LAGraph_Init(msg);
    OK (LAGraph_Random_Init (msg)) ;
    printf("\n");
    for(int k = 0;;k++){
        if (strlen(files[k].matrix_file) == 0)
            break;
        snprintf (filename, LEN, LG_DATA_DIR "%s", files[k].matrix_file) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;

        // check if the pattern is symmetric - if it isn't make it.
        OK (LAGraph_Cached_OutDegree (G, msg)) ;
        OK (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;
        
        GrB_Vector Iset;
        double tsimple = LAGraph_WallClockTime ( ) ;
        uint64_t seed = 122 ^ (uint64_t)tsimple << 2044542;
        OK(LAGraph_IsolateSets(&Iset,G,seed,msg));
        GxB_print(Iset,5);
        tsimple = LAGraph_WallClockTime ( ) - tsimple ;
        printf(" time: %f\n",tsimple);


    }
}
TEST_LIST = {
    {"Isolate Set", test_IsolateSets},
    {NULL, NULL}
};
