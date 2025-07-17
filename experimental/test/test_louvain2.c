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
    const double mod;
} matrix_info;

const matrix_info files[] = {

    {"comm0.mtx", 0.357142857142857},
    {"res1.mtx", 0.0},
    {"karate2.mtx", .42},
    {"",-1}  
};
//Store matrix by row
void test_Louvain(void){
    LAGraph_Init(msg);
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
        OK (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;

        if (G->is_symmetric_structure == LAGraph_FALSE)
        {
            printf("This matrix is not symmetric. \n");
            // make the adjacency matrix symmetric
            OK (LAGraph_Cached_AT (G, msg)) ;
            OK (GrB_eWiseAdd (G->A, NULL, NULL, GrB_LOR, G->A, G->AT, NULL)) ;
            G->is_symmetric_structure = true ;
            // consider the graph as directed
            G->kind = LAGraph_ADJACENCY_DIRECTED ;
        }
        else
        {
            G->kind = LAGraph_ADJACENCY_UNDIRECTED ;
        }
        GrB_Matrix S;
        double tsimple = LAGraph_WallClockTime ( ) ;
        OK(LAGraph_Louvain2(&S,G,msg));
        tsimple = LAGraph_WallClockTime ( ) - tsimple ;
        printf(" time: %f\n",tsimple);


    }
}
TEST_LIST = {
    {"Louvain", test_Louvain},
    {NULL, NULL}
};