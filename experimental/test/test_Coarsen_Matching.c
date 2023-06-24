#include <stdio.h>
#include <acutest.h>
#include "LAGraphX.h"
#include "LAGraph_test.h"
#include "LG_Xtest.h"
#include "LG_internal.h"

char msg [LAGRAPH_MSG_LEN] ;

/*
- list of graphs (generate randomly, similar method as matching test)
- do coarsening independently (use another formulation)
    - How to do matching in independent method?
    - Can try to call matching w/ seed in independent method
- independent method should write coarsened graph to file
    - alternate methods? Possible to directly call?
    - main thing: need to be able to create 
- do coarsening using Coarsen_Matching in test
- MMRead from file from independent method, check equality of graphs
*/

void test_Coarsen_Matching () {

}

TEST_LIST = {
    {"Coarsen_Matching", test_Coarsen_Matching},
    {NULL, NULL}
} ;