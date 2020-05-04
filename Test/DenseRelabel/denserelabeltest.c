//------------------------------------------------------------------------------
// LAGraph/Test/DenseRelabel/denserelabeltest.c: test program for LAGraph_dense_relabel
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

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// Contributed by Marton Elekes and Gabor Szarnyas, BME

// Usage:
//
// denserelabeltest

#include "LAGraph.h"

#include <inttypes.h>

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&MMapping) ;          \
    GrB_free (&VMapping) ;          \
}

int main(int argc, char **argv) {

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Info info;
    GrB_Matrix MMapping = NULL;
    GrB_Vector VMapping = NULL;

    LAGRAPH_OK (LAGraph_init());

    GrB_Index big_id = ((GrB_Index) 1) << 48;
    GrB_Index identifiers[] = {42, 0, big_id, 1};
    GrB_Index nids = sizeof(identifiers) / sizeof(identifiers[0]);

    for (size_t i = 0; i < nids; ++i) {
        printf("%" PRIu64 "\n", identifiers[i]);
    }

    LAGRAPH_OK(LAGraph_dense_relabel(&MMapping, &VMapping, identifiers, nids));

    LAGRAPH_OK(GxB_fprint(MMapping, GxB_COMPLETE, stdout));
    LAGRAPH_OK(GxB_fprint(VMapping, GxB_COMPLETE, stdout));


    LAGRAPH_FREE_ALL;
    LAGraph_finalize();
}
