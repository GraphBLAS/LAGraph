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

/**
 * Code is based on the algorithm described in the following paper
 * Zhang, Azad, Hu. FastSV: FastSV: A Distributed-Memory Connected Component Algorithm with Fast Convergence (SIAM PP20)
 **/

#define LAGRAPH_FREE_ALL

#include "LAGraph.h"

// mask = NULL, accumulator = GrB_MIN_UINT64, descriptor = NULL
static GrB_Info Reduce_assign (GrB_Vector w, GrB_Vector src, GrB_Index *index, GrB_Index nLocs)
{
    GrB_Index nw, ns;
    GrB_Vector_nvals(&nw, w);
    GrB_Vector_nvals(&ns, src);
    GrB_Index *mem = (GrB_Index*) malloc(sizeof(GrB_Index) * nw * 3);
    GrB_Index *ind = mem, *sval = mem + nw, *wval = sval + nw;
    GrB_Vector_extractTuples(ind, wval, &nw, w);
    GrB_Vector_extractTuples(ind, sval, &ns, src);
    for (GrB_Index i = 0; i < nLocs; i++)
        if (sval[i] < wval[index[i]])
            wval[index[i]] = sval[i];
    GrB_Vector_clear(w);
    GrB_Vector_build(w, ind, wval, nw, GrB_PLUS_UINT64);
    free(mem);
    return GrB_SUCCESS;
}

GrB_Info LAGraph_cc_fastsv
(
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A,           // input matrix
    bool sanitize           // if true, ensure A is symmetric
)
{
    GrB_Info info;
    GrB_Index n ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;
    //GrB_Index nnz ;
    //LAGRAPH_OK (GrB_Matrix_nvals (&nnz, A)) ;
    //printf ("number of nodes: %g\n", (double) n) ;
    //printf ("number of edges: %g\n", (double) nnz) ;

    GrB_Matrix S = NULL;
    if (sanitize)
    {
        GrB_Descriptor desc = NULL ;
        LAGRAPH_OK(GrB_Descriptor_new(&desc)) ;
        LAGRAPH_OK(GrB_Descriptor_set(desc, GrB_INP1, GrB_TRAN)) ;

        LAGRAPH_OK (GrB_Matrix_new (&S, GrB_BOOL, n, n)) ;
        LAGRAPH_OK (GrB_eWiseAdd (S, NULL, NULL, GrB_LOR, A, A, desc)) ;
        LAGRAPH_FREE(desc) ;
    }
    else
    {
        // Use the input as-is, and assume it is binary and symmetric
        S = A ;
    }

    // vectors
    GrB_Vector f, dup, mngp;
    GrB_Vector mod, gp;
    LAGRAPH_OK (GrB_Vector_new(&f,   GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new(&mod, GrB_UINT64, n));
    // temporary arrays
    GrB_Index *I = malloc(sizeof(GrB_Index) * n);
    GrB_Index *V = malloc(sizeof(GrB_Index) * n);
    // prepare vectors
    for (GrB_Index i = 0; i < n; i++)
        I[i] = V[i] = i;
    LAGRAPH_OK (GrB_Vector_build (f, I, V, n, GrB_PLUS_UINT64));
    LAGRAPH_OK (GrB_Vector_dup (&gp,  f));
    LAGRAPH_OK (GrB_Vector_dup (&dup, f));
    LAGRAPH_OK (GrB_Vector_dup (&mngp,f));
    // semiring & monoid
    GrB_Monoid Min, Add;
    GrB_Semiring sel2ndMin;
    LAGRAPH_OK (GrB_Monoid_new (&Min, GrB_MIN_UINT64, n));
    LAGRAPH_OK (GrB_Semiring_new (&sel2ndMin, Min, GrB_SECOND_UINT64));
    LAGRAPH_OK (GrB_Monoid_new (&Add, GrB_PLUS_UINT64, (GrB_Index) 0));
    // main computation
    GrB_Index diff = n;
    while (diff != 0) {
        // hooking & shortcutting
        LAGRAPH_OK (GrB_mxv (mngp, 0, GrB_MIN_UINT64, sel2ndMin, S, gp, 0));
        LAGRAPH_OK (Reduce_assign (f, mngp, V, n));
        LAGRAPH_OK (GrB_eWiseMult (f, 0, 0, GrB_MIN_UINT64, f, mngp, 0));
        LAGRAPH_OK (GrB_eWiseMult (f, 0, 0, GrB_MIN_UINT64, f, gp, 0));
        // calculate grandparent
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, f));
        LAGRAPH_OK (GrB_extract (gp, 0, 0, f, V, n, 0));
        // check termination
        LAGRAPH_OK (GrB_eWiseMult (mod, 0, 0, GxB_ISNE_UINT64, dup, gp, 0));
        LAGRAPH_OK (GrB_reduce (&diff, 0, Add, mod, 0));
        LAGRAPH_OK (GrB_assign (dup, 0, 0, gp, GrB_ALL, 0, 0));
    }
    *result = f;
    // free
    free(I);
    free(V);

    LAGRAPH_OK (GrB_free (&gp));
    LAGRAPH_OK (GrB_free (&mngp));
    LAGRAPH_OK (GrB_free (&dup));
    LAGRAPH_OK (GrB_free (&mod));
    LAGRAPH_OK (GrB_free (&Add));
    LAGRAPH_OK (GrB_free (&Min));
    LAGRAPH_OK (GrB_free (&sel2ndMin));
    return GrB_SUCCESS;
}

