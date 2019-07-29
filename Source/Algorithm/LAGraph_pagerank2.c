//------------------------------------------------------------------------------
// LAGraph_pagerank2: pagerank using a real semiring
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

// LAGraph_pagerank2: Alternative PageRank implementation using a real
// semiring. Not to be confused with the dpagerank2.c algorithm, included
// in the SuiteSparse:GrapBLAS demos.
//
// This algorithm follows the specification given in the LDBC Graphalytics
// benchmark, see https://github.com/ldbc/ldbc_graphalytics_docs/
//
// Contributed by Gabor Szarnyas and Balint Hegyi, Budapest University of
// Technology and Economics (with accented characters: G\'{a}bor Sz\'{a}rnyas
// and B\'{a}lint Hegyi, using LaTeX syntax).

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL { \
    GrB_Descriptor_free(&transpose_desc); \
    GrB_Descriptor_free(&invmask_desc); \
    GrB_Matrix_free(&A); \
    GrB_Vector_free(&A_out); \
    GrB_Vector_free(&nondangling_mask); \
    GrB_Vector_free(&importance_vec); \
    GrB_Vector_free(&dangling_vec); \
    GrB_Vector_free(pr_next); \
};

GrB_Info LAGraph_pagerank2(
    GrB_Vector *result,
    GrB_Matrix A,
    double damping_factor,
    unsigned long iteration_num
) {
    GrB_Info info;
    GrB_Index n;

    GrB_Descriptor invmask_desc;
    GrB_Descriptor transpose_desc;
    GrB_Vector A_out;
    GrB_Vector nondangling_mask;

    GrB_Vector importance_vec;
    GrB_Vector dangling_vec;

    GrB_Vector *pr_prev = NULL;
    GrB_Vector *pr_next = NULL;

    LAGRAPH_OK(GrB_Matrix_nrows(&n, A))
    GrB_Index nvals;
    LAGRAPH_OK(GrB_Matrix_nvals(&nvals, A))

    // Make a complement descriptor

    LAGRAPH_OK(GrB_Descriptor_new(&invmask_desc))
    LAGRAPH_OK(GrB_Descriptor_set(invmask_desc, GrB_MASK, GrB_SCMP))

    // Make a transpose descriptor

    LAGRAPH_OK(GrB_Descriptor_new(&transpose_desc))
    LAGRAPH_OK(GrB_Descriptor_set(transpose_desc, GrB_INP0, GrB_TRAN))
    LAGRAPH_OK(GrB_Descriptor_set(transpose_desc, GrB_OUTP, GrB_REPLACE))

    //
    // Matrix A row sum
    //
    // Stores the outbound degrees of all vertices
    //
    LAGRAPH_OK(GrB_Vector_new(&A_out, GrB_UINT64, n))
    LAGRAPH_OK(GrB_Matrix_reduce_Monoid(
        A_out,
        GrB_NULL,
        GrB_NULL,
        GxB_PLUS_UINT64_MONOID,
        A,
        GrB_NULL
    ))

    //
    // Non-dangling_vec vector determination
    //
    // These vertices the ones which have outgoing edges. In further operations,
    // this mask can be negated to select dangling_vec vertices.
    //
    LAGRAPH_OK(GrB_Vector_new(&nondangling_mask, GrB_BOOL, n))
    LAGRAPH_OK(GrB_Matrix_reduce_Monoid(
        nondangling_mask,
        GrB_NULL,
        GrB_NULL,
        GxB_LOR_BOOL_MONOID,
        A,
        GrB_NULL
    ))

    //
    // Iteration
    //

    // Result vectors

    LAGRAPH_OK(GrB_Vector_new(pr_prev, GrB_FP64, n))
    LAGRAPH_OK(GrB_Vector_new(pr_next, GrB_FP64, n))

    // Fill result vector with initial value (1 / |V|)
    LAGRAPH_OK(GrB_Vector_assign_FP64(
        *pr_prev,
        GrB_NULL,
        GrB_NULL,
        (1.0 / n),
        GrB_ALL,
        n,
        GrB_NULL
    ))

    LAGRAPH_OK(GrB_Vector_new(&importance_vec, GrB_FP64, n))
    LAGRAPH_OK(GrB_Vector_new(&dangling_vec, GrB_FP64, n))

    // Teleport value
    const double teleport = (1 - damping_factor) / n;

    pr_prev = (GrB_Vector *) malloc(sizeof(GrB_Vector));
    if (pr_prev == NULL) {
        LAGRAPH_ERROR("pr_prev allocation error, pointer is NULL", GrB_NULL_POINTER)
    }
    pr_prev = (GrB_Vector *) malloc(sizeof(GrB_Vector));
    if (pr_prev == NULL) {
        LAGRAPH_ERROR("pr_next allocation error, pointer is NULL", GrB_NULL_POINTER)
    }

    pr_prev = (GrB_Vector *) malloc(sizeof(GrB_Vector));
    pr_next = (GrB_Vector *) malloc(sizeof(GrB_Vector));

    for (int i = 0; i < iteration_num; i++) {
        //
        // Importance calculation
        //

        // Divide previous PageRank with number of outbound edges
        LAGRAPH_OK(GrB_eWiseMult_Vector_BinaryOp(
            importance_vec,
            nondangling_mask,
            GrB_NULL,
            GrB_DIV_FP64,
            *pr_prev,
            A_out,
            GrB_NULL
        ))

        // Multiply importance with damping factor
        LAGRAPH_OK(GrB_Vector_assign_FP64(
            importance_vec,
            nondangling_mask,
            GrB_TIMES_FP64,
            damping_factor,
            GrB_ALL,
            n,
            GrB_NULL
        ))

        // Calculate total PR of all inbound vertices
        LAGRAPH_OK(GrB_mxv(
            importance_vec,
            nondangling_mask,
            GrB_NULL,
            GxB_PLUS_TIMES_FP64,
            A,
            importance_vec,
            transpose_desc
        ))

        //
        // Dangling calculation
        //

        // Extract all the dangling PR entries from the previous result
        LAGRAPH_OK(GrB_Vector_extract(
            dangling_vec,
            nondangling_mask,
            GrB_NULL,
            *pr_prev,
            GrB_ALL,
            n,
            invmask_desc
        ))

        // Sum the previous PR values together
        double dangling_sum;
        LAGRAPH_OK(GrB_Vector_reduce_FP64(
            &dangling_sum,
            GrB_NULL,
            GxB_PLUS_FP64_MONOID,
            dangling_vec,
            GrB_NULL
        ))
        dangling_sum *= (damping_factor / n);

        //
        // PageRank summarization
        // Add teleport, importance_vec, and dangling_vec components together
        //
        LAGRAPH_OK(GrB_Vector_assign_FP64(
            *pr_next,
            GrB_NULL,
            GrB_NULL,
            (teleport + dangling_sum),
            GrB_ALL,
            n,
            GrB_NULL
        ))
        LAGRAPH_OK(GrB_eWiseAdd_Vector_Monoid(
            *pr_next,
            GrB_NULL,
            GrB_NULL,
            GxB_PLUS_FP64_MONOID,
            *pr_next,
            importance_vec,
            GrB_NULL
        ))

        // Swap the new and old vectors
        GrB_Vector *temp = pr_next;
        pr_next = pr_prev;
        pr_prev = temp;
    }

    // Because of the swap, the previous result will be the actual output
    result = pr_prev;
}
