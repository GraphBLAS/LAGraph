
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
 * Azad, Buluç. LACC: a linear-algebraic algorithm for finding connected components in distributed memory (IPDPS 2019)
 **/

#define LAGRAPH_FREE_ALL            \
{                                   \
}

#include "LAGraph.h"

GrB_Info CondHook(GrB_Matrix A, GrB_Vector parents, GrB_Vector stars)
{
    GrB_Info info;
    GrB_Index n ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;

    // Create (Sel2nd,Min)semiring
    GrB_Monoid Min = NULL ;                // Min monoid
    GrB_Semiring Sel2ndMin = NULL ;          // (Sel2nd,Min)semiring
    LAGRAPH_OK (GrB_Monoid_new (&Min, GrB_MIN_UINT64, (uint64_t)UINT_MAX) );
    LAGRAPH_OK (GrB_Semiring_new (&Sel2ndMin, Min, GrB_SECOND_UINT64) );

    // identify minNeighborparent for star vertices
    // hook stores minNeighborparent
    GrB_Vector hook = NULL;
    LAGRAPH_OK (GrB_Vector_new (&hook, GrB_UINT64, n));
    LAGRAPH_OK (GrB_mxv (hook, stars, NULL, Sel2ndMin,  A, parents, NULL)) ;
    //GxB_Vector_fprint(hook, "--Hook --", GxB_SHORT, stderr);

    
    GrB_Vector mask = NULL ;
    LAGRAPH_OK (GrB_Vector_new (&mask, GrB_BOOL, n));
    LAGRAPH_OK (GrB_eWiseMult(mask,NULL, NULL, GxB_ISLT_UINT64, hook,parents, NULL));
    //GxB_Vector_fprint(mask, "---- Mask ------", GxB_SHORT, stderr);
    
    GrB_Descriptor descReplace = NULL ;
    LAGRAPH_OK (GrB_Descriptor_new (&descReplace) );
    LAGRAPH_OK (GrB_Descriptor_set (descReplace, GrB_OUTP, GrB_REPLACE)) ;
    LAGRAPH_OK (GrB_assign(hook, mask, NULL, hook, GrB_ALL, n, descReplace));
    //GxB_Vector_fprint(hook, "---- Hook ------", GxB_SHORT, stderr);
    //GrB_eWiseMult(hook,mask, NULL, GrB_MIN_UINT64, hook,parents, descReplace);
    //parents of hooks
    GrB_Vector hookP = NULL;
    LAGRAPH_OK (GrB_Vector_new (&hookP, GrB_UINT64, n));
    LAGRAPH_OK (GrB_eWiseMult(hookP,NULL, NULL, GrB_SECOND_UINT64, hook,parents, NULL));
    //GxB_Vector_fprint(hookP, "---- hookP ------", GxB_SHORT, stderr);

    // extract values in hookP
    GrB_Index nhooks;
    LAGRAPH_OK (GrB_Vector_nvals(&nhooks, hook));
    GrB_Index *nzid = LAGraph_malloc (nhooks, sizeof (GrB_Index));
    GrB_Index *p = LAGraph_malloc (nhooks, sizeof (GrB_Index));
    if (nzid == NULL || p == NULL)
    {
        LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
    }
    LAGRAPH_OK (GrB_Vector_extractTuples(nzid, p, &nhooks, hookP));

    //a dense vector of hooks for assignments
    GrB_Vector hook_dense = NULL;
    LAGRAPH_OK (GrB_Vector_new (&hook_dense, GrB_UINT64, nhooks));
    LAGRAPH_OK (GrB_extract (hook_dense, NULL, NULL, hook,  nzid, nhooks, NULL)) ;
    //GxB_Vector_fprint(hook_dense, "---- hook_dense ------", GxB_SHORT, stderr);

    // update grandparents of hooks
    LAGRAPH_OK (GrB_assign (parents, NULL, NULL, hook_dense, p,  nhooks, NULL)) ;
    //GxB_Vector_fprint(parents, "---- parents****** ------", GxB_SHORT, stderr);
    
    LAGRAPH_FREE (nzid);
    LAGRAPH_FREE (p);
    GrB_free (&hook_dense) ;
    GrB_free (&hookP) ;
    GrB_free (&descReplace) ;
    GrB_free (&mask) ;
    GrB_free (&hook) ;
    GrB_free (&Min) ;
    GrB_free (&Sel2ndMin) ;
    
    return (GrB_SUCCESS) ;
}


GrB_Info UnCondHook(GrB_Matrix A, GrB_Vector parents, GrB_Vector stars)
{
    GrB_Info info;
    GrB_Index n ;
    GrB_Matrix_nrows (&n, A) ;

    // Create (Sel2nd,Min)semiring
    GrB_Monoid Min = NULL ;                // Min monoid
    GrB_Semiring Sel2ndMin = NULL ;          // (Sel2nd,Min)semiring
    LAGRAPH_OK (GrB_Monoid_new (&Min, GrB_MIN_UINT64, (uint64_t)UINT_MAX) );
    LAGRAPH_OK (GrB_Semiring_new (&Sel2ndMin, Min, GrB_SECOND_UINT64)) ;


    // extract parents of nonstar vertices
    GrB_Vector pNonstars = NULL ;
    LAGRAPH_OK (GrB_Vector_new (&pNonstars, GrB_UINT64, n));
    GrB_Descriptor descNonstars = NULL ;
    LAGRAPH_OK (GrB_Descriptor_new (&descNonstars)) ;
    LAGRAPH_OK (GrB_Descriptor_set (descNonstars, GrB_MASK, GrB_SCMP)) ;     // invert the mask
    LAGRAPH_OK (GrB_extract (pNonstars, stars, NULL, parents,  GrB_ALL, 0, descNonstars)) ;
    //GxB_Vector_fprint(pNonstars, "--pNonstars--", GxB_SHORT, stderr);


    // identify minNeighborparent for star vertices
    // hook stores minNeighborparent
    GrB_Vector hook = NULL;
    LAGRAPH_OK (GrB_Vector_new (&hook, GrB_UINT64, n));
    LAGRAPH_OK (GrB_mxv (hook, stars, NULL, Sel2ndMin,  A, pNonstars, NULL)) ;
    //GxB_Vector_fprint(hook, "--Hook --", GxB_SHORT, stderr);

    //parents of hooks
    GrB_Vector hookP = NULL;
    LAGRAPH_OK (GrB_Vector_new (&hookP, GrB_UINT64, n));
    LAGRAPH_OK (GrB_eWiseMult(hookP,NULL, NULL, GrB_SECOND_UINT64, hook,parents, NULL));
    //GxB_Vector_fprint(hookP, "---- hookP ------", GxB_SHORT, stderr);


    // extract values in hookP
    GrB_Index nhooks;
    LAGRAPH_OK (GrB_Vector_nvals(&nhooks, hook));
    GrB_Index *nzid = LAGraph_malloc (nhooks, sizeof (GrB_Index));
    GrB_Index *p = LAGraph_malloc (nhooks, sizeof (GrB_Index));
    if (nzid == NULL || p == NULL)
    {
        LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
    }
    LAGRAPH_OK (GrB_Vector_extractTuples(nzid, p, &nhooks, hookP));

    //a dense vector of hooks for assignments
    GrB_Vector hook_dense = NULL;
    LAGRAPH_OK (GrB_Vector_new (&hook_dense, GrB_UINT64, nhooks));
    LAGRAPH_OK (GrB_extract (hook_dense, NULL, NULL, hook,  nzid, nhooks, NULL)) ;
    //GxB_Vector_fprint(hook_dense, "---- hook_dense ------", GxB_SHORT, stderr);

    // update grandparents of hooks
    LAGRAPH_OK (GrB_assign (parents, NULL, NULL, hook_dense, p,  nhooks, NULL)) ;
    //GxB_Vector_fprint(parents, "---- parents ------", GxB_SHORT, stderr);
    
    LAGRAPH_FREE (nzid);
    LAGRAPH_FREE (p);
    GrB_free (&hook_dense) ;
    GrB_free (&hookP) ;
    GrB_free (&descNonstars) ;
    GrB_free (&pNonstars) ;
    GrB_free (&hook) ;
    GrB_free (&Min) ;
    GrB_free (&Sel2ndMin) ;
    
    return (GrB_SUCCESS) ;
}


GrB_Info GrandParents(GrB_Vector parents, GrB_Vector grandParents)
{
    GrB_Info info;
    GrB_Index n;
    LAGRAPH_OK (GrB_Vector_size(&n, parents));

    // extract parents for indexing
    GrB_Index *index = LAGraph_malloc (n, sizeof (GrB_Index));
    GrB_Index *p = LAGraph_malloc (n, sizeof (GrB_Index));
    if (index == NULL || p == NULL)
    {
        LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
    }
    LAGRAPH_OK (GrB_Vector_extractTuples(index, p, &n, parents));

    LAGRAPH_OK (GrB_extract (grandParents, NULL, NULL, parents,  p, n, NULL)) ;
    //GxB_Vector_fprint(grandParents, "---- grandParents ------", GxB_SHORT, stderr);

    LAGRAPH_FREE (index);
    LAGRAPH_FREE (p);
    
    return (GrB_SUCCESS) ;
}
GrB_Info Shortcut(GrB_Vector parents)
{
    GrB_Info info;
    GrB_Index n;
    LAGRAPH_OK (GrB_Vector_size(&n, parents));

    //grandparents
    GrB_Vector grandParents = NULL;
    LAGRAPH_OK (GrB_Vector_new (&grandParents, GrB_UINT64, n));
    GrandParents(parents, grandParents);

    // replace parents with grandparents
    LAGRAPH_OK (GrB_assign (parents, NULL, NULL, grandParents, GrB_ALL,  0, NULL)) ;
    
    GrB_free (&grandParents) ;
    return (GrB_SUCCESS) ;
}


GrB_Info StarCheck(GrB_Vector parents, GrB_Vector stars)
{
    GrB_Info info;
    GrB_Index n;
    LAGRAPH_OK (GrB_Vector_size(&n, parents));

    // every vertex is a star
    LAGRAPH_OK (GrB_assign (stars, NULL, NULL, true, GrB_ALL,  0, NULL)) ;

    //grandparents
    GrB_Vector grandParents = NULL;
    LAGRAPH_OK (GrB_Vector_new (&grandParents, GrB_UINT64, n));
    GrandParents(parents, grandParents);
    //GxB_Vector_fprint(grandParents, "---- grandParents ------", GxB_SHORT, stderr);

    //identify vertices whose parents and grandparents are different
    GrB_Vector ns = NULL;
    LAGRAPH_OK (GrB_Vector_new (&ns, GrB_BOOL, n));
    GrB_Vector nsGrandParents = NULL;
    LAGRAPH_OK (GrB_Vector_new (&nsGrandParents, GrB_UINT64, n));
    LAGRAPH_OK (GrB_eWiseMult(ns,NULL, NULL, GrB_NE_UINT64, grandParents,parents, NULL));
    LAGRAPH_OK (GrB_extract (nsGrandParents, ns, NULL, grandParents,  GrB_ALL, 0, NULL)) ;
    //GxB_Vector_fprint(nsGrandParents, "---- p!=gp ------", GxB_SHORT, stderr);


    // extract indices and values for assign
    GrB_Index nNonstars;
    LAGRAPH_OK (GrB_Vector_nvals(&nNonstars, nsGrandParents));
    GrB_Index *vertex = LAGraph_malloc (nNonstars, sizeof (GrB_Index));
    GrB_Index *gp = LAGraph_malloc (nNonstars, sizeof (GrB_Index));
    if (vertex == NULL || gp == NULL)
    {
        LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
    }
    LAGRAPH_OK (GrB_Vector_extractTuples(vertex, gp, &nNonstars, nsGrandParents));

    LAGRAPH_OK (GrB_assign (stars, NULL, NULL, false,  vertex,  nNonstars, NULL)) ;
    LAGRAPH_OK (GrB_assign (stars, NULL, NULL, false,  gp,  nNonstars, NULL)) ;


    // extract indices and values for assign
    GrB_Index *v = LAGraph_malloc (n, sizeof (GrB_Index));
    GrB_Index *p = LAGraph_malloc (n, sizeof (GrB_Index));
    if (v == NULL || p == NULL)
    {
        LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
    }
    LAGRAPH_OK (GrB_Vector_extractTuples(v, p, &n, parents));

    GrB_Vector starsf = NULL;
    LAGRAPH_OK (GrB_Vector_new (&starsf, GrB_BOOL, n));
    LAGRAPH_OK (GrB_extract (starsf, NULL, NULL, stars,  p, n, NULL)) ;
    LAGRAPH_OK (GrB_assign (stars, NULL, NULL, starsf, GrB_ALL,  0, NULL)) ;

    LAGRAPH_FREE (vertex);
    LAGRAPH_FREE (gp);
    LAGRAPH_FREE (v);
    LAGRAPH_FREE (p);
    
    GrB_free (&grandParents) ;
    GrB_free (&ns) ;
    GrB_free (&nsGrandParents) ;
    GrB_free (&starsf) ;
    
    return (GrB_SUCCESS) ;
}



GrB_Info CountCC(GrB_Vector parents, GrB_Index* countcc)
{
    
    GrB_Info info ;
    GrB_Index n;
    LAGRAPH_OK (GrB_Vector_size(&n, parents));

    // extract parents for indexing
    GrB_Index *v = LAGraph_malloc (n, sizeof (GrB_Index));
    GrB_Index *p = LAGraph_malloc (n, sizeof (GrB_Index));
    GrB_Index *temp = LAGraph_malloc (n, sizeof (GrB_Index));
    LAGRAPH_OK (GrB_Vector_extractTuples(v, p, &n, parents));
    
    for (int64_t i = 0 ; i < n ; i++)
    {
        temp[i] = 0;
    }
    
    GrB_Index ncc= (uint64_t)0;
    for (int64_t i = 0 ; i < n ; i++)
    {
        if(temp[p[i]] == 0)
        {
            temp[p[i]] = 1;
            ncc ++;
        }
            
    }
    *countcc = ncc;
    /*
     // This code is giving memory error
     // Hence, I did a pure C implementation above
    GrB_Vector cc = NULL;
    LAGRAPH_OK (GrB_Vector_new (&cc, GrB_UINT64, n));
    LAGRAPH_OK (GrB_assign (cc, NULL, NULL, (uint64_t)1, p,  n, NULL)) ;
   
    GrB_Monoid sum = NULL ;
    LAGRAPH_OK (GrB_Monoid_new (&sum, GrB_PLUS_UINT64, (uint64_t)0)) ;
    GrB_Index ncc= (uint64_t)0;
    LAGRAPH_OK (GrB_reduce (&ncc, NULL, sum, cc, NULL)) ;

    LAGRAPH_FREE (v);
    LAGRAPH_FREE (p);
    GrB_free(&cc);
    GrB_free(&sum);

    printf("number of components: %ld\n", ncc);
     */
    LAGRAPH_FREE (v);
    LAGRAPH_FREE (p);
    LAGRAPH_FREE (temp);
    return (GrB_SUCCESS) ;
}


GrB_Info LAGraph_cc
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

    GrB_Vector stars = NULL ;
    GrB_Vector parents = NULL ;
    LAGRAPH_OK (GrB_Vector_new (&stars, GrB_BOOL, n));
    LAGRAPH_OK (GrB_Vector_new (&parents, GrB_UINT64, n));
    for (int64_t i = 0 ; i < n ; i++)
    {
        LAGRAPH_OK (GrB_Vector_setElement (stars, true, i)) ;
        LAGRAPH_OK (GrB_Vector_setElement (parents, i, i) );
    }

    GrB_Vector pchange = NULL;
    LAGRAPH_OK (GrB_Vector_new (&pchange, GrB_BOOL, n));
    bool change = true ;
    GrB_Monoid Lor = NULL ;
    LAGRAPH_OK (GrB_Monoid_new (&Lor, GrB_LOR, (bool) false)) ;

    GrB_Vector parents1 = NULL ;
    LAGRAPH_OK (GrB_Vector_new (&parents1, GrB_UINT64, n));

    while(change)
    {
        LAGRAPH_OK (GrB_Vector_dup(&parents1, parents));
        CondHook(S, parents, stars);
        StarCheck(parents, stars);
        UnCondHook(S, parents, stars);
        Shortcut(parents);
        StarCheck(parents, stars);
        LAGRAPH_OK (GrB_eWiseMult(pchange,NULL, NULL, GrB_NE_UINT64, parents1,parents, NULL));
        LAGRAPH_OK (GrB_reduce (&change, NULL, Lor, pchange, NULL)) ;
    }

    GrB_Index ncc;
    CountCC(parents, &ncc);
    //printf("number of components: %ld\n", ncc);
    *result = parents;
    
    
    GrB_free(&stars);
    GrB_free(&pchange);
    GrB_free(&Lor);
    GrB_free(&parents1);
    return GrB_SUCCESS;
}


