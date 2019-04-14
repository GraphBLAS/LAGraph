#include "LAGraph.h"

/**
 * Code is based on the algorithm described in the following paper
 * Azad, Buluç. LACC: a linear-algebraic algorithm for finding connected components in distributed memory (IPDPS 2019)
 **/

void CondHook(GrB_Matrix A, GrB_Vector parents, GrB_Vector stars)
{
    GrB_Index n ;
    GrB_Matrix_nrows (&n, A) ;
    
    // Create (Sel2nd,Min)semiring
    GrB_Monoid Min = NULL ;                // Min monoid
    GrB_Semiring Sel2ndMin = NULL ;          // (Sel2nd,Min)semiring
    GrB_Monoid_new (&Min, GrB_MIN_UINT64, (uint64_t)UINT_MAX) ;
    GrB_Semiring_new (&Sel2ndMin, Min, GrB_SECOND_UINT64) ;  
    
    // identify minNeighborparent for star vertices
    // hook stores minNeighborparent
    GrB_Vector hook = NULL;
    GrB_Vector_new (&hook, GrB_UINT64, n);
    GrB_mxv (hook, stars, NULL, Sel2ndMin,  A, parents, NULL) ;
    //GxB_Vector_fprint(hook, "--Hook --", GxB_SHORT, stderr);
    
    // only keep vertices whose minNeighborparent is smaller than its own parent
    GrB_eWiseMult(hook,NULL, NULL, GrB_MIN_UINT64, hook,parents, NULL);
    //GxB_Vector_fprint(hook, "---- Hook ------", GxB_SHORT, stderr);
    
    //parents of hooks
    GrB_Vector hookP = NULL;
    GrB_Vector_new (&hookP, GrB_UINT64, n);
    GrB_eWiseMult(hookP,NULL, NULL, GrB_SECOND_UINT64, hook,parents, NULL);
    //GxB_Vector_fprint(hookP, "---- hookP ------", GxB_SHORT, stderr);
    
    // extract values in hookP
    GrB_Index nhooks;
    GrB_Vector_nvals(&nhooks, hook);
    GrB_Index *nzid = LAGraph_malloc (nhooks, sizeof (GrB_Index));
    GrB_Index *p = LAGraph_malloc (nhooks, sizeof (GrB_Index));
    GrB_Vector_extractTuples(nzid, p, &nhooks, hookP);
    
    //a dense vector of hooks for assignments
    GrB_Vector hook_dense = NULL;
    GrB_Vector_new (&hook_dense, GrB_UINT64, nhooks);
    GrB_extract (hook_dense, NULL, NULL, hook,  nzid, nhooks, NULL) ;
    //GxB_Vector_fprint(hook_dense, "---- hook_dense ------", GxB_SHORT, stderr);
    
    // update grandparents of hooks
    GrB_assign (parents, NULL, NULL, hook_dense, p,  nhooks, NULL) ;
    //GxB_Vector_fprint(parents, "---- parents ------", GxB_SHORT, stderr);
}


void UnCondHook(GrB_Matrix A, GrB_Vector parents, GrB_Vector stars)
{
    GrB_Index n ;
    GrB_Matrix_nrows (&n, A) ;
    
    // Create (Sel2nd,Min)semiring
    GrB_Monoid Min = NULL ;                // Min monoid
    GrB_Semiring Sel2ndMin = NULL ;          // (Sel2nd,Min)semiring
    GrB_Monoid_new (&Min, GrB_MIN_UINT64, (uint64_t)UINT_MAX) ;
    GrB_Semiring_new (&Sel2ndMin, Min, GrB_SECOND_UINT64) ;
    
    
    // extract parents of nonstar vertices
    GrB_Vector pNonstars = NULL ;
    GrB_Vector_new (&pNonstars, GrB_UINT64, n);
    GrB_Descriptor descNonstars = NULL ;
    GrB_Descriptor_new (&descNonstars) ;
    GrB_Descriptor_set (descNonstars, GrB_MASK, GrB_SCMP) ;     // invert the mask
    GrB_extract (pNonstars, stars, NULL, parents,  GrB_ALL, 0, descNonstars) ;
    //GxB_Vector_fprint(pNonstars, "--pNonstars--", GxB_SHORT, stderr);
    
    
    // identify minNeighborparent for star vertices
    // hook stores minNeighborparent
    GrB_Vector hook = NULL;
    GrB_Vector_new (&hook, GrB_UINT64, n);
    GrB_mxv (hook, stars, NULL, Sel2ndMin,  A, pNonstars, NULL) ;
    //GxB_Vector_fprint(hook, "--Hook --", GxB_SHORT, stderr);
    
    //parents of hooks
    GrB_Vector hookP = NULL;
    GrB_Vector_new (&hookP, GrB_UINT64, n);
    GrB_eWiseMult(hookP,NULL, NULL, GrB_SECOND_UINT64, hook,parents, NULL);
    //GxB_Vector_fprint(hookP, "---- hookP ------", GxB_SHORT, stderr);
    
    
    // exatrct values in hookP
    GrB_Index nhooks;
    GrB_Vector_nvals(&nhooks, hook);
    GrB_Index *nzid = LAGraph_malloc (nhooks, sizeof (GrB_Index));
    GrB_Index *p = LAGraph_malloc (nhooks, sizeof (GrB_Index));
    GrB_Vector_extractTuples(nzid, p, &nhooks, hookP);
    
    //a dense vector of hooks for assignments
    GrB_Vector hook_dense = NULL;
    GrB_Vector_new (&hook_dense, GrB_UINT64, nhooks);
    GrB_extract (hook_dense, NULL, NULL, hook,  nzid, nhooks, NULL) ;
    //GxB_Vector_fprint(hook_dense, "---- hook_dense ------", GxB_SHORT, stderr);
    
    // update grandparents of hooks
    GrB_assign (parents, NULL, NULL, hook_dense, p,  nhooks, NULL) ;
    //GxB_Vector_fprint(parents, "---- parents ------", GxB_SHORT, stderr);
}


void GrandParents(GrB_Vector parents, GrB_Vector grandParents)
{
    GrB_Index n;
    GrB_Vector_size(&n, parents);
    
    // extract parents for indexing
    GrB_Index *index = LAGraph_malloc (n, sizeof (GrB_Index));
    GrB_Index *p = LAGraph_malloc (n, sizeof (GrB_Index));
    GrB_Vector_extractTuples(index, p, &n, parents);
    
    GrB_extract (grandParents, NULL, NULL, parents,  p, n, NULL) ;
    //GxB_Vector_fprint(grandParents, "---- grandParents ------", GxB_SHORT, stderr);
    
    LAGRAPH_FREE (index);
    LAGRAPH_FREE (p);
}
void Shortcut(GrB_Vector parents)
{
    GrB_Index n;
    GrB_Vector_size(&n, parents);
    
    //grandparents
    GrB_Vector grandParents = NULL;
    GrB_Vector_new (&grandParents, GrB_UINT64, n);
    GrandParents(parents, grandParents);
    //GxB_Vector_fprint(grandParents, "---- grandParents ------", GxB_SHORT, stderr);
    
    // replace parents with grandparents
    GrB_assign (parents, NULL, NULL, grandParents, GrB_ALL,  0, NULL) ;
    //GxB_Vector_fprint(parents, "---- parents ------", GxB_SHORT, stderr);
}


void StarCheck(GrB_Vector parents, GrB_Vector stars)
{
    GrB_Index n;
    GrB_Vector_size(&n, parents);
    
    // every vertex is a star
    GrB_assign (stars, NULL, NULL, true, GrB_ALL,  0, NULL) ;
    
    //grandparents
    GrB_Vector grandParents = NULL;
    GrB_Vector_new (&grandParents, GrB_UINT64, n);
    GrandParents(parents, grandParents);
    //GxB_Vector_fprint(grandParents, "---- grandParents ------", GxB_SHORT, stderr);
    
    //identify vertices whose parents and grandparents are different
    GrB_Vector ns = NULL;
    GrB_Vector_new (&ns, GrB_BOOL, n);
    GrB_Vector nsGrandParents = NULL;
    GrB_Vector_new (&nsGrandParents, GrB_UINT64, n);
    GrB_eWiseMult(ns,NULL, NULL, GrB_NE_UINT64, grandParents,parents, NULL);
    GrB_extract (nsGrandParents, ns, NULL, grandParents,  GrB_ALL, 0, NULL) ;
    //GxB_Vector_fprint(nsGrandParents, "---- p!=gp ------", GxB_SHORT, stderr);
    
    
    // extract indices and values for assign
    GrB_Index nNonstars;
    GrB_Vector_nvals(&nNonstars, nsGrandParents);
    GrB_Index *vertex = LAGraph_malloc (nNonstars, sizeof (GrB_Index));
    GrB_Index *gp = LAGraph_malloc (nNonstars, sizeof (GrB_Index));
    GrB_Vector_extractTuples(vertex, gp, &nNonstars, nsGrandParents);
    
    GrB_assign (stars, NULL, NULL, false,  vertex,  nNonstars, NULL) ;
    GrB_assign (stars, NULL, NULL, false,  gp,  nNonstars, NULL) ;
      
    
    // extract indices and values for assign
    GrB_Index *v = LAGraph_malloc (n, sizeof (GrB_Index));
    GrB_Index *p = LAGraph_malloc (n, sizeof (GrB_Index));
    GrB_Vector_extractTuples(v, p, &n, parents);
    
    
    GrB_Vector starsf = NULL;
    GrB_Vector_new (&starsf, GrB_BOOL, n);
    GrB_extract (starsf, NULL, NULL, stars,  p, n, NULL) ;
    GrB_assign (stars, NULL, NULL, starsf, GrB_ALL,  0, NULL) ;
    
    LAGRAPH_FREE (vertex);
    LAGRAPH_FREE (gp);
    LAGRAPH_FREE (v);
    LAGRAPH_FREE (p);
}



void CountCC(GrB_Vector parents)
{
    GrB_Index n;
    GrB_Vector_size(&n, parents);
    // exatrct parents for indexing
    GrB_Index *v = LAGraph_malloc (n, sizeof (GrB_Index));
    GrB_Index *p = LAGraph_malloc (n, sizeof (GrB_Index));
    GrB_Vector_extractTuples(v, p, &n, parents);
    GrB_Vector cc = NULL;
    GrB_Vector_new (&cc, GrB_UINT64, n);
    
    GrB_assign (cc, NULL, NULL, (uint64_t)1, p,  n, NULL) ;
    
    //GxB_Vector_fprint(cc, "---- cc ------", GxB_SHORT, stderr);
    
    GrB_Monoid sum = NULL ;
    GrB_Monoid_new (&sum, GrB_PLUS_UINT64, (uint64_t)0) ;
    GrB_Index ncc= (uint64_t)0;
    GrB_reduce (&ncc, NULL, sum, cc, NULL) ;
    
    LAGRAPH_FREE (v);
    LAGRAPH_FREE (p);
    
    printf("Number of clusters: %ld\n", ncc);
}

// TODO: in progress (nothing returned yet...)
GrB_Info LAGraph_lacc(GrB_Matrix A)
{
    double tic [2], t ;
    GrB_Info info;    
    
    GrB_Index n ;
    GrB_Index nnz ;
    GrB_Matrix_nrows (&n, A) ;
    GrB_Matrix_nvals (&nnz, A) ;
    printf ("number of nodes: %g\n", (double) n) ;
    printf ("number of edges: %g\n", (double) nnz) ;
    
    GrB_Vector stars = NULL ;
    GrB_Vector parents = NULL ;
    GrB_Vector_new (&stars, GrB_BOOL, n);
    GrB_Vector_new (&parents, GrB_UINT64, n);
    for (int64_t i = 0 ; i < n ; i++)
    {
        GrB_Vector_setElement (stars, true, i) ;
        GrB_Vector_setElement (parents, i, i) ;
    }
    
    GrB_Vector pchange = NULL;
    GrB_Vector_new (&pchange, GrB_BOOL, n);
    bool change = true ;
    GrB_Monoid Lor = NULL ;
    GrB_Monoid_new (&Lor, GrB_LOR, (bool) false) ;
    
    GrB_Vector parents1 = NULL ;
    GrB_Vector_new (&parents1, GrB_UINT64, n);
    
    while(change)
    {
        GrB_Vector_dup(&parents1, parents);
        CondHook(A, parents, stars);
        StarCheck(parents, stars);
        UnCondHook(A, parents, stars);
        Shortcut(parents);
        
        
        GrB_eWiseMult(pchange,NULL, NULL, GrB_NE_UINT64, parents1,parents, NULL);
        GrB_reduce (&change, NULL, Lor, pchange, NULL) ;
        //GxB_Vector_fprint(parents1, "---- parents1 ------", GxB_SHORT, stderr);
        //GxB_Vector_fprint(parents, "---- parents ------", GxB_SHORT, stderr);
        //GxB_Vector_fprint(pchange, "---- change ------", GxB_SHORT, stderr);
    }
    
    CountCC(parents);

    return (GrB_SUCCESS) ;
}


