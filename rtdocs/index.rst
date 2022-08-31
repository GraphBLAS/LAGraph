LAGraph Documentation
=====================

The LAGraph library is a collection of high level graph algorithms
based on the GraphBLAS C API.  These algorithms construct
graph algorithms expressed *in the language of linear algebra*.
Graphs are expressed as matrices, and the operations over
these matrices are generalized through the use of a
semiring algebraic structure.

LAGraph is available at https://github.com/GraphBLAS/LAGraph.
LAGraph requires SuiteSparse:GraphBLAS, available at 
https://github.com/DrTimothyAldenDavis/GraphBLAS.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   core
   algorithms
   utils
   experimental
   installation


Example Usage
-------------

.. code-block:: C

    #include "LAGraph.h"
    #define LAGRAPH_FREE_ALL            \
        /* free everything */           \
        LAGraph_Delete (&G, msg) ;      \
        GrB_free (&centrality) ;

    int main (void)
    {
        // initialize LAGraph
        LAGraph_Init (msg) ;
        GrB_Matrix A = NULL ;
        GrB_Vector centrality = NULL ;
        int niters = 0 ;

        // create the karate graph and compute the out-degree of all nodes
        FILE *f = fopen ("LAGraph/data/karate.mtx", "r") ;
        LAGRAPH_TRY (LAGraph_MMRead (&A, f, msg)) ;
        fclose (f) ;
        LAGRAPH_TRY (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
        LAGRAPH_TRY (LAGraph_Cached_OutDegree (G, msg)) ;

        // compute the pagerank
        LAGRAPH_TRY (LAGr_PageRank (&centrality, &niters, G, 0.85, 1e-4, 100, msg)) ;

        // print the result
        LAGRAPH_TRY (LAGraph_Vector_print (centrality, LAGRAPH_COMPLETE, stdout, msg)) ;

        // free the graph, the pagerank, and finish LAGraph
        LAGRAPH_FREE_ALL ;
        LAGraph_Finalize (msg) ;
    }


:ref:`genindex`
