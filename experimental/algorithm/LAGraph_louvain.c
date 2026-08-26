//------------------------------------------------------------------------------
// LAGraph_louvain.c: Louvain method
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2026 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Roi Lipman and Gabriel Gomez, FalkorDB

//------------------------------------------------------------------------------

#include <LAGraph.h>
#include <LAGraphX.h>
#include <LG_internal.h>

#undef  LG_FREE_ALL
#define LG_FREE_ALL          \
{                            \
	GrB_free(&V);            \
	GrB_free(&cont);         \
}

// construct C, C[i:] = [i]
// and N, N[i] = i's community id
static GrB_Info _initCommunities
(
	GrB_Matrix *C,        // C [i] = i
	uint64_t **I,         // C's indices array
	GrB_Index node_count  // node count
)
{
	//--------------------------------------------------------------------------
	// initialize C
	//--------------------------------------------------------------------------

	// create full vector V[i] = 1
	char *msg = NULL ;
	GrB_Vector V = NULL ;
	GxB_Container cont = NULL ;

	GRB_TRY (GrB_Vector_new (&V, GrB_BOOL, node_count)) ;
	GRB_TRY (GrB_assign (V, NULL, NULL, true, GrB_ALL, node_count, NULL)) ;

	// create diagonal matrix C from V
	GRB_TRY (GrB_Matrix_diag (C, V, 0)) ;
	GRB_TRY (GrB_free (&V)) ;

	// C should be CSR with 64 bit indicies
	GRB_TRY (GrB_set (*C, 64,           GxB_COLINDEX_INTEGER_HINT   )) ;
	GRB_TRY (GrB_set (*C, GxB_SPARSE,   GxB_SPARSITY_CONTROL        )) ;
	GRB_TRY (GrB_set (*C, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT)) ;

	//--------------------------------------------------------------------------
	// get a handle to C's indicies array
	//--------------------------------------------------------------------------

	// unload C into a container
	GRB_TRY (GxB_Container_new (&cont)) ;

	// GrB_Matrix -> GxB_Container
	GRB_TRY (GxB_unload_Matrix_into_Container (*C, cont, NULL)) ;

	int handling          ;
	GrB_Type type         ;
	uint64_t n, X_memsize ;

	// unload I
	GRB_TRY (GxB_Vector_unload (cont->i, (void**)I, &type, &n, &X_memsize,
				&handling, NULL)) ;

	// load I, mark it as READONLY
	GRB_TRY (GxB_Vector_load (cont->i, (void **)I, type, n, X_memsize,
				GxB_IS_READONLY, NULL)) ;

	// GxB_Container -> GrB_Matrix
	GRB_TRY (GxB_load_Matrix_from_Container (*C, cont, NULL)) ;
	GRB_TRY (GrB_free (&cont)) ;

	LG_FREE_ALL ;

	return GrB_SUCCESS ;
}

#undef  LG_FREE_WORK
#define LG_FREE_WORK                               \
{                                                  \
	GrB_free(&C);                                  \
	GrB_free(&D);                                  \
	GrB_free(&it);                                 \
	GrB_free(&Ni);                                 \
	GrB_free(&desc);                               \
	LAGraph_Free((void**)&degree, NULL);           \
	LAGraph_Free((void**)&community_degree, NULL); \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL      \
{                        \
	LG_FREE_WORK ;       \
}

// computes the modularity contribution of attaching a node of degree
// i_degree, with kin edges into a community, to that community -- where
// sigma_tot is the total degree of the community's members, *excluding*
// the node itself if it happens to already be a member (the caller is
// responsible for that adjustment; see base_gain)
static inline double gain
(
	uint64_t i_degree,   // node i's degree
	uint64_t kjn,        // number of edges connecting node i to community j
	uint64_t sigma_tot,  // total degree of community j members
	double M2            // number of edges * 2
)
{
	return (2.0 * (double) kjn) / M2 -
		((double) sigma_tot * (double) i_degree) / (M2 * M2) ;
}

// compute clustering by running the Louvain algorithm against the graph's
// adjacency matrix A
GrB_Info LAGraph_louvain
(
	GrB_Vector *com,  // output communities
	LAGraph_Graph G,  // graph adjacency matrix
	int itermax,      // max number of modularity improvements sweeps per level
	int levelmax,     // max number of modularity improve and cluster condense
	float e,          // min change in modularity considered an improvement
	char *msg         // error message
)
{
	//--------------------------------------------------------------------------
	// check inputs
	//--------------------------------------------------------------------------

	if (com == NULL || G == NULL || msg == NULL) {
		return (GrB_NULL_POINTER) ;
	}

	GrB_Matrix C           = NULL ;  // map between node to community id
	GrB_Vector D           = NULL ;  // map between nodes to community ids
	GxB_Iterator it        = NULL ;
	GrB_Vector Ni          = NULL ;  // node i's neighbors
    GrB_Descriptor desc    = NULL ;

	uint64_t *degree           = NULL ;  // node degree
	uint64_t *community_degree = NULL ;  // sigma tot

    // find out if graph is symmetric, compute G->out_degree, and G->nself_edges
    LG_TRY (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;
    LG_TRY (LAGraph_Cached_OutDegree (G, msg)) ;

    LG_TRY (LAGraph_Cached_NSelfEdges (G, msg)) ;
    LG_ASSERT_MSG (G->nself_edges == 0, GrB_INVALID_VALUE,
			"G->nself_edges must be zero") ;

	GrB_Matrix A = G->A ;

	// TODO: make sure `A` is row-wise

	GrB_Index nrows ;
	GrB_Index ncols ;
	GRB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
	GRB_TRY (GrB_Matrix_ncols (&ncols, A)) ;

	// expecting a square matrix
	LG_TRY (LAGraph_CheckGraph (G, msg)) ;
	LG_ASSERT_MSG (nrows == ncols, LAGRAPH_INVALID_GRAPH,
            "adjacency matrix must be square") ;

	GrB_Type t ;
	GRB_TRY (GxB_Matrix_type  (&t, A)) ;
	LG_ASSERT_MSG (t == GrB_BOOL, LAGRAPH_INVALID_GRAPH, "A must be boolean") ;

	GrB_Index  nvals ;
	GRB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
	double M2 = (double) nvals ;

	//--------------------------------------------------------------------------
	// compute nodes degree
	//--------------------------------------------------------------------------

	GRB_TRY (GrB_Vector_new (&D, GrB_UINT64, nrows)) ;
	GRB_TRY (GrB_assign (D, NULL, NULL, 0, GrB_ALL, nrows, NULL)) ;
	GRB_TRY (GrB_assign (D, NULL, GrB_PLUS_UINT64, G->out_degree, GrB_ALL,
				nrows, NULL)) ;

	// TODO: enable once we support weighted graphs
	//GRB_TRY (GrB_mxv (D, NULL, GrB_PLUS_UINT64, GxB_PLUS_PAIR_UINT64, A, D,
				//NULL)) ;

	// unpack D to an array for direct access
	uint64_t degree_size ;
	GRB_TRY (GxB_Vector_unpack_Full (D, (void**)&degree, &degree_size, NULL,
				NULL)) ;
	GRB_TRY (GrB_free (&D)) ;

	// community_degree [i] == degree [i]
	// as each node is in its own community
    LG_TRY (LAGraph_Malloc ((void **) &community_degree, 1, degree_size, msg)) ;
	memcpy (community_degree, degree, degree_size) ;

	// initialize C
	uint64_t *I = NULL ;
	GRB_TRY   (_initCommunities (&C, &I, nrows)) ;
	LG_ASSERT (I != NULL, GrB_NULL_POINTER)  ;

	GRB_TRY (GrB_Vector_new   (&Ni, GrB_UINT64, ncols)) ;
    GRB_TRY (GxB_Iterator_new (&it)) ;

	GRB_TRY (GrB_Descriptor_new (&desc)) ;
	GRB_TRY (GrB_set (desc, GrB_STRUCTURE,   GrB_MASK_FIELD)) ;
	GRB_TRY (GrB_set (desc, GxB_USE_INDICES, GxB_ROWINDEX_LIST)) ;

	bool improved ;
	double modularity_gain ;

	// start sweep
	do
	{
		modularity_gain = 0 ;

		// for each node
		for (GrB_Index i = 0 ; i < nrows ; i++) {
			// determine i's current community
			uint64_t ic = I [i] ; // i's community

			//------------------------------------------------------------------
			// get a set of communities i can migrate to
			//------------------------------------------------------------------

			// get i's neighbors; A[i:]
			GRB_TRY (GrB_Col_extract (Ni, NULL, NULL, A, GrB_ALL, ncols,
						i, GrB_DESC_T0)) ;

			// get neighbors communities
			// Ni * C = X[i]=j i community ID, j #neighbors in the ith community
			GRB_TRY (GrB_vxm (Ni, NULL, NULL, GxB_PLUS_PAIR_UINT64, Ni, C, NULL)) ;

			GRB_TRY (GxB_Vector_Iterator_attach (it, Ni, NULL)) ;
			GrB_Info info = GxB_Vector_Iterator_seek (it, 0) ;

			//------------------------------------------------------------------
			// compute base score removing i from its current community
			//------------------------------------------------------------------

			// i's degree
			uint64_t i_degree = degree [i] ;

			// number of edges connecting i to its current community
			uint64_t kin = 0 ;
			GRB_TRY (GrB_Vector_extractElement (&kin, Ni, ic)) ;

			// adjusted sigma tot
			uint64_t sigma_tot = community_degree [ic] ;
			//LG_ASSERT (sigma_tot >= i_degree) ;
			sigma_tot -= i_degree ;

			double   base_gain  = gain (i_degree, kin, sigma_tot, M2) ;
			double   max_modularity = base_gain ;
			uint64_t best_community = ic ;

			while (info != GxB_EXHAUSTED)
			{
				// candidate community
				GrB_Index jc = GxB_Vector_Iterator_getIndex (it) ;

				// number of edges from node i to candidate community
				uint64_t kjn = GxB_Iterator_get_UINT64 (it) ;

				// move to the next entry in Ni
				info = GxB_Vector_Iterator_next (it) ;

				// skip i's community
				if (jc == ic)
				{
					continue ;
				}

				// total degree of nodes in community
				sigma_tot = community_degree [jc] ;

				// modularity gain of i joining the candidate community
				double g = gain (i_degree, kjn, sigma_tot, M2) ;

				if (g > max_modularity)
				{
					best_community = jc ;
					max_modularity = g  ;
				}
			}

			if (best_community != ic) {
				// migrate i to its new community j
				I [i] = best_community ;

				//--------------------------------------------------------------
				// update community degree
				// community_degree [ic]             -= degree [i]
				// community_degree [best_community] += degree [i]
				//--------------------------------------------------------------

				community_degree [ic]             -= i_degree ;
				community_degree [best_community] += i_degree ;

				// accumulate modularity
				// TODO: at the end of the sweep compute modularity from A and C
				// compare that againt the initial modularity before the sweep
				modularity_gain += max_modularity - base_gain ;
			}
		}

		improved = (modularity_gain > e) ;
	} while (improved && --itermax > 0) ;


	for (uint l = 0 ; l < levelmax ; l++)
	{
		//----------------------------------------------------------------------
		// pick a representative for each community
		//----------------------------------------------------------------------

		GrB_Vector R ;  // representatives R [i] = n; i community ID, n node ID
		GRB_TRY (GrB_Vector_new (&R, GrB_UINT64, nrows)) ;
		GRB_TRY (GrB_Matrix_reduce_Monoid (R, NULL, NULL, GrB_MAX_MONOID_UINT64,
					C, GrB_DESC_T0)) ;

		//----------------------------------------------------------------------
		// map representative to its members
		//----------------------------------------------------------------------

		// MAP [:i] nodes represented by i
		// MAP = MAP * (C * Rdiag)

		//----------------------------------------------------------------------
		// compute A
		//----------------------------------------------------------------------

		// A = CT * C; A[i,j] = x community i is connected to community j with x
		// different connections
	}

	//--------------------------------------------------------------------------
	// set output
	//--------------------------------------------------------------------------

	GRB_TRY (GrB_Vector_new  (com, GrB_UINT64, nrows)) ;
	GRB_TRY (GxB_Vector_load (*com, (void **)(&I), GrB_UINT64, nrows,
				sizeof (uint64_t) * nrows, GrB_DEFAULT, NULL)) ;

	LG_FREE_WORK ;

	return GrB_SUCCESS ;
}

