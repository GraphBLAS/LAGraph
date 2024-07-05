//------------------------------------------------------------------------------
// LAGraph_MaximumMatching: maximum matching between nodes of disjoint sets
//                          in bipartite graphs
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Christina Koutsou, Aristotle University of Thessaloniki

// add paper

//------------------------------------------------------------------------------

// add explanation of paper

#include "LG_internal.h"
#include "LAGraphX.h"

//------------------------------------------------------------------------------
// the Vertex tuple: (parentC, rootC)
//------------------------------------------------------------------------------

typedef struct
{
    uint64_t parentC;
    uint64_t rootC;
} vertex;

// repeat the typedef as a string, to give to GraphBLAS
#define VERTEX_DEFN      \
    "typedef struct "    \
    "{ "                 \
    "uint64_t parentC; " \
    "uint64_t rootC; "   \
    "} "                 \
    "vertex; "

void *initFrontier(vertex *z, void *x, uint64_t i, uint64_t j, const void *y)
{
    z->parentC = i;
    z->rootC = i;
}

#define INIT_FRONTIER_DEFN                                                           \
    "void *initFrontier(vertex *z, void *x, uint64_t i, uint64_t j, const void *y) " \
    "{ "                                                                             \
    "z->parentC = i; "                                                               \
    "z->rootC = i; "                                                                 \
    "} "

void *minparent(vertex *z, vertex *x, vertex *y)
{
    *z = x->parentC < y->parentC ? *x : *y;
}

#define MIN_PARENT_DEFN                                 \
    "void *minparent(vertex *z, vertex *x, vertex *y) " \
    "{ "                                                \
    "*z = x->parentC < y->parentC ? *x : *y; "          \
    "} "

// FIXME: revise GraphBLAS so we can tell it that the select2nd operator
// does not use the 'x' input.
void *select2nd(vertex *z, bool *x, vertex *y)
{
    z->parentC = y->parentC;
    z->rootC = y->rootC;
}

#define SELECT_2ND_DEFN                               \
    "void *select2nd(vertex *z, bool *x, vertex *y) " \
    "{ "                                              \
    "z->parentC = y->parentC; "                       \
    "z->rootC = y->rootC;"                            \
    "} "

void *keepParents(uint64_t *z, vertex *x)
{
    *z = x->parentC;
}

#define KEEP_PARENTS_DEFN                        \
    "void *keepParents(uint64_t *z, vertex *x) " \
    "{ "                                         \
    "*z = x->parentC; "                          \
    "} "

void *keepRoots(uint64_t *z, vertex *x)
{
    *z = x->rootC;
}

#define KEEP_ROOTS_DEFN                        \
    "void *keepRoots(uint64_t *z, vertex *x) " \
    "{ "                                       \
    "*z = x->rootC; "                          \
    "} "

void *buildfCTuples(vertex *z, uint64_t *x, uint64_t i, uint64_t j, const void *y)
{
    z->parentC = i;
    z->rootC = *x;
}

#define BUILT_FC_TUPLES_DEFN                                                              \
    "void *buildfCTuples(vertex *z, uint64_t *x, uint64_t i, uint64_t j, const void *y) " \
    "{ "                                                                                  \
    "z->parentC = i; "                                                                    \
    "z->rootC = *x; "                                                                     \
    "} "

void *vertexTypecast(vertex *z, uint64_t *x)
{
    z->parentC = *x;
    z->rootC = *x;
}

#define VERTEX_TYPECAST_DEFN                        \
    "void *vertexTypecast(vertex *z, uint64_t *x) " \
    "{ "                                            \
    "z->parentC = *x; "                             \
    "z->rootC = *x; "                               \
    "} "

void *setParentsMates(vertex *z, vertex *x, vertex *y)
{
    z->parentC = y->parentC;
    z->rootC = x->rootC;
};

#define SET_PARENTS_MATES_DEFN                                \
    "void *setParentsMates(vertex *z, vertex *x, vertex *y) " \
    "{ "                                                      \
    "z->parentC = y->parentC; "                               \
    "z->rootC = x->rootC; "                                   \
    "} "

//------------------------------------------------------------------------------
// invert_noduplicates
//------------------------------------------------------------------------------

// invert_noduplicates "inverts" an input vector by swapping its row indices
// and its values, returning the result in an output vector.
//
// For example, for the indices/values of an input vector (in) with 5 entries
// and length 100:
//
//      indices: 0  3  5 42 99
//      values:  4 98  1  3 12
//
// on output, the out vector will contain:
//
//      indices: 4 98  1  3 12
//      values:  0  3  5 42 99
//
// The output vector will normally be jumbled since the values will not appear
// in any particular order.  The method assumes that the input values are in
// range 0 to n-1 where n = length(out), and that no values in the input vector
// are duplicated.  Both the in vector and out vector must have the same type
// (GrB_UINT64).  The lengths of the two vectors need not be the same, so long
// as the indices remain in range.  Results are undefined if these conditions
// do not hold.
//
// The in and out vectors may be aliased.  If not aliased, the input vector is
// cleared of all entries on output.  If in and out are aliased, then the
// inversion is performed in-place.
//
// In SuiteSparse:GraphBLAS, this method takes O(1) time if the in vector is in
// CSC (sparse, by column) format.  Otherwise it can take O(e) time if e =
// nvals(in), because the unpack below will convert the in vector to CSC and
// then unpack it.

#undef  LG_FREE_ALL
#define LG_FREE_ALL                         \
    {                                       \
        LAGraph_Free ((void *) &I, NULL) ;  \
        LAGraph_Free ((void *) &X, NULL) ;  \
    }

static GrB_Info invert_noduplicates
(
    GrB_Vector out, // output, contents (except type and size) ignored on input
    GrB_Vector in,  // input vector, empty on output (unless in == out)
    char *msg
)
{
    // the values in the input vector must not contain any duplicates (this is not checked).
    // the output vector will normally be returned in a jumbled state.
    bool jumbled = 1;
    GrB_Index *I = NULL ;
    GrB_Index *X = NULL ;
    GrB_Index IBytes = 0, XBytes = 0;
    uint64_t nvals = 0 ;
    #if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_Vector_unpack_CSC(in, (GrB_Index **)&I, (void **)&X, &IBytes, &XBytes, NULL, &nvals, &jumbled, NULL));
    GRB_TRY(GxB_Vector_pack_CSC(out,  (GrB_Index **)&X, (void **)&I,  XBytes,  IBytes, NULL,  nvals, true, NULL));
    #else
    // vanilla case using extractTuples and build:
    allocate I and X
    GrB_extractTuples (I, X, in, ...)
    GrB_build (out, X, I, ...)
    // free I and X:
    LG_FREE_ALL ;
    #endif
}

//------------------------------------------------------------------------------
// LAGraph_MaximumMatching
//------------------------------------------------------------------------------

#undef  LG_FREE_WORK
#define LG_FREE_WORK                  \
    {                                 \
        GrB_free(&pathC);             \
        GrB_free(&parentsR);          \
        GrB_free(&Vertex);            \
        GrB_free(&frontierC);         \
        GrB_free(&frontierR);         \
        GrB_free(&initFrontierOp);    \
        GrB_free(&I);                 \
        GrB_free(&MinParent);         \
        GrB_free(&MinParent_Monoid);  \
        GrB_free(&Select2ndOp);       \
        GrB_free(&MinParent_2nd_Semiring);          \
        GrB_free(&getParentsOp);      \
        GrB_free(&getRootsOp);        \
        GrB_free(&parentsUpdate);     \
        GrB_free(&ufrontierR);        \
        GrB_free(&mateR);             \
        GrB_free(&rootsufR);          \
        GrB_free(&pathUpdate);        \
        GrB_free(&rootufRIndexes);    \
        GrB_free(&rootsfR);           \
        GrB_free(&rootfRIndexes);     \
        GrB_free(&buildfCTuplesOp);   \
        GrB_free(&vertexTypecastOp);  \
        GrB_free(&setParentsMatesOp); \
        GrB_free(&vr);                \
        GrB_free(&pathCopy);          \
        GrB_free(&currentMatesR);     \
    }

#undef  LG_FREE_ALL
#define LG_FREE_ALL           \
    {                         \
        LG_FREE_WORK;         \
        GrB_free(&mateCcopy); \
    }

int LAGraph_MaximumMatching(
    // output/input:

    GrB_Vector *mateC, // mateC(j) = i : Column j of the C subset is matched to row i of the R subset

    // output:
    // GrB_Vector *mateC_handle,    // ignored on input 
    // input:
    // GrB_Vector mateC_init,       // input only, not modified, ignored if NULL

    // input:
    GrB_Matrix A, // input adjacency matrix, TODO: this should be a LAGraph of a BIPARTITE kind
    char *msg)
{

    // GrB_set (GrB_GLOBAL, true, GxB_BURBLE) ;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Vector pathC = NULL;    // make this bitmap, if dense I would have to give all the entries and make the matrix 1-based
    GrB_Vector parentsR = NULL; // parents of row nodes that are reachable from paths of the initial column frontier
    GrB_Type Vertex = NULL;
    GrB_Vector frontierC = NULL;
    GrB_Vector frontierR = NULL;
    GrB_IndexUnaryOp initFrontierOp = NULL;
    GrB_Vector I = NULL; // dense vector of 1's
    GrB_BinaryOp MinParent = NULL;
    GrB_Monoid MinParent_Monoid = NULL;
    GrB_BinaryOp Select2ndOp = NULL;
    GrB_Semiring MinParent_2nd_Semiring = NULL;
    GrB_UnaryOp getParentsOp = NULL;
    GrB_UnaryOp getRootsOp = NULL;
    GrB_Vector parentsUpdate = NULL; // unmatched rows of R frontier
    GrB_Vector ufrontierR = NULL;    // unmatched rows of R frontier
    GrB_Vector mateR = NULL;         // mateR(i) = j : Row i of the R subset is matched to column j of the C subset
    GrB_Vector rootsufR = NULL;
    GrB_Vector pathUpdate = NULL;
    GrB_Vector rootufRIndexes = NULL;
    GrB_Vector rootsfR = NULL;
    GrB_Vector rootfRIndexes = NULL;
    GrB_IndexUnaryOp buildfCTuplesOp = NULL;
    GrB_UnaryOp vertexTypecastOp = NULL;
    GrB_BinaryOp setParentsMatesOp = NULL;
    GrB_Vector vr = NULL;
    GrB_Vector pathCopy = NULL;
    GrB_Vector currentMatesR = NULL;

    // FIXME: no need for mateCcopy, just use mateC
    GrB_Vector mateCcopy = *mateC;

    LG_CLEAR_MSG;

//  LG_ASSERT (mateC_handle != NULL, GrB_NULL_POINTER) ;
//  LG_ASSERT (A != NULL, GrB_NULL_POINTER) ;
//  (*mateC_handle) = NULL ;

//  GrB_Vector mateC = NULL :

    uint64_t ncols = 0;
    GRB_TRY(GrB_Matrix_ncols(&ncols, A));

    uint64_t nrows = 0;
    GRB_TRY(GrB_Matrix_nrows(&nrows, A));

#if 0
    GRB_TRY (GrB_Vector_new (&mateC, GrB_UINT64, ncols)) ;
    if (mateC_init != NULL)
    {
        // mateC = (uint64_t) mateC_init
        GRB_TRY (GrB_assign (mateC, NULL, NULL, mateC_init, GrB_ALL, ncols, NULL)) ;
    }
#endif

    GRB_TRY(GrB_Vector_new(&pathC, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&parentsR, GrB_UINT64, nrows));

    GRB_TRY(GxB_Type_new(&Vertex, sizeof(vertex), "vertex", VERTEX_DEFN));

    GRB_TRY(GrB_Vector_new(&frontierC, Vertex, ncols));

    GRB_TRY(GrB_Vector_new(&frontierR, Vertex, nrows));

    GRB_TRY(GxB_IndexUnaryOp_new(&initFrontierOp, (void *)initFrontier, Vertex, GrB_BOOL, GrB_BOOL, "initFrontier", INIT_FRONTIER_DEFN));

    GRB_TRY(GrB_Vector_new(&I, GrB_BOOL, ncols));
    GRB_TRY(GrB_Vector_assign_BOOL(I, NULL, NULL, 1, GrB_ALL, ncols, NULL)); // pack with GrB_ALL as indexes?

    GRB_TRY(GxB_BinaryOp_new(&MinParent, (void *)minparent, Vertex, Vertex, Vertex, "minparent", MIN_PARENT_DEFN));
    vertex infinityParent = {GrB_INDEX_MAX + 1, 0};
    GRB_TRY(GrB_Monoid_new_UDT(&MinParent_Monoid, MinParent, &infinityParent));

    GRB_TRY(GxB_BinaryOp_new(&Select2ndOp, (void *)select2nd,
                             Vertex, GrB_BOOL, Vertex, "select2nd", SELECT_2ND_DEFN));

    GRB_TRY(GrB_Semiring_new(&MinParent_2nd_Semiring, MinParent_Monoid, Select2ndOp));

    GRB_TRY(GxB_UnaryOp_new(&getParentsOp, (void *)keepParents, GrB_UINT64, Vertex, "keepParents", KEEP_PARENTS_DEFN));

    GRB_TRY(GxB_UnaryOp_new(&getRootsOp, (void *)keepRoots, GrB_UINT64, Vertex, "keepRoots", KEEP_ROOTS_DEFN));

    GRB_TRY(GrB_Vector_new(&parentsUpdate, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&ufrontierR, Vertex, nrows));

    GRB_TRY(GrB_Vector_new(&mateR, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&rootsufR, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&pathUpdate, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&rootufRIndexes, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&rootsfR, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&rootfRIndexes, GrB_UINT64, ncols));

    GRB_TRY(GxB_IndexUnaryOp_new(&buildfCTuplesOp, (void *)buildfCTuples, Vertex, GrB_UINT64, GrB_BOOL, "buildfCTuples", BUILT_FC_TUPLES_DEFN));

    GRB_TRY(GxB_UnaryOp_new(&vertexTypecastOp, (void *)vertexTypecast, Vertex, GrB_UINT64, "vertexTypecast", VERTEX_TYPECAST_DEFN));

    GRB_TRY(GxB_BinaryOp_new(&setParentsMatesOp, (void *)setParentsMates, Vertex, Vertex, Vertex, "setParentsMates", SET_PARENTS_MATES_DEFN));

    GRB_TRY(GrB_Vector_new(&vr, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&pathCopy, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&currentMatesR, GrB_UINT64, nrows));

    uint64_t npath = 0;
    bool y = 0; // see if I can get rid of this

    // FIXME: 80 characters wide? or at least something small than 243
    uint64_t nmatched = 0;
    GRB_TRY(GrB_Vector_nvals(&nmatched, mateCcopy));
    if (nmatched)
    {

        // FIXME: make this a helper function {
        // mateR = invert (mateC), but do not clear the input
        GrB_Index *J, *X; // unpack allocates space for these lists
        GrB_Index Jbytes = 0, Xbytes = 0;
        bool jumbledMateC = 0;
        GRB_TRY(GxB_Vector_unpack_CSC(mateCcopy, (GrB_Index **)&J, (void **)&X, &Jbytes, &Xbytes, NULL, &nmatched, &jumbledMateC, NULL)); // mateC and mateR do not have duplicates, so the order doesn't matter
        GRB_TRY(GrB_Vector_clear(mateR));                                                                                                 // clear mateR first as a prerequisite of the build method
        GRB_TRY(GrB_Vector_build_UINT64(mateR, X, J, nmatched, NULL));                                                                    // build does not take ownership of the lists J and X, but only copies them,
                                                                                                                                          // these lists will be given again to mateC
                                                                                                                                          // mateC has no duplicates in the values list, so mateR doesn't need to handle dups
        GRB_TRY(GxB_Vector_pack_CSC(mateCcopy, (GrB_Index **)&J, (void **)&X, Jbytes, Xbytes, NULL, nmatched, jumbledMateC, NULL));
        // } end of helper function

    }
    /* debug
    GxB_Vector_fprint(mateR, "mateR", GxB_COMPLETE, stdout);
    */

    do
    {
        GRB_TRY(GrB_Vector_clear(parentsR));
        // for every col j not matched, assign f(j) = VERTEX(j,j)
        GRB_TRY(GrB_Vector_apply_IndexOp_UDT(frontierC, mateCcopy, NULL, initFrontierOp, I, &y, GrB_DESC_RSC));

        /* debug
        GrB_Index C[ncols];
        vertex *V = malloc(ncols * sizeof(vertex));
        GrB_Vector_extractTuples_UDT(C, V, &ncols, frontierC);
        for (int k = 0; k < ncols; k++)
        {
            printf("\nfc (%d) = (%ld, %ld)", (int)C[k], V[k].parentC, V[k].rootC);
        }
        GxB_Vector_fprint(mateCcopy, "mateC", GxB_COMPLETE, stdout);
        */

        uint64_t nfC = 0;

        do
        {
            // STEPS 1,2: Explore neighbors of column frontier (one step of BFS),
            // keeping only unvisited rows in the frontierR
            GRB_TRY(GrB_mxv(frontierR, parentsR, NULL, MinParent_2nd_Semiring, A, frontierC, GrB_DESC_RSC));

            // STEPS 3,4: Select matched and unmatched row vertices

            // set parents of row frontier
            GRB_TRY(GrB_Vector_apply(parentsR, frontierR, NULL, getParentsOp, frontierR, GrB_DESC_S)); // use input as mask to only update or insert parents without deleting the ones not updated

            // select unmatched rows of the R frontier
            GRB_TRY(GrB_Vector_assign(ufrontierR, mateR, NULL, frontierR, GrB_ALL, nrows, GrB_DESC_RSC));
            // select matched rows of the R frontier
            GRB_TRY(GrB_Vector_assign(frontierR, mateR, NULL, frontierR, GrB_ALL, nrows, GrB_DESC_RS));
            // keep only mates of rows in frontierR
            GRB_TRY(GrB_Vector_assign(currentMatesR, frontierR, NULL, mateR, GrB_ALL, nrows, GrB_DESC_RS));

            /* debug
            uint64_t nvals = 0;
            GxB_Vector_fprint(currentMatesR, "currentMatesR", GxB_COMPLETE, stdout);
            GrB_Index *R = (GrB_Index *)malloc(nrows * sizeof(GrB_Index));
            vertex *VR = (vertex *)malloc(nrows * sizeof(vertex));
            GrB_Vector_nvals(&nvals, frontierR);
            GrB_Vector_extractTuples_UDT(R, VR, &nrows, frontierR);
            for (int k = 0; k < nrows; k++)
            {
                printf("\nfr (%d) = (%ld, %ld)", (int)R[k], VR[k].parentC, VR[k].rootC);
            }
            GxB_Vector_fprint(parentsR, "pr", GxB_COMPLETE, stdout);
            */

            uint64_t nUfR = 0, nfR = 0;
            GRB_TRY(GrB_Vector_nvals(&nUfR, ufrontierR));
            GRB_TRY(GrB_Vector_nvals(&nfR, frontierR));

            if (nUfR)
            {
                // STEP 5:

                // get roots of unmatched row nodes in the R frontier
                GRB_TRY(GrB_Vector_apply(rootsufR, NULL, NULL, getRootsOp, ufrontierR, NULL));

                // FIXME: make this a helper function {
                // pathUpdate = invert (rootsufR), but need to handle duplicates
                GrB_Index *IrootsufR, *VrootsufR;
                GrB_Index Ibytes = 0, Valbytes = 0;
                GRB_TRY(GxB_Vector_unpack_CSC(rootsufR, (GrB_Index **)&IrootsufR, (void **)&VrootsufR, &Ibytes, &Valbytes, NULL, &nUfR, NULL, NULL)); // sorted indices so we keep the min child
                GRB_TRY(GrB_Vector_clear(pathUpdate));
                GRB_TRY(GrB_Vector_build_UINT64(pathUpdate, VrootsufR, IrootsufR, nUfR, GrB_FIRST_UINT64));  // useful to handle duplicates
                // build copies the lists so they need to be freed:
                LG_TRY(LAGraph_Free((void **)&IrootsufR, msg));
                LG_TRY(LAGraph_Free((void **)&VrootsufR, msg));
                // } end of helper function

                GRB_TRY(GrB_Vector_assign(pathC, pathUpdate, NULL, pathUpdate, GrB_ALL, ncols, GrB_DESC_S)); // update path without deleting the values not updated
                                                                                                             // when GrB_ALL is used, ni is the number of rows of the vector

                // STEP 6:
                GRB_TRY(GrB_Vector_clear(rootfRIndexes));

                if (nfR)
                {
                    // get roots of row nodes in the current R frontier
                    GRB_TRY(GrB_Vector_apply(rootsfR, NULL, NULL, getRootsOp, frontierR, NULL));

                    /* debug
                    GxB_Vector_fprint(rootsfR, "rootsfR", GxB_COMPLETE, stdout);
                    */

                    // FIXME: make this a helper function {
                    GrB_Index *VmatesfR, *VrootsfR, *dummy;
                    GrB_Index nRootsfR = 0;
                    GrB_Index n_dummy = 1, bytes_dummy = 0;
                    // keep mates of the R frontier (ordered indices)
                    GRB_TRY(GxB_Vector_unpack_CSC(currentMatesR, (GrB_Index **)&dummy, (void **)&VmatesfR, &bytes_dummy, &Valbytes, NULL, &nfR, NULL, NULL));
                    LG_TRY(LAGraph_Free((void **)&dummy, msg));
                    // keep roots of the R frontier (ordered indices)
                    GRB_TRY(GxB_Vector_unpack_CSC(rootsfR, (GrB_Index **)&dummy, (void **)&VrootsfR, &bytes_dummy, &Ibytes, NULL, &nRootsfR, NULL, NULL));
                    LG_TRY(LAGraph_Free((void **)&dummy, msg));
                    GRB_TRY(GrB_Vector_clear(rootfRIndexes));
                    GRB_TRY(GrB_Vector_build_UINT64(rootfRIndexes, VrootsfR, VmatesfR, nRootsfR, GrB_FIRST_UINT64)); // rootfRIndexes(j) = i, where i is the col mate of the first row
                                                                                                                     // included in the current R frontier with a col root of j
                    LG_TRY(LAGraph_Free((void **)&VmatesfR, msg));
                    LG_TRY(LAGraph_Free((void **)&VrootsfR, msg));
                    // } end of helper function

                    // keep only col roots that are not included in ufR
                    GRB_TRY(GrB_Vector_assign(rootfRIndexes, pathUpdate, NULL, rootfRIndexes, GrB_ALL, ncols, GrB_DESC_RSC));

                    // rootfRIndexes = invert (rootfRIndexes), so that
                    // rootfRIndexes(i) = j, where (i,j) = (parentC, rootC) of the new frontier C
                    LAGRAPH_TRY (invert_noduplicates (rootfRIndexes, rootfRIndexes, msg)) ;
                }

                // STEP 7b: when ufrontierR is not empty
                // build tuple of (parentC, rootC)
                GRB_TRY(GrB_Vector_apply_IndexOp_UDT(frontierC, NULL, NULL, buildfCTuplesOp, rootfRIndexes, &y, NULL));

                /* debug
                GrB_Index C[ncols];
                vertex *V = malloc(ncols * sizeof(vertex));
                GrB_Vector_extractTuples_UDT(C, V, &ncols, frontierC);
                for (int k = 0; k < ncols; k++)
                {
                    printf("\nfc (%d) = (%ld, %ld)", (int)C[k], V[k].parentC, V[k].rootC);
                }
                */
            }
            else
            {
                // STEP 7: when ufrontierR is empty

                // typecast mateR to ensure domain match with frontier R and apply op on frontier to set parents to mates
                GRB_TRY(GrB_Vector_apply(frontierR, NULL, setParentsMatesOp, vertexTypecastOp, currentMatesR, NULL)); // fR(i) = (column mate of i, rootC)  // add the structural mask
                // invert fr

                // FIXME: make this a helper function {
                GrB_Index *VmatesfR, *VfR, *dummy;
                GrB_Index bytes_dummy = 0, Vmatesbytes = 0, VfRBytes = 0, nfR = 0;
                GRB_TRY(GxB_Vector_unpack_CSC(currentMatesR, (GrB_Index **)&dummy, (void **)&VmatesfR, &bytes_dummy, &Vmatesbytes, NULL, &nfR, NULL, NULL)); // currentMatesR already contains only the rows of fR
                LG_TRY(LAGraph_Free((void **)&dummy, msg));
                GRB_TRY(GxB_Vector_unpack_CSC(frontierR, (GrB_Index **)&dummy, (void **)&VfR, &bytes_dummy, &VfRBytes, NULL, &nfR, NULL, NULL));
                LG_TRY(LAGraph_Free((void **)&dummy, msg));
                // assign to fC
                GRB_TRY(GxB_Vector_pack_CSC(frontierC, (GrB_Index **)&VmatesfR, (void **)&VfR, Vmatesbytes, VfRBytes, NULL, nfR, true, NULL)); // the values are not ordered,
                                                                                                                                               // so the indices of the inverted fR are jumbled
                // } end of helper function

            }

            GRB_TRY(GrB_Vector_nvals(&nfC, frontierC));

        } while (nfC);

        //----------------------------------------------------------------------
        // STEP 8: Augment matching by all augmenting paths discovered in this phase
        //----------------------------------------------------------------------

        GRB_TRY(GrB_Vector_nvals(&npath, pathC));
        uint64_t npathCopy = npath;
//      GrB_Index *Ipath, *Xpath;
//      GrB_Index IpathBytes = 0, XpathBytes = 0;
        while (npath)
        {
            // vr = invert (pathC), leaving pathC empty
            // pathC doesn't have dup values as it stems from an invertion
            LAGRAPH_TRY (invert_noduplicates (vr, pathC, msg)) ;

            /* debug
            GxB_Vector_fprint(vr, "vr", GxB_COMPLETE, stdout);
            GxB_Vector_fprint(parentsR, "parentsR", GxB_COMPLETE, stdout);
            */

            // assign parents of rows to rows
            GRB_TRY(GrB_Vector_assign(vr, vr, NULL, parentsR, GrB_ALL, nrows, GrB_DESC_S)); // update the values of vr (descriptor needed to use mask's structure and not values)

            /* debug
            GxB_Vector_fprint(vr, "vr with updated parents", GxB_COMPLETE, stdout);
            */

            // update mateR:  mateR<vr> = vr
            GRB_TRY(GrB_Vector_assign(mateR, vr, NULL, vr, GrB_ALL, nrows, GrB_DESC_S));

            // pathC = invert (vr), leaving vr empty (vr has no duplicates)
            LAGRAPH_TRY (invert_noduplicates (pathC, vr, msg)) ;

            /* debug
            GxB_Vector_fprint(pathC, "pathC", GxB_COMPLETE, stdout);
            */

            // keep a copy of the previous row matches of the matched cols that will alter mates
            GRB_TRY(GrB_Vector_assign(pathCopy, pathC, NULL, mateCcopy, GrB_ALL, ncols, GrB_DESC_RS));

            /* debug
            GxB_Vector_fprint(pathCopy, "pathCopy", GxB_COMPLETE, stdout);
            */

            // update mateC
            GRB_TRY(GrB_Vector_assign(mateCcopy, pathC, NULL, pathC, GrB_ALL, ncols, GrB_DESC_S));
            // swap path and pathCopy
            GrB_Vector temp = pathC;
            pathC = pathCopy;
            pathCopy = temp;

            GRB_TRY(GrB_Vector_nvals(&npath, pathC));

            /* debug
            GxB_Vector_fprint(mateCcopy, "mateC", GxB_COMPLETE, stdout);
            */
        }

        npath = npathCopy;
    } while (npath); // only in the first and last iteration should this condition be false

    /* debug
    GxB_Vector_fprint(mateCcopy, "mateC", GxB_COMPLETE, stdout);
    */

    (*mateC) = mateCcopy;
    LG_FREE_WORK;

    // GrB_set (GrB_GLOBAL, false, GxB_BURBLE) ;
//  (*mateC_handle) = mateC ;
    return (GrB_SUCCESS);
}
