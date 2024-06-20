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

// #define LG_FREE_WORK \
//     {
// }

// #define LG_FREE_ALL \
//     {               \
//         // LG_FREE_WORK ;                  \
//     // GrB_free (centrality) ;
// }

#include "LG_internal.h"

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

void *setParentsMates(vertex *z, uint64_t *x)
{
    z->parentC = *x;
};

#define SET_PARENTS_MATES_DEFN                       \
    "void *setParentsMates(vertex *z, uint64_t *x) " \
    "{ "                                             \
    "z->parentC = *x; "                              \
    "} "

int LAGraph_MaximumMatching(
    // output/input:
    GrB_Vector *mateC, // mateC(j) = i : Column j of the C subset is matched to row i of the R subset
    // input:
    LAGraph_Graph G, // input graph
    char *msg)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    LG_TRY(LAGraph_CheckGraph(G, msg));

    GrB_Matrix A = G->A;
    uint64_t ncols = 0;
    GRB_TRY(GrB_Matrix_ncols(&ncols, A));

    uint64_t nrows = 0;
    GRB_TRY(GrB_Matrix_nrows(&nrows, A));

    GrB_Vector pathC = NULL; // make this bitmap, if dense I would have to give all the entries and make the matrix 1-based
    GRB_TRY(GrB_Vector_new(&pathC, GrB_UINT64, ncols));

    GrB_Vector parentsR = NULL; // parents of row nodes that are reachable from paths of the initial column frontier
    GRB_TRY(GrB_Vector_new(&parentsR, GrB_UINT64, nrows));

    GrB_Type Vertex = NULL;
    GRB_TRY(GxB_Type_new(&Vertex, sizeof(vertex), "vertex", VERTEX_DEFN));

    GrB_Vector frontierC = NULL;
    GRB_TRY(GrB_Vector_new(&frontierC, Vertex, ncols));
    GrB_Vector frontierR = NULL;
    GRB_TRY(GrB_Vector_new(&frontierR, Vertex, nrows));

    GrB_IndexUnaryOp initFrontierOp = NULL;
    GRB_TRY(GxB_IndexUnaryOp_new(&initFrontierOp, (void *)initFrontier, Vertex, GrB_BOOL, GrB_BOOL, "initFrontier", INIT_FRONTIER_DEFN));

    uint64_t nvals = 0;
    bool y = 0; // see if I can get rid of this

    GrB_Vector I = NULL; // dense vector of 1's
    GRB_TRY(GrB_Vector_new(&I, GrB_BOOL, ncols));
    GRB_TRY(GrB_Vector_assign_INT32(I, NULL, NULL, 1, GrB_ALL, ncols, NULL));

    GrB_BinaryOp MinParent = NULL;
    GRB_TRY(GxB_BinaryOp_new(&MinParent, (void *)minparent, Vertex, Vertex, Vertex, "minparent", MIN_PARENT_DEFN));
    vertex infinityParent = {GrB_INDEX_MAX + 1, 0};
    GrB_Monoid AddMonoid;
    GRB_TRY(GrB_Monoid_new_UDT(&AddMonoid, MinParent, &infinityParent));

    GrB_BinaryOp MultMonoid;
    GRB_TRY(GxB_BinaryOp_new(&MultMonoid, (void *)select2nd,
                             Vertex, GrB_BOOL, Vertex, "select2nd", SELECT_2ND_DEFN));
    GrB_Semiring semiring = NULL;
    GRB_TRY(GrB_Semiring_new(&semiring, AddMonoid, MultMonoid));

    GrB_UnaryOp getParentsOp = NULL;
    GRB_TRY(GxB_UnaryOp_new(&getParentsOp, (void *)keepParents, GrB_UINT64, Vertex, "keepParents", KEEP_PARENTS_DEFN));

    GrB_UnaryOp getRootsOp = NULL;
    GRB_TRY(GxB_UnaryOp_new(&getRootsOp, (void *)keepRoots, GrB_UINT64, Vertex, "keepRoots", KEEP_ROOTS_DEFN));

    GrB_Vector ufrontierR = NULL; // unmatched rows of R frontier
    GRB_TRY(GrB_Vector_new(&ufrontierR, Vertex, nrows));

    GrB_Vector mateR = NULL; // mateR(i) = j : Row i of the R subset is matched to column j of the C subset
    GRB_TRY(GrB_Vector_new(&mateR, GrB_UINT64, nrows));

    GrB_Vector rootsufR = NULL;
    GRB_TRY(GrB_Vector_new(&rootsufR, GrB_UINT64, nrows));

    GrB_Vector pathUpdate = NULL;
    GRB_TRY(GrB_Vector_new(&pathUpdate, GrB_UINT64, ncols));

    GrB_Vector rootufRIndexes = NULL;
    GRB_TRY(GrB_Vector_new(&rootufRIndexes, GrB_UINT64, ncols));

    GrB_Vector rootsfR = NULL;
    GRB_TRY(GrB_Vector_new(&rootsfR, GrB_UINT64, nrows));

    GrB_Vector rootfRIndexes = NULL;
    GRB_TRY(GrB_Vector_new(&rootfRIndexes, GrB_UINT64, ncols));

    GrB_IndexUnaryOp buildfCTuplesOp = NULL;
    GRB_TRY(GxB_IndexUnaryOp_new(&buildfCTuplesOp, (void *)buildfCTuples, Vertex, GrB_UINT64, GrB_BOOL, "buildfCTuples", BUILT_FC_TUPLES_DEFN));

    GrB_UnaryOp setParentsMatesOp = NULL;
    GRB_TRY(GxB_UnaryOp_new(&setParentsMatesOp, (void *)setParentsMates, Vertex, GrB_UINT64, "setParentsMates", SET_PARENTS_MATES_DEFN));

    do
    {
        GRB_TRY(GrB_Vector_clear(pathC));
        GRB_TRY(GrB_Vector_clear(parentsR));
        // for every col j not matched, assign f(j) = VERTEX(j,j)
        GRB_TRY(GrB_Vector_apply_IndexOp_UDT(frontierC, *mateC, NULL, initFrontierOp, I, &y, GrB_DESC_RSC));

        /* debug
         GrB_Index C[ncols];
         vertex *V = malloc(ncols * sizeof(vertex));
         GrB_Vector_extractTuples_UDT(C, V, &ncols, frontierC);
         for (int k = 0; k < ncols; k++)
         {
             printf("\nfc (%d) = (%ld, %ld)", (int)C[k], V[k].parentC, V[k].rootC);
         }
        GxB_Vector_fprint(*mateC, "mateC", GxB_COMPLETE, stdout);
        */

        uint64_t nmatched = 0;
        GRB_TRY(GrB_Vector_nvals(&nmatched, *mateC));
        GrB_Index *J = (GrB_Index *)malloc(nmatched * sizeof(GrB_Index));
        GrB_Index *X = (GrB_Index *)malloc(nmatched * sizeof(GrB_Index));
        GrB_Index Jbytes = 0, Xbytes = 0;
        GRB_TRY(GxB_Vector_unpack_CSC(*mateC, (GrB_Index **)&J, (void **)&X, &Jbytes, &Xbytes, NULL, &nmatched, NULL, NULL));
        GRB_TRY(GrB_Vector_clear(mateR));                                          // clear mateR first as a prerequisite of the build method
        GRB_TRY(GrB_Vector_build_UINT64(mateR, X, J, nmatched, GrB_FIRST_UINT64)); // build does not take ownership of the lists J and X, but only copies them
        GRB_TRY(GxB_Vector_pack_CSC(*mateC, (GrB_Index **)&J, (void **)&X, Jbytes, Xbytes, NULL, nmatched, NULL, NULL));
        free(J);
        free(X);

        /* debug
        GxB_Vector_fprint(mateR, "mateR", GxB_COMPLETE, stdout);
        */

        do
        {
            // perform one step of BFS from C nodes and keep only unvisited rows
            GRB_TRY(GrB_mxv(frontierR, parentsR, NULL, semiring, A, frontierC, GrB_DESC_RSC));
            // set parents of row frontier // does it erase the previous values
            GRB_TRY(GrB_Vector_apply(parentsR, NULL, NULL, getParentsOp, frontierR, NULL));

            // select unmatched rows of the R frontier
            GRB_TRY(GrB_Vector_assign(ufrontierR, mateR, NULL, frontierR, GrB_ALL, nrows, GrB_DESC_RSC));
            // select matched rows of the R frontier
            GRB_TRY(GrB_Vector_assign(frontierR, mateR, NULL, frontierR, GrB_ALL, nrows, GrB_DESC_RS));
            // keep only mates of rows in frontierR
            GRB_TRY(GrB_Vector_assign(mateR, frontierR, NULL, mateR, GrB_ALL, nrows, GrB_DESC_RS));

            /* debug
            GrB_Index R[nrows];
            vertex *VR = malloc(nrows * sizeof(vertex));
            GrB_Vector_extractTuples_UDT(R, VR, &nrows, frontierR);
            for (int k = 0; k < nrows; k++)
            {
                printf("\nfr (%d) = (%ld, %ld)", (int)R[k], VR[k].parentC, VR[k].rootC);
            }
            GxB_Vector_fprint(parentsR, "pr", GxB_COMPLETE, stdout);
            */
            uint64_t nUfR = 0;
            GRB_TRY(GrB_Vector_nvals(&nUfR, ufrontierR));

            if (nUfR)
            {
                // get roots of unmatched row nodes in the R frontier
                GRB_TRY(GrB_Vector_apply(rootsufR, NULL, NULL, getRootsOp, ufrontierR, NULL));

                GrB_Index *IrootsufR = (GrB_Index *)malloc(nUfR * sizeof(GrB_Index));
                GrB_Index *VrootsufR = (GrB_Index *)malloc(nUfR * sizeof(GrB_Index));
                GrB_Index Ibytes = 0, Valbytes = 0;
                GRB_TRY(GxB_Vector_unpack_CSC(rootsufR, (GrB_Index **)&IrootsufR, (void **)&VrootsufR, &Ibytes, &Valbytes, NULL, &nUfR, NULL, NULL)); // does it have space afterwards or should I pack again?
                GRB_TRY(GrB_Vector_clear(pathUpdate));
                GRB_TRY(GrB_Vector_build_UINT64(pathUpdate, VrootsufR, IrootsufR, nUfR, GrB_FIRST_UINT64)); // useful to handle duplicates
                GRB_TRY(GrB_Vector_assign(pathC, NULL, NULL, pathUpdate, GrB_ALL, nUfR, NULL));
                free(IrootsufR);
                free(VrootsufR);

                // get roots of row nodes in the current R frontier
                GRB_TRY(GrB_Vector_apply(rootsfR, NULL, NULL, getRootsOp, frontierR, NULL));

                GrB_Index *VmatesfR = (GrB_Index *)malloc(nrows * sizeof(GrB_Index));
                GrB_Index *VrootsfR = (GrB_Index *)malloc(nrows * sizeof(GrB_Index));
                GrB_Index nfR = 0, nRootsfR = 0;
                GrB_Index *dummy; // should allocate space for this erither way?
                GrB_Index n_dummy = 1, bytes_dummy = 0;
                // keep mates of the R frontier (ordered indices)
                GRB_TRY(GxB_Vector_unpack_CSC(mateR, (GrB_Index **)&dummy, (void **)&VmatesfR, &bytes_dummy, &Valbytes, NULL, &nfR, NULL, NULL));
                // keep roots of the R frontier (ordered indices)
                GRB_TRY(GxB_Vector_unpack_CSC(rootsfR, (GrB_Index **)&dummy, (void **)&VrootsfR, &bytes_dummy, &Ibytes, NULL, &nRootsfR, NULL, NULL));
                GRB_TRY(GrB_Vector_clear(rootfRIndexes));
                GRB_TRY(GrB_Vector_build_UINT64(rootfRIndexes, VrootsfR, VmatesfR, nRootsfR, GrB_FIRST_UINT64)); // rootfRIndexes(j) = i, where i is the col mate of the first row
                                                                                                                 // included in the current R frontier with a col root of j
                // keep only col roots that are not included in ufR
                GRB_TRY(GrB_Vector_assign(rootfRIndexes, pathUpdate, NULL, rootfRIndexes, GrB_ALL, ncols, GrB_DESC_RSC));
                free(VmatesfR);
                free(VrootsfR);

                GrB_Index *IrootfRIndexes = (GrB_Index *)malloc(ncols * sizeof(GrB_Index));
                GrB_Index *VrootfRIndexes = (GrB_Index *)malloc(ncols * sizeof(GrB_Index));
                GrB_Index nRootfRIndexes = 0;
                GRB_TRY(GxB_Vector_unpack_CSC(rootfRIndexes, (GrB_Index **)&IrootfRIndexes, (void **)&VrootfRIndexes, &Ibytes, &Valbytes, NULL, &nRootfRIndexes, NULL, NULL));
                GRB_TRY(GxB_Vector_pack_CSC(rootfRIndexes, (GrB_Index **)&VrootfRIndexes, (void **)&IrootfRIndexes, Valbytes, Ibytes, NULL, nRootfRIndexes, NULL, NULL)); // rootfRIndexes(i) = j,
                                                                                                                                                                          // where (i,j) = (parentC, rootC) of the new frontier C
                free(IrootfRIndexes);
                free(VrootfRIndexes);

                // build tuple of (parentC, rootC)
                GRB_TRY(GrB_Vector_apply_IndexOp_UDT(frontierC, NULL, NULL, buildfCTuplesOp, rootfRIndexes, &y, NULL));

                // /* debug
                GrB_Index C[ncols];
                vertex *V = malloc(ncols * sizeof(vertex));
                GrB_Vector_extractTuples_UDT(C, V, &ncols, frontierC);
                for (int k = 0; k < ncols; k++)
                {
                    printf("\nfc (%d) = (%ld, %ld)", (int)C[k], V[k].parentC, V[k].rootC);
                }
                // */
            }
            else
            {
                // apply op on frontier to set parents to mates
                GRB_TRY(GrB_Vector_apply(frontierR, NULL, NULL, setParentsMatesOp, mateR, NULL)); // fR(i) = (column mate of i, rootC)
                // invert fr
                GrB_Index *VmatesfR = (GrB_Index *)malloc(nrows * sizeof(GrB_Index));
                GrB_Index *VfR = (GrB_Index *)malloc(nrows * sizeof(GrB_Index));
                GrB_Index *dummy;
                GrB_Index bytes_dummy = 0, Vmatesbytes = 0, VfRBytes = 0, nfR = 0;
                GRB_TRY(GxB_Vector_unpack_CSC(mateR, (GrB_Index **)&dummy, (void **)&VmatesfR, &bytes_dummy, &Vmatesbytes, NULL, &nfR, NULL, NULL)); // mateR already only contains the rows of fR
                GRB_TRY(GxB_Vector_unpack_CSC(frontierR, (GrB_Index **)&dummy, (void **)&VfR, &bytes_dummy, &VfRBytes, NULL, &nfR, NULL, NULL));
                GRB_TRY(GxB_Vector_pack_CSC(frontierR, (GrB_Index **)&VmatesfR, (void **)&VfR, Vmatesbytes, VfRBytes, NULL, nfR, NULL, NULL));
                free(VmatesfR);
                free(VfR);
                // assign to fC
                GRB_TRY(GrB_Vector_assign(frontierC, NULL, NULL, frontierR, GrB_ALL, ncols, NULL));
            }

            GRB_TRY(GrB_Vector_nvals(&nvals, frontierC));

        } while (nvals);

        LG_ASSERT_MSG(1 == 0, GrB_INVALID_VALUE, "dummy");

        GRB_TRY(GrB_Vector_nvals(&nvals, pathC));
    } while (nvals); // only in the first and last iteration should this condition be false

    return (GrB_SUCCESS);
}