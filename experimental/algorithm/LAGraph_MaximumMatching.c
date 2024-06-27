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
        GrB_free(&AddMonoid);         \
        GrB_free(&MultOp);            \
        GrB_free(&semiring);          \
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
        GrB_free(&ur);                \
        GrB_free(&pathCopy);          \
        GrB_free(&currentMatesR);     \
    }

#define LG_FREE_ALL           \
    {                         \
        LG_FREE_WORK;         \
        GrB_free(&mateCcopy); \
    }

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

    GrB_Vector pathC = NULL;    // make this bitmap, if dense I would have to give all the entries and make the matrix 1-based
    GrB_Vector parentsR = NULL; // parents of row nodes that are reachable from paths of the initial column frontier
    GrB_Type Vertex = NULL;
    GrB_Vector frontierC = NULL;
    GrB_Vector frontierR = NULL;
    GrB_IndexUnaryOp initFrontierOp = NULL;
    GrB_Vector I = NULL; // dense vector of 1's
    GrB_BinaryOp MinParent = NULL;
    GrB_Monoid AddMonoid = NULL;
    GrB_BinaryOp MultOp = NULL;
    GrB_Semiring semiring = NULL;
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
    GrB_Vector ur = NULL;
    GrB_Vector pathCopy = NULL;
    GrB_Vector currentMatesR = NULL;

    GrB_Vector mateCcopy = *mateC;

    LG_CLEAR_MSG;

    LG_TRY(LAGraph_CheckGraph(G, msg));

    GrB_Matrix A = G->A;
    uint64_t ncols = 0;
    GRB_TRY(GrB_Matrix_ncols(&ncols, A));

    uint64_t nrows = 0;
    GRB_TRY(GrB_Matrix_nrows(&nrows, A));

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
    GRB_TRY(GrB_Monoid_new_UDT(&AddMonoid, MinParent, &infinityParent));

    GRB_TRY(GxB_BinaryOp_new(&MultOp, (void *)select2nd,
                             Vertex, GrB_BOOL, Vertex, "select2nd", SELECT_2ND_DEFN));

    GRB_TRY(GrB_Semiring_new(&semiring, AddMonoid, MultOp));

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

    GRB_TRY(GrB_Vector_new(&ur, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&pathCopy, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&currentMatesR, GrB_UINT64, nrows));

    uint64_t npath = 0;
    bool y = 0; // see if I can get rid of this

    do
    {
        GRB_TRY(GrB_Vector_clear(pathC));
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

        uint64_t nmatched = 0;
        GRB_TRY(GrB_Vector_nvals(&nmatched, mateCcopy));
        if (nmatched)
        {
            GrB_Index *J, *X; // unpack allocates space for these lists
            GrB_Index Jbytes = 0, Xbytes = 0;
            bool jumbledMateC = 0;
            GRB_TRY(GxB_Vector_unpack_CSC(mateCcopy, (GrB_Index **)&J, (void **)&X, &Jbytes, &Xbytes, NULL, &nmatched, &jumbledMateC, NULL)); // mateC and mateR do not have duplicates, so the order doesn't matter
            GRB_TRY(GrB_Vector_clear(mateR));                                                                                                 // clear mateR first as a prerequisite of the build method
            GRB_TRY(GrB_Vector_build_UINT64(mateR, X, J, nmatched, NULL));                                                                    // build does not take ownership of the lists J and X, but only copies them,
                                                                                                                                              // these lists will be given again to mateC
                                                                                                                                              // mateC has no duplicates in the values list, so mateR doesn't need to handle dups
            GRB_TRY(GxB_Vector_pack_CSC(mateCcopy, (GrB_Index **)&J, (void **)&X, Jbytes, Xbytes, NULL, nmatched, true, NULL));
        }

        /* debug
        GxB_Vector_fprint(mateR, "mateR", GxB_COMPLETE, stdout);
        */

        uint64_t nfC = 0;

        do
        {
            // perform one step of BFS from C nodes and keep only unvisited rows
            GRB_TRY(GrB_mxv(frontierR, parentsR, NULL, semiring, A, frontierC, GrB_DESC_RSC));
            // set parents of row frontier
            GRB_TRY(GrB_Vector_apply(parentsUpdate, NULL, NULL, getParentsOp, frontierR, NULL));                // previous values are erased
            GRB_TRY(GrB_Vector_assign(parentsR, NULL, GrB_SECOND_UINT64, parentsUpdate, GrB_ALL, nrows, NULL)); // update parents without deleting the ones not updated
                                                                                                                // when GrB_ALL is used, ni is the number of rows of the vector

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
                // get roots of unmatched row nodes in the R frontier
                GRB_TRY(GrB_Vector_apply(rootsufR, NULL, NULL, getRootsOp, ufrontierR, NULL));

                GrB_Index *IrootsufR, *VrootsufR;
                GrB_Index Ibytes = 0, Valbytes = 0;
                GRB_TRY(GxB_Vector_unpack_CSC(rootsufR, (GrB_Index **)&IrootsufR, (void **)&VrootsufR, &Ibytes, &Valbytes, NULL, &nUfR, NULL, NULL)); // sorted indices so we keep the min child
                GRB_TRY(GrB_Vector_clear(pathUpdate));
                GRB_TRY(GrB_Vector_build_UINT64(pathUpdate, VrootsufR, IrootsufR, nUfR, GrB_FIRST_UINT64));   // useful to handle duplicates
                GRB_TRY(GrB_Vector_assign(pathC, NULL, GrB_SECOND_UINT64, pathUpdate, GrB_ALL, ncols, NULL)); // update path without deleting the values not updated
                LG_TRY(LAGraph_Free((void **)&IrootsufR, msg));                                               // build copies the lists so they need to be freed
                LG_TRY(LAGraph_Free((void **)&VrootsufR, msg));

                // get roots of row nodes in the current R frontier
                GRB_TRY(GrB_Vector_apply(rootsfR, NULL, NULL, getRootsOp, frontierR, NULL));

                /* debug
                GxB_Vector_fprint(rootsfR, "rootsfR", GxB_COMPLETE, stdout);
                */

                if (nfR)
                {
                    GrB_Index *VmatesfR, *VrootsfR, *dummy;
                    GrB_Index nRootsfR = 0;
                    GrB_Index n_dummy = 1, bytes_dummy = 0;
                    // keep mates of the R frontier (ordered indices)
                    GRB_TRY(GxB_Vector_unpack_CSC(currentMatesR, (GrB_Index **)&dummy, (void **)&VmatesfR, &bytes_dummy, &Valbytes, NULL, &nfR, NULL, NULL));
                    // keep roots of the R frontier (ordered indices)
                    GRB_TRY(GxB_Vector_unpack_CSC(rootsfR, (GrB_Index **)&dummy, (void **)&VrootsfR, &bytes_dummy, &Ibytes, NULL, &nRootsfR, NULL, NULL));
                    GRB_TRY(GrB_Vector_clear(rootfRIndexes));
                    GRB_TRY(GrB_Vector_build_UINT64(rootfRIndexes, VrootsfR, VmatesfR, nRootsfR, GrB_FIRST_UINT64)); // rootfRIndexes(j) = i, where i is the col mate of the first row
                                                                                                                     // included in the current R frontier with a col root of j
                    // keep only col roots that are not included in ufR
                    GRB_TRY(GrB_Vector_assign(rootfRIndexes, pathUpdate, NULL, rootfRIndexes, GrB_ALL, ncols, GrB_DESC_RSC));
                    LG_TRY(LAGraph_Free((void **)&VmatesfR, msg));
                    LG_TRY(LAGraph_Free((void **)&VrootsfR, msg));

                    GrB_Index *IrootfRIndexes, *VrootfRIndexes;
                    GrB_Index nRootfRIndexes = 0;
                    bool jumbledRoots = 1;
                    GRB_TRY(GxB_Vector_unpack_CSC(rootfRIndexes, (GrB_Index **)&IrootfRIndexes, (void **)&VrootfRIndexes, &Ibytes, &Valbytes, NULL, &nRootfRIndexes, &jumbledRoots, NULL)); // no need to sort them
                    GRB_TRY(GxB_Vector_pack_CSC(rootfRIndexes, (GrB_Index **)&VrootfRIndexes, (void **)&IrootfRIndexes, Valbytes, Ibytes, NULL, nRootfRIndexes, true, NULL));               // rootfRIndexes(i) = j,
                                                                                                                                                                                            // where (i,j) = (parentC, rootC) of the new frontier C
                }

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
                // typecast mateR to ensure domain match with frontier R and apply op on frontier to set parents to mates
                GRB_TRY(GrB_Vector_apply(frontierR, NULL, setParentsMatesOp, vertexTypecastOp, currentMatesR, NULL)); // fR(i) = (column mate of i, rootC)  // add the structural mask
                // invert fr
                GrB_Index *VmatesfR, *VfR, *dummy;
                GrB_Index bytes_dummy = 0, Vmatesbytes = 0, VfRBytes = 0, nfR = 0;
                GRB_TRY(GxB_Vector_unpack_CSC(currentMatesR, (GrB_Index **)&dummy, (void **)&VmatesfR, &bytes_dummy, &Vmatesbytes, NULL, &nfR, NULL, NULL)); // currentMatesR already contains only the rows of fR
                GRB_TRY(GxB_Vector_unpack_CSC(frontierR, (GrB_Index **)&dummy, (void **)&VfR, &bytes_dummy, &VfRBytes, NULL, &nfR, NULL, NULL));
                GRB_TRY(GxB_Vector_pack_CSC(frontierR, (GrB_Index **)&VmatesfR, (void **)&VfR, Vmatesbytes, VfRBytes, NULL, nfR, true, NULL)); // the values are not ordered,
                                                                                                                                               // so the indices of the inverted fR are jumbled
                // assign to fC
                GRB_TRY(GrB_Vector_assign(frontierC, NULL, NULL, frontierR, GrB_ALL, ncols, GrB_DESC_RS));
            }

            GRB_TRY(GrB_Vector_nvals(&nfC, frontierC));

        } while (nfC);

        GRB_TRY(GrB_Vector_nvals(&npath, pathC));
        uint64_t npathCopy = npath;
        GrB_Index *Ipath, *Xpath;
        GrB_Index IpathBytes = 0, XpathBytes = 0;
        while (npath)
        {
            // invert pathC
            bool jumbledPathC = 1;
            GRB_TRY(GxB_Vector_unpack_CSC(pathC, (GrB_Index **)&Ipath, (void **)&Xpath, &IpathBytes, &XpathBytes, NULL, &npath, &jumbledPathC, NULL)); // pathC doesn't have dup values as it stems from an invertion
            GRB_TRY(GxB_Vector_pack_CSC(ur, (GrB_Index **)&Xpath, (void **)&Ipath, XpathBytes, IpathBytes, NULL, npath, true, NULL));                  // ur is already empty because in the previous iteration it was unpacked

            /* debug
            GxB_Vector_fprint(ur, "ur", GxB_COMPLETE, stdout);
            GxB_Vector_fprint(parentsR, "parentsR", GxB_COMPLETE, stdout);
            */

            // assign parents of rows to rows
            GRB_TRY(GrB_Vector_assign(ur, ur, NULL, parentsR, GrB_ALL, nrows, GrB_DESC_S)); // update the values of ur (descriptor needed to use mask's structure and not values)

            /* debug
            GxB_Vector_fprint(ur, "ur with updated parents", GxB_COMPLETE, stdout);
            */

            // invert ur
            bool jumbledUR = 1;
            GRB_TRY(GxB_Vector_unpack_CSC(ur, (GrB_Index **)&Ipath, (void **)&Xpath, &IpathBytes, &XpathBytes, NULL, &npath, &jumbledUR, NULL));
            GRB_TRY(GxB_Vector_pack_CSC(pathC, (GrB_Index **)&Xpath, (void **)&Ipath, XpathBytes, IpathBytes, NULL, npath, true, NULL));

            /* debug
            GxB_Vector_fprint(pathC, "pathC", GxB_COMPLETE, stdout);
            */

            // keep a copy of the previous row matches of the matched cols that will alter mates
            GRB_TRY(GrB_Vector_assign(pathCopy, pathC, NULL, mateCcopy, GrB_ALL, ncols, GrB_DESC_RS));

            /* debug
            GxB_Vector_fprint(pathCopy, "pathCopy", GxB_COMPLETE, stdout);
            */

            // update mateC
            GRB_TRY(GrB_Vector_assign(mateCcopy, NULL, GrB_SECOND_UINT64, pathC, GrB_ALL, ncols, NULL));
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

    return (GrB_SUCCESS);
}