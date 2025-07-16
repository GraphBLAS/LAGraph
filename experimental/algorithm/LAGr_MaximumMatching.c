//------------------------------------------------------------------------------
// LAGr_MaximumMatching: maximum matching between nodes of disjoint sets
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

//------------------------------------------------------------------------------

// This implmentation is based on the algorithm described in the following
// paper: "Distributed-Memory Algorithms for Maximum Cardinality Matching in
// Bipartite Graphs" by A. Azad and A. Buluç.
//
// In brief, this algorithm explores alternate paths by reversing an already
// traversed path to see if one of the already matched columns encountered in
// the path has at least one free children to be matched with instead. If so, it
// concedes its previous match to another previously matched parent or to the
// unmatched root of the path.

// More detailed explanation of the algorithm and its implementation can be
// found in the LAGraph "papers" folder, in file
// "MaximumMatching_Report.pdf".

// This is an Advanced method, G->A is required as input instead of G due to
// lack of a Bipartite graph kind in GraphBLAS.

//------------------------------------------------------------------------------

#include "LAGraphX.h"
#include "LG_internal.h"

//------------------------------------------------------------------------------
// the Vertex tuple: (parentC, rootC)
//------------------------------------------------------------------------------

typedef struct
{
    uint64_t parentC;
    uint64_t rootC;
} LG_MM_vertex;

// repeat the typedef as a string, to give to GraphBLAS
#define VERTEX_DEFN         \
"typedef struct         \n" \
"{                      \n" \
"    uint64_t parentC;  \n" \
"    uint64_t rootC;    \n" \
"} LG_MM_vertex;"

void LG_MM_initFrontier(LG_MM_vertex *z, const bool *x, uint64_t i, uint64_t j, const bool *y)
{
    z->parentC = i;
    z->rootC = i;
}

#define INIT_FRONTIER_DEFN  \
"void LG_MM_initFrontier(LG_MM_vertex *z, const bool *x, uint64_t i, uint64_t j, const bool *y) \n" \
"{                      \n" \
"    z->parentC = i;    \n" \
"    z->rootC = i;      \n" \
"}"

void LG_MM_minparent(LG_MM_vertex *z, const LG_MM_vertex *x, const LG_MM_vertex *y)
{
    *z = x->parentC < y->parentC ? *x : *y;
}

#define MIN_PARENT_DEFN                                             \
"void LG_MM_minparent(LG_MM_vertex *z, const LG_MM_vertex *x, const LG_MM_vertex *y)    \n" \
"{                                                              \n" \
"    *z = x->parentC < y->parentC ? *x : *y;                    \n" \
"}"

// TODO: revise GraphBLAS so we can tell it that the LG_MM_select2nd operator
// does not use the 'x' input.
void LG_MM_select2nd(LG_MM_vertex *z, const bool *x, const LG_MM_vertex *y)
{
    z->parentC = y->parentC;
    z->rootC = y->rootC;
}

#define SELECT_2ND_DEFN                                             \
"void LG_MM_select2nd(LG_MM_vertex *z, const bool *x, const LG_MM_vertex *y)      \n" \
"{                                                              \n" \
"    z->parentC = y->parentC;                                   \n" \
"    z->rootC = y->rootC;                                       \n" \
"}"

void LG_MM_select1st(LG_MM_vertex *z, const LG_MM_vertex *x, const bool *y)
{
    z->parentC = x->parentC;
    z->rootC = x->rootC;
}

#define SELECT_1ST_DEFN                                             \
"void LG_MM_select1st(LG_MM_vertex *z, const LG_MM_vertex *x, const bool *y)      \n" \
"{                                                              \n" \
"    z->parentC = x->parentC;                                   \n" \
"    z->rootC = x->rootC;                                       \n" \
"}"

void LG_MM_keepParents(uint64_t *z, const LG_MM_vertex *x) { *z = x->parentC; }

#define KEEP_PARENTS_DEFN                                           \
"void LG_MM_keepParents(uint64_t *z, const LG_MM_vertex *x) { *z = x->parentC; }\n"

void LG_MM_keepRoots(uint64_t *z, const LG_MM_vertex *x) { *z = x->rootC; }

#define KEEP_ROOTS_DEFN                                                        \
"void LG_MM_keepRoots(uint64_t *z, const LG_MM_vertex *x) { *z = x->rootC; }\n"

void LG_MM_buildfCTuples(LG_MM_vertex *z, const uint64_t *x, uint64_t i, uint64_t j,
                   const bool *y)
{
    z->parentC = i;
    z->rootC = *x;
}

#define BUILT_FC_TUPLES_DEFN                                        \
"void LG_MM_buildfCTuples(LG_MM_vertex *z, const uint64_t *x, uint64_t i, uint64_t j,   \n" \
"                   const bool *y)                              \n" \
"{                                                              \n" \
"    z->parentC = i;                                            \n" \
"    z->rootC = *x;                                             \n" \
"}"

void LG_MM_vertexTypecast(LG_MM_vertex *z, const uint64_t *x)
{
    z->parentC = *x;
    z->rootC = *x;
}

#define VERTEX_TYPECAST_DEFN                                        \
"void LG_MM_vertexTypecast(LG_MM_vertex *z, const uint64_t *x)              \n" \
"{                                                              \n" \
"    z->parentC = *x;                                           \n" \
"    z->rootC = *x;                                             \n" \
"}"

void LG_MM_setParentsMates(LG_MM_vertex *z, const LG_MM_vertex *x, const LG_MM_vertex *y)
{
    z->parentC = y->parentC;
    z->rootC = x->rootC;
}

#define SET_PARENTS_MATES_DEFN                                          \
"void LG_MM_setParentsMates(LG_MM_vertex *z, const LG_MM_vertex *x, const LG_MM_vertex *y)  \n" \
"{                                                                  \n" \
"    z->parentC = y->parentC;                                       \n" \
"    z->rootC = x->rootC;                                           \n" \
"}"

//------------------------------------------------------------------------------
// invert
//------------------------------------------------------------------------------

// This function "inverts" an input vector by swapping its row indices
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
// range 0 to n-1 where n = length(out). The values in the input vector
// may be duplicated and this argument of the function must be set accordingly.
// Both the in vector and out vector must have the same type (GrB_UINT64).  The
// lengths of the two vectors need not be the same, so long as the indices
// remain in range.  Results are undefined if these conditions do not hold.
//
// The in and out vectors may be aliased.  If not aliased, the input vector is
// cleared of all entries on output.  If in and out are aliased, then the
// inversion is performed in-place.
//
// In SuiteSparse:GraphBLAS, this method takes O(1) time if the in vector is in
// CSC (sparse, by column) format.  Otherwise it can take O(e) time if e =
// nvals(in), because the unpack below will convert the in vector to CSC and
// then unpack it.

#undef LG_FREE_ALL
#define LG_FREE_ALL                                                            \
    {                                                                          \
        LAGraph_Free((void *)&I, NULL);                                        \
        LAGraph_Free((void *)&X1, NULL);                                       \
        LAGraph_Free((void *)&X2, NULL);                                       \
    }

static inline GrB_Info invert_nondestructive(
    GrB_Vector out, // input/output.  On input, only the size and type are
                    // kept; any entries in the 'out' vector cleared.  It is
                    // then replaced with the inversion of the input vector.
    GrB_Vector in,  // input vector, not modified.  There must be no duplicate
                    // values in the input vector.
    char *msg)
{
    bool jumbled = 1;
    GrB_Index *I = NULL;
    GrB_Index *X1 = NULL;
    GrB_Index *X2 = NULL; // not used
    GrB_Index IBytes = 0, XBytes = 0;
    uint64_t nvals = 0;

    // the input and output vectors cannot be the same vector
    ASSERT(in != out);

    // All input/output vectors must be of type GrB_UINT64.
#if LAGRAPH_SUITESPARSE
    GRB_TRY(
        GxB_Vector_unpack_CSC(in, (GrB_Index **)&I, (void **)&X1, &IBytes,
                              &XBytes, NULL, &nvals, &jumbled,
                              NULL)); // the output and input should have no
                                      // duplicates, so the order doesn't matter
#else
    GRB_TRY(GrB_Vector_nvals(&nvals, in));
    LG_TRY(LAGraph_Malloc((void **)&I, nvals, sizeof(GrB_Index), msg));
    LG_TRY(LAGraph_Malloc((void **)&X1, nvals, sizeof(GrB_Index), msg));
    GRB_TRY(GrB_Vector_wait(in, GrB_MATERIALIZE));
    GRB_TRY(GrB_Vector_extractTuples_UINT64(
        I, X1, &nvals, in)); // the output and input should have no
                             // duplicates, so the order doesn't matter
#endif

    GRB_TRY(GrB_Vector_clear(out)); // clear the output first as a prerequisite
                                    // of the build method
    GRB_TRY(GrB_Vector_build_UINT64(
        out, X1, I, nvals,
        NULL)); // build does not take ownership of the lists I and X,
                // but only copies them, these lists will be given
                // again to the input
                // the input should have no duplicates in the
                // values list, so dups are not handled
#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_Vector_pack_CSC(in, (GrB_Index **)&I, (void **)&X1, IBytes,
                                XBytes, false, nvals, jumbled, NULL));
#endif
    return (GrB_SUCCESS);
}

static inline GrB_Info
invert(GrB_Vector out,  // input/output.  Same as invert_nondescructive above.
       GrB_Vector in,   // input vector, empty on output (unless in == out)
       const bool dups, // flag that indicates if duplicates exist in the input
                        // vector's values
       char *msg)
{
    // The input and output vectors can be the same vector
    // that is, in == out is OK.
    // All input/output vectors must be of type GrB_UINT64.

    // the output vector will normally be returned in a jumbled state
    bool jumbled = dups ? 0 : 1; // if there are duplicates, we want the indices
                                 // to be ordered so we can keep the min child
    GrB_Index *I = NULL;         // unpack allocates space for these lists
    GrB_Index *X1 = NULL;
    GrB_Index *X2 = NULL; // not used
    GrB_Index IBytes = 0, XBytes = 0;
    uint64_t nvals = 0;

#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_Vector_unpack_CSC(in, (GrB_Index **)&I, (void **)&X1, &IBytes,
                                  &XBytes, NULL, &nvals, &jumbled, NULL));
#else
    // vanilla case using extractTuples and build:
    // allocate I and X for GrB_extractTuples
    GRB_TRY(GrB_Vector_nvals(&nvals, in));
    LG_TRY(LAGraph_Malloc((void **)&I, nvals, sizeof(GrB_Index), msg));
    LG_TRY(LAGraph_Malloc((void **)&X1, nvals, sizeof(GrB_Index), msg));
    GRB_TRY(GrB_Vector_wait(in, GrB_MATERIALIZE));
    GRB_TRY(GrB_Vector_extractTuples_UINT64(
        I, X1, &nvals, in)); // the output and input should have no
                             // duplicates, so the order doesn't matter
    GRB_TRY(GrB_Vector_clear(in));
#endif
    if (!dups)
    {
#if LAGRAPH_SUITESPARSE
        GRB_TRY(GxB_Vector_pack_CSC(out, (GrB_Index **)&X1, (void **)&I, XBytes,
                                    IBytes, false, nvals, true, NULL));
#else
        GRB_TRY(GrB_Vector_clear(out));
        // GrB_MIN_UINT64 is used instead of first because first is an extension
        GRB_TRY(GrB_Vector_build_UINT64(out, X1, I, nvals, GrB_MIN_UINT64));
        // build copies the lists so they need to be freed in LG_FREE_ALL
        LG_FREE_ALL;
#endif
    }
    else
    {
        GRB_TRY(GrB_Vector_clear(out));
#if LAGRAPH_SUITESPARSE
        GRB_TRY(GrB_Vector_build_UINT64(out, X1, I, nvals, GrB_FIRST_UINT64));
#else
        GRB_TRY(GrB_Vector_build_UINT64(out, X1, I, nvals, GrB_MIN_UINT64));
#endif
        // build copies the lists so they need to be freed in LG_FREE_ALL
        LG_FREE_ALL;
    }
    return (GrB_SUCCESS);
}

static inline GrB_Info
invert_2(GrB_Vector out,  // input/output
         GrB_Vector in1,  // input vector, empty on output (unless in1 == out)
         GrB_Vector in2,  // input vector, empty on output (unless in2 == out)
         const bool dups, // flag that indicates if duplicates exist in the
                          // input vector's values
         const bool udt,  // type of the first input, in1, and out vectors
         char *msg)
{
    // The input vectors cannot be aliased.  However in1==out or in2==out is
    // OK.  The two input vectors must have the same # of entries.
    // All input/output vectors must be of type GrB_UINT64.
    ASSERT(in1 != in2);

    GrB_Index *I = NULL;
    GrB_Index *X1 = NULL;
    GrB_Index *X2 = NULL;
    GrB_Index IBytes = 0, X1Bytes = 0, X2Bytes = 0;
    uint64_t nvals1 = 0, nvals2 = 0;

#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_Vector_unpack_CSC(in1, (GrB_Index **)&I, (void **)&X1, &IBytes,
                                  &X1Bytes, NULL, &nvals1, NULL, NULL));
    LAGraph_Free((void *)&I, NULL);
    GRB_TRY(GxB_Vector_unpack_CSC(in2, (GrB_Index **)&I, (void **)&X2, &IBytes,
                                  &X2Bytes, NULL, &nvals2, NULL, NULL));
    ASSERT(nvals1 == nvals2);
    if (!dups)
    {
        LAGraph_Free((void *)&I, NULL);
        GRB_TRY(GxB_Vector_pack_CSC(out, (GrB_Index **)&X2, (void **)&X1,
                                    X2Bytes, X1Bytes, false, nvals2, true,
                                    NULL)); // the values are not ordered,
                                            // so the indices are jumbled
    }
    else
    {
        GRB_TRY(GrB_Vector_clear(out));
        GRB_TRY(GrB_Vector_build_UINT64(out, X2, X1, nvals2, GrB_FIRST_UINT64));
        LG_FREE_ALL;
    }
#else
    // vanilla case using extractTuples and build:
    if (dups)
    {
        // invert in2 to eliminate dups and get the final indices of out (we
        // cannot invert into out if out is of type udt, because we would have a
        // domain mismatch with in2)
        GrB_Vector in_rev = NULL;
        GrB_Index out_size = 0;
        GRB_TRY(GrB_Vector_size(&out_size, out));
        GRB_TRY(GrB_Vector_new(&in_rev, GrB_UINT64, out_size));
        LG_TRY(invert(in_rev, in2, true, msg));
        GRB_TRY(GrB_Vector_nvals(
            &nvals2,
            in_rev)); // out may have less entries than in2 and, thus, in1
        LG_TRY(LAGraph_Malloc((void **)&I, nvals2, sizeof(GrB_Index), msg));
        LG_TRY(LAGraph_Malloc((void **)&X2, nvals2, sizeof(GrB_Index), msg));
        // get indices of out
        GRB_TRY(GrB_Vector_extractTuples_UINT64(
            X2, I, &nvals2,
            in_rev)); // X2 are the final values of in2 and will be the indices
                      // of out. I are the indices of in2 (which are a subset of
                      // the indices of in1)
        GRB_TRY(GrB_Vector_free(&in_rev));
    }
    else
    {
        GRB_TRY(GrB_Vector_nvals(&nvals2, in2));
        LG_TRY(LAGraph_Malloc((void **)&I, nvals2, sizeof(GrB_Index), msg));
        LG_TRY(LAGraph_Malloc((void **)&X2, nvals2, sizeof(GrB_Index), msg));
        GRB_TRY(GrB_Vector_wait(in2, GrB_MATERIALIZE));
        GRB_TRY(GrB_Vector_extractTuples_UINT64(
            I, X2, &nvals2, in2)); // X2 will be the indices of out and I are
                                   // the indices of in2 (which are the indices
                                   // of in1 as well)
    }

    GRB_TRY(GrB_Vector_clear(out));
    GRB_TRY(GrB_Vector_wait(in1, GrB_MATERIALIZE));
    if (udt)
    {
        LG_MM_vertex *X1_V = NULL;
        LG_TRY(LAGraph_Malloc((void **)&X1_V, nvals2, sizeof(LG_MM_vertex), msg));
        for (uint64_t i = 0; i < nvals2; i++)
        {
            GRB_TRY(GrB_Vector_extractElement_UDT(
                (void *)&X1_V[i], in1,
                (GrB_Index)I[i])); // value I[i] corresponds to pos X2[i] in
                                   // out due to the tuple extraction and X[i]
                                   // will replace it
        }
        GRB_TRY(GrB_Vector_clear(in1));
        GRB_TRY(GrB_Vector_build_UDT(out, X2, (void *)X1_V, nvals2,
                                     NULL)); // no dups are left
        LAGRAPH_TRY(LAGraph_Free((void **)&X1_V, msg));
    }
    else
    {
        LG_TRY(LAGraph_Malloc((void **)&X1, nvals2, sizeof(GrB_UINT64), msg));
        for (GrB_Index i = 0; i < nvals2; i++)
        {
            GRB_TRY(GrB_Vector_extractElement_UINT64(&X1[i], in1, I[i]));
        }
        GRB_TRY(GrB_Vector_clear(in1));
        GRB_TRY(GrB_Vector_build_UINT64(out, X2, X1, nvals2, NULL));
    }
    LG_FREE_ALL;
#endif
    return (GrB_SUCCESS);
}

//------------------------------------------------------------------------------
// LAGr_MaximumMatching
//------------------------------------------------------------------------------

#undef LG_FREE_WORK
#define LG_FREE_WORK                                                           \
    {                                                                          \
        GrB_free(&pathC);                                                      \
        GrB_free(&parentsR);                                                   \
        GrB_free(&Vertex);                                                     \
        GrB_free(&frontierC);                                                  \
        GrB_free(&frontierR);                                                  \
        GrB_free(&initFrontierOp);                                             \
        GrB_free(&I);                                                          \
        GrB_free(&MinParent);                                                  \
        GrB_free(&MinParent_Monoid);                                           \
        GrB_free(&Select2ndOp);                                                \
        GrB_free(&Select1stOp);                                                \
        GrB_free(&MinParent_2nd_Semiring);                                     \
        GrB_free(&MinParent_1st_Semiring);                                     \
        GrB_free(&getParentsOp);                                               \
        GrB_free(&getRootsOp);                                                 \
        GrB_free(&parentsUpdate);                                              \
        GrB_free(&ufrontierR);                                                 \
        GrB_free(&mateC);                                                      \
        GrB_free(&mateR);                                                      \
        GrB_free(&rootsufR);                                                   \
        GrB_free(&pathUpdate);                                                 \
        GrB_free(&rootufRIndexes);                                             \
        GrB_free(&rootsfR);                                                    \
        GrB_free(&rootfRIndexes);                                              \
        GrB_free(&buildfCTuplesOp);                                            \
        GrB_free(&vertexTypecastOp);                                           \
        GrB_free(&setParentsMatesOp);                                          \
        GrB_free(&vr);                                                         \
        GrB_free(&pathCopy);                                                   \
        GrB_free(&currentMatesR);                                              \
    }

#undef LG_FREE_ALL
#define LG_FREE_ALL                                                            \
    {                                                                          \
        LG_FREE_WORK;                                                          \
        GrB_free(&mateC);                                                      \
    }

int LAGr_MaximumMatching(
    // outputs:
    GrB_Vector
        *mateC_handle, // mateC(j) = i : Column j of the C subset is matched
                       // to row i of the R subset (ignored on input)
    GrB_Vector *mateR_handle, // mateR(i) = j : Row i of the R subset is matched
                              // to column j of the C subset (ignored on input)
    // inputs:
    GrB_Matrix A,  // input adjacency matrix, TODO: this should be a LAGraph
                   // of a BIPARTITE kind
    GrB_Matrix AT, // transpose of the input adjacency matrix, necessary to
                   // perform push-pull optimization
                   // if NULL, the push-pull optimization is not performed
    GrB_Vector mate_init, // input only, not modified, ignored if NULL
    bool col_init, // flag to indicate if the initial matching is provided from
                   // the columns' or from the rows' perspective, ignored if
                   // mate_init is NULL
    char *msg)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Vector pathC = NULL;    // make this bitmap/sparse,
                                // if dense I would have to give
                                //  all the entries and make the matrix 1-based
    GrB_Vector parentsR = NULL; // parents of row nodes that are reachable
                                // from paths of the initial column frontier
    GrB_Type Vertex = NULL;
    GrB_Vector frontierC = NULL;
    GrB_Vector frontierR = NULL;
    GrB_IndexUnaryOp initFrontierOp = NULL;
    GrB_Vector I = NULL; // dense vector of 1's
    GrB_BinaryOp MinParent = NULL;
    GrB_Monoid MinParent_Monoid = NULL;
    GrB_BinaryOp Select2ndOp = NULL;
    GrB_BinaryOp Select1stOp = NULL;
    GrB_Semiring MinParent_2nd_Semiring = NULL;
    GrB_Semiring MinParent_1st_Semiring = NULL;
    GrB_UnaryOp getParentsOp = NULL;
    GrB_UnaryOp getRootsOp = NULL;
    GrB_Vector parentsUpdate = NULL; // unmatched rows of R frontier
    GrB_Vector ufrontierR = NULL;    // unmatched rows of R frontier
    GrB_Vector mateR = NULL;         // mateR(i) = j : Row i of the R subset is
                                     // matched to column j of the C subset
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

    GrB_Vector mateC = NULL;

    LG_CLEAR_MSG;

    LG_ASSERT_MSG(mateC_handle != NULL || mateR_handle != NULL,
                  GrB_NULL_POINTER, "At least one output must be provided\n");
    LG_ASSERT_MSG(A != NULL || AT != NULL, GrB_NULL_POINTER,
                  "A matrix is NULL\n");

    if (mateC_handle != NULL)
    {
        (*mateC_handle) = NULL;
    }
    if (mateR_handle != NULL)
    {
        (*mateR_handle) = NULL;
    }

    bool do_pushpull = (AT != NULL) && (A != NULL);

    uint64_t ncols = 0;
    uint64_t nrows = 0;

    if (A != NULL)
    {
        GRB_TRY(GrB_Matrix_ncols(&ncols, A));
        GRB_TRY(GrB_Matrix_nrows(&nrows, A));
    }
    else
    {
        GRB_TRY(GrB_Matrix_nrows(&ncols, AT));
        GRB_TRY(GrB_Matrix_ncols(&nrows, AT));
    }

    GRB_TRY(GrB_Vector_new(&pathC, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&parentsR, GrB_UINT64, nrows));

#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_Type_new(&Vertex, sizeof(LG_MM_vertex), "LG_MM_vertex", VERTEX_DEFN));
#else
    GRB_TRY(GrB_Type_new(&Vertex, sizeof(LG_MM_vertex)));
#endif

    GRB_TRY(GrB_Vector_new(&frontierC, Vertex, ncols));

    GRB_TRY(GrB_Vector_new(&frontierR, Vertex, nrows));

#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_IndexUnaryOp_new(&initFrontierOp, (void *)LG_MM_initFrontier, Vertex,
                                 GrB_BOOL, GrB_BOOL, "LG_MM_initFrontier",
                                 INIT_FRONTIER_DEFN));
#else
    GRB_TRY(GrB_IndexUnaryOp_new(&initFrontierOp, (void *)LG_MM_initFrontier, Vertex,
                                 GrB_BOOL, GrB_BOOL));
#endif

    GRB_TRY(GrB_Vector_new(&I, GrB_BOOL, ncols));
    GRB_TRY(GrB_Vector_assign_BOOL(I, NULL, NULL, 1, GrB_ALL, ncols, NULL));

#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_BinaryOp_new(&MinParent, (void *)LG_MM_minparent, Vertex, Vertex,
                             Vertex, "LG_MM_minparent", MIN_PARENT_DEFN));
#else
    GRB_TRY(GrB_BinaryOp_new(&MinParent, (void *)LG_MM_minparent, Vertex, Vertex,
                             Vertex));
#endif
    LG_MM_vertex infinityParent = {GrB_INDEX_MAX + 1, 0};
    GRB_TRY(GrB_Monoid_new_UDT(&MinParent_Monoid, MinParent, &infinityParent));

#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_BinaryOp_new(&Select2ndOp, (void *)LG_MM_select2nd, Vertex, GrB_BOOL,
                             Vertex, "LG_MM_select2nd", SELECT_2ND_DEFN));
#else
    GRB_TRY(GrB_BinaryOp_new(&Select2ndOp, (void *)LG_MM_select2nd, Vertex, GrB_BOOL,
                             Vertex));
#endif

    GRB_TRY(GrB_Semiring_new(&MinParent_2nd_Semiring, MinParent_Monoid,
                             Select2ndOp));

#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_BinaryOp_new(&Select1stOp, (void *)LG_MM_select1st, Vertex, Vertex,
                             GrB_BOOL, "LG_MM_select1st", SELECT_1ST_DEFN));
#else
    GRB_TRY(GrB_BinaryOp_new(&Select1stOp, (void *)LG_MM_select1st, Vertex, Vertex,
                             GrB_BOOL));
#endif

    GRB_TRY(GrB_Semiring_new(&MinParent_1st_Semiring, MinParent_Monoid,
                             Select1stOp));

#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_UnaryOp_new(&getParentsOp, (void *)LG_MM_keepParents, GrB_UINT64,
                            Vertex, "LG_MM_keepParents", KEEP_PARENTS_DEFN));
#else
    GRB_TRY(GrB_UnaryOp_new(&getParentsOp, (void *)LG_MM_keepParents, GrB_UINT64,
                            Vertex));
#endif

#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_UnaryOp_new(&getRootsOp, (void *)LG_MM_keepRoots, GrB_UINT64, Vertex,
                            "LG_MM_keepRoots", KEEP_ROOTS_DEFN));
#else
    GRB_TRY(
        GrB_UnaryOp_new(&getRootsOp, (void *)LG_MM_keepRoots, GrB_UINT64, Vertex));
#endif

    GRB_TRY(GrB_Vector_new(&parentsUpdate, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&ufrontierR, Vertex, nrows));

    GRB_TRY(GrB_Vector_new(&rootsufR, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&pathUpdate, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&rootufRIndexes, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&rootsfR, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&rootfRIndexes, GrB_UINT64, ncols));

#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_IndexUnaryOp_new(&buildfCTuplesOp, (void *)LG_MM_buildfCTuples,
                                 Vertex, GrB_UINT64, GrB_BOOL, "LG_MM_buildfCTuples",
                                 BUILT_FC_TUPLES_DEFN));
#else
    GRB_TRY(GrB_IndexUnaryOp_new(&buildfCTuplesOp, (void *)LG_MM_buildfCTuples,
                                 Vertex, GrB_UINT64, GrB_BOOL));
#endif

#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_UnaryOp_new(&vertexTypecastOp, (void *)LG_MM_vertexTypecast, Vertex,
                            GrB_UINT64, "LG_MM_vertexTypecast",
                            VERTEX_TYPECAST_DEFN));
#else
    GRB_TRY(GrB_UnaryOp_new(&vertexTypecastOp, (void *)LG_MM_vertexTypecast, Vertex,
                            GrB_UINT64));
#endif

#if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_BinaryOp_new(&setParentsMatesOp, (void *)LG_MM_setParentsMates,
                             Vertex, Vertex, Vertex, "LG_MM_setParentsMates",
                             SET_PARENTS_MATES_DEFN));
#else
    GRB_TRY(GrB_BinaryOp_new(&setParentsMatesOp, (void *)LG_MM_setParentsMates,
                             Vertex, Vertex, Vertex));
#endif

    GRB_TRY(GrB_Vector_new(&vr, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&pathCopy, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&currentMatesR, GrB_UINT64, nrows));

    uint64_t npath = 0;
    bool y = 0;

    GRB_TRY(GrB_Vector_new(&mateC, GrB_UINT64, ncols));
    GRB_TRY(GrB_Vector_new(&mateR, GrB_UINT64, nrows));

    if (mate_init != NULL)
    {
        uint64_t nmatched = 0;
        GRB_TRY(GrB_Vector_nvals(&nmatched, mate_init));
        if (nmatched)
        {
            if (col_init)
            {
                // mateC = (uint64_t) mate_init
                GRB_TRY(GrB_assign(mateC, NULL, NULL, mate_init, GrB_ALL, ncols,
                                   NULL));
                // mateR = invert (mateC), but do not clear the input
                LAGRAPH_TRY(invert_nondestructive(mateR, mateC, msg));
            }
            else
            {
                // mateR = (uint64_t) mate_init
                GRB_TRY(GrB_assign(mateR, NULL, NULL, mate_init, GrB_ALL, nrows,
                                   NULL));
                // mateC = invert (mateR), but do not clear the input
                LAGRAPH_TRY(invert_nondestructive(mateC, mateR, msg));
            }
        }
    }

    do
    {
        GRB_TRY(GrB_Vector_clear(parentsR));
        // for every col j not matched, assign f(j) = VERTEX(j,j)
        GRB_TRY(GrB_Vector_apply_IndexOp_BOOL(
            frontierC, mateC, NULL, initFrontierOp, I, true, GrB_DESC_RSC));

        uint64_t nfC = 0;

        do
        {
            //----------------------------------------------------------------------
            // STEPS 1,2: Explore neighbors of column frontier (one step of
            // BFS), keeping only unvisited rows in the frontierR
            //----------------------------------------------------------------------
            if (do_pushpull)
            {
                int32_t kind;
                LAGRAPH_TRY(LG_GET_FORMAT_HINT(frontierC, &kind));
                if (kind == LG_BITMAP || kind == LG_FULL)
                {
                    // the frontierC vector is bitmap or full
                    // pull (vector's values are pulled into A)
                    GRB_TRY(GrB_mxv(frontierR, parentsR, NULL,
                                    MinParent_2nd_Semiring, A, frontierC,
                                    GrB_DESC_RSC));
                }
                else
                {
                    // the frontierC vector is sparse or hypersparse
                    // push (vector's values are pushed to A)
                    GRB_TRY(GrB_vxm(frontierR, NULL, NULL,
                                    MinParent_1st_Semiring, frontierC, AT,
                                    GrB_DESC_R));
                    GRB_TRY(GrB_Vector_assign(frontierR, parentsR, NULL,
                                              frontierR, GrB_ALL, nrows,
                                              GrB_DESC_RSC));
                }
            }
            else
            {
                if (A != NULL)
                {
                    // Only the pull method can be used if AT is not given
                    GRB_TRY(GrB_mxv(frontierR, parentsR, NULL,
                                    MinParent_2nd_Semiring, A, frontierC,
                                    GrB_DESC_RSC));
                }
                else
                {
                    // Only the push method can be used if A is not given
                    GRB_TRY(GrB_vxm(frontierR, NULL, NULL,
                                    MinParent_1st_Semiring, frontierC, AT,
                                    GrB_DESC_R));
                    GRB_TRY(GrB_Vector_assign(frontierR, parentsR, NULL,
                                              frontierR, GrB_ALL, nrows,
                                              GrB_DESC_RSC));
                }
            }

            //----------------------------------------------------------------------
            // STEPS 3,4: Select univisited, matched and unmatched row vertices
            //----------------------------------------------------------------------
            // set parents of row frontier
            GRB_TRY(GrB_Vector_apply(
                parentsR, frontierR, NULL, getParentsOp, frontierR,
                GrB_DESC_S)); // use input as mask to only update or insert
                              // parents without deleting the ones not
                              // updated

            // select unmatched rows of the R frontier
            GRB_TRY(GrB_Vector_assign(ufrontierR, mateR, NULL, frontierR,
                                      GrB_ALL, nrows, GrB_DESC_RSC));
            // select matched rows of the R frontier
            GRB_TRY(GrB_Vector_assign(frontierR, mateR, NULL, frontierR,
                                      GrB_ALL, nrows, GrB_DESC_RS));

            // keep only mates of rows in frontierR
            GRB_TRY(GrB_Vector_assign(currentMatesR, frontierR, NULL, mateR,
                                      GrB_ALL, nrows, GrB_DESC_RS));

            uint64_t nUfR = 0, nfR = 0;
            GRB_TRY(GrB_Vector_nvals(&nUfR, ufrontierR));
            GRB_TRY(GrB_Vector_nvals(&nfR, frontierR));

            if (nUfR)
            {
                //----------------------------------------------------------------------
                // STEP 5: Store endpoints of newly discovered augmenting paths
                //----------------------------------------------------------------------
                // get roots of unmatched row nodes in the R frontier
                GRB_TRY(GrB_Vector_apply(rootsufR, NULL, NULL, getRootsOp,
                                         ufrontierR, NULL));

                // pathUpdate = invert (rootsufR), but need to handle
                // duplicates
                LAGRAPH_TRY(invert(pathUpdate, rootsufR, /*dups=*/true, msg));

                GRB_TRY(GrB_Vector_assign(
                    pathC, pathUpdate, NULL, pathUpdate, GrB_ALL, ncols,
                    GrB_DESC_S)); // update path without deleting the values
                                  // not updated when GrB_ALL is used, ni is
                                  // the number of rows of the vector

                //----------------------------------------------------------------------
                // STEP 6: Prune vertices in trees yielding augmenting paths
                //----------------------------------------------------------------------
                GRB_TRY(GrB_Vector_clear(rootfRIndexes));

                if (nfR)
                {
                    // get roots of row nodes in the current R frontier
                    GRB_TRY(GrB_Vector_apply(rootsfR, NULL, NULL, getRootsOp,
                                             frontierR, NULL));

                    // keep mates and roots of the R frontier (ordered indices)
                    LAGRAPH_TRY(
                        invert_2(rootfRIndexes, currentMatesR, rootsfR,
                                 /*dups=*/true, /*udt=*/false,
                                 msg)); // rootsfRIndexes(j) = i, where i
                                        // is the col mate of the first
                                        // row included in the current R
                                        // frontier with a col root of j

                    // keep only col roots that are not included in ufR
                    GRB_TRY(GrB_Vector_assign(rootfRIndexes, pathUpdate, NULL,
                                              rootfRIndexes, GrB_ALL, ncols,
                                              GrB_DESC_RSC));

                    //----------------------------------------------------------------------
                    // STEP 7a (ufrontierR not empty): Move values in the
                    // correct positions for the C frontier
                    //----------------------------------------------------------------------
                    // rootfRIndexes = invert (rootfRIndexes), so that:
                    // if LAGRAPH_SUITESPARSE: rootfRIndexes(i) = j,
                    // where (i,j) = (parentC, rootC) of the new frontier C
                    // else: rootfRIndexes(i) = j, where i is the matched child
                    // of the current frontier R that stems from root path of j
                    LAGRAPH_TRY(invert(rootfRIndexes, rootfRIndexes,
                                       /*dups=*/false, msg));
                }

                //----------------------------------------------------------------------
                // STEP 7b (ufrontierR not empty): Build tuple of (parentC,
                // rootC)
                //----------------------------------------------------------------------
                GRB_TRY(GrB_Vector_clear(frontierC));
                GRB_TRY(GrB_Vector_apply_IndexOp_BOOL(frontierC, NULL, NULL,
                                                     buildfCTuplesOp,
                                                     rootfRIndexes, true, NULL));
            }
            else
            {
                //----------------------------------------------------------------------
                // STEP 7a (ufrontierR is empty): Set parents of the R frontier
                // to their mates
                //----------------------------------------------------------------------
                // typecast mateR to ensure domain match with frontier R and
                // apply op on frontier to set parents to mates
                GRB_TRY(
                    GrB_Vector_apply(frontierR, NULL, setParentsMatesOp,
                                     vertexTypecastOp, currentMatesR,
                                     NULL)); // fR(i) = (column mate of i,
                                             // rootC) add the structural mask

                //----------------------------------------------------------------------
                // STEP 7b (ufrontierR is empty): Move values in the correct
                // positions for the C frontier
                //----------------------------------------------------------------------
                // invert fr and assign to fC
                // (currentMatesR already contains only the rows of fR)
                LAGRAPH_TRY(invert_2(frontierC, frontierR, currentMatesR,
                                     /*dups=*/false, /*udt=*/true, msg));
            }

            GRB_TRY(GrB_Vector_nvals(&nfC, frontierC));

        } while (nfC);

        //----------------------------------------------------------------------
        // STEP 8: Augment matching by all augmenting paths discovered in
        // this phase
        //----------------------------------------------------------------------
        GRB_TRY(GrB_Vector_nvals(&npath, pathC));
        uint64_t npathCopy = npath;
        while (npath)
        {
            // vr = invert (pathC), leaving pathC empty
            // pathC doesn't have dup values as it stems from an invertion
            LAGRAPH_TRY(invert(vr, pathC, /*dups=*/false, msg));

            // assign parents of rows to rows
            GRB_TRY(GrB_Vector_assign(
                vr, vr, NULL, parentsR, GrB_ALL, nrows,
                GrB_DESC_S)); // update the values of vr (descriptor needed
                              // to use mask's structure and not values)

            // update mateR:  mateR<vr> = vr
            GRB_TRY(GrB_Vector_assign(mateR, vr, NULL, vr, GrB_ALL, nrows,
                                      GrB_DESC_S));

            // pathC = invert (vr), leaving vr empty
            LAGRAPH_TRY(invert(pathC, vr, /*dups=*/false, msg));

            // keep a copy of the previous row matches of the matched cols
            // that will alter mates
            GRB_TRY(GrB_Vector_assign(pathCopy, pathC, NULL, mateC, GrB_ALL,
                                      ncols, GrB_DESC_RS));

            // update mateC
            GRB_TRY(GrB_Vector_assign(mateC, pathC, NULL, pathC, GrB_ALL, ncols,
                                      GrB_DESC_S));

            // At this point, mateR and mateC are in sync.  That is, they
            // are inversions of each other (mateR == invert(mateC) and
            // mateC == invert (mateR) both hold).

            // swap path and pathCopy
            GrB_Vector tmp = pathC;
            pathC = pathCopy;
            pathCopy = tmp;

            GRB_TRY(GrB_Vector_nvals(&npath, pathC));
        }

        npath = npathCopy;
    } while (npath); // only in the first and last iteration should this
                     // condition be false

    // return result
    if (mateC_handle != NULL)
    {
        (*mateC_handle) = mateC;
        mateC = NULL ;
    }
    if (mateR_handle != NULL)
    {
        (*mateR_handle) = mateR;
        mateR = NULL ;
    }
    LG_FREE_WORK;

    return (GrB_SUCCESS);
}
