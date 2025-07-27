//------------------------------------------------------------------------------
// LAGraph_RpqMatrix: regular path query algortithm
//------------------------------------------------------------------------------

#define LG_FREE_WORK \
    {                \
    }

#define LG_FREE_ALL   \
    {                 \
        LG_FREE_WORK; \
    }

#include "LG_internal.h"
#include "LAGraphX.h"
#include <assert.h>

#define OK(s)                                               \
    {                                                       \
        GrB_Info info = s;                                  \
        if (!(info == GrB_SUCCESS))                         \
        {                                                   \
            printf("GraphBLAS error: %d\n", info);          \
            fprintf(stderr, "GraphBLAS error: %d\n", info); \
        }                                                   \
    }

// typedef enum RpqMatrixOp
// {
//     RPQ_MATRIX_OP_LABEL,
//     RPQ_MATRIX_OP_LOR,
//     RPQ_MATRIX_OP_CONCAT,
//     RPQ_MATRIX_OP_KLEENE,
//     RPQ_MATRIX_OP_KLEENE_L,
//     RPQ_MATRIX_OP_KLEENE_R,
// } RpqMatrixOp;

// typedef struct RpqMatrixPlan
// {
//     RpqMatrixOp op;
//     struct RpqMatrixPlan *lhs;
//     struct RpqMatrixPlan *rhs;
//     GrB_Matrix mat;
//     GrB_Matrix res_mat;
// } RpqMatrixPlan;

static GrB_Semiring sr;
static GrB_Monoid op;

GrB_Info LAGraph_RpqMatrix(RpqMatrixPlan *plan, char *msg);

static GrB_Info LAGraph_RpqMatrixLor(RpqMatrixPlan *plan, char *msg)
{
    LG_ASSERT(plan != NULL, GrB_NULL_POINTER);
    LG_ASSERT(plan->op == RPQ_MATRIX_OP_LOR, GrB_INVALID_VALUE);
    LG_ASSERT(plan->res_mat == NULL, GrB_INVALID_VALUE);

    RpqMatrixPlan *lhs = plan->lhs;
    RpqMatrixPlan *rhs = plan->rhs;

    LG_ASSERT(lhs != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER);

    OK(LAGraph_RPQMatrix(lhs, msg));
    OK(LAGraph_RPQMatrix(rhs, msg));

    GrB_Matrix lhs_mat = lhs->res_mat;
    GrB_Matrix rhs_mat = rhs->res_mat;

    GRB_TRY(GrB_eWiseAdd(plan->res_mat, GrB_NULL, GrB_NULL,
                         op, lhs_mat, rhs_mat, GrB_DESC_R));

    return (GrB_SUCCESS);
}

static GrB_Info LAGraph_RpqMatrixConcat(RpqMatrixPlan *plan, char *msg)
{

    LG_ASSERT(plan != NULL, GrB_NULL_POINTER);
    LG_ASSERT(plan->op == RPQ_MATRIX_OP_CONCAT, GrB_INVALID_VALUE);
    LG_ASSERT(plan->res_mat == NULL, GrB_INVALID_VALUE);

    RpqMatrixPlan *lhs = plan->lhs;
    RpqMatrixPlan *rhs = plan->rhs;

    GrB_Index nvalsA, nvalsB;
    GrB_Matrix_nvals(&nvalsA, lhs->mat);
    GrB_Matrix_nvals(&nvalsB, rhs->mat);
    fprintf(stderr, "\nDEBUG: A:%lu and B:%lu in Concat\n", nvalsA, nvalsB);

    LG_ASSERT(lhs != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER);

    OK(LAGraph_RPQMatrix(lhs, msg));
    OK(LAGraph_RPQMatrix(rhs, msg));

    GrB_Matrix lhs_mat = lhs->res_mat;
    GrB_Matrix rhs_mat = rhs->res_mat;

    GrB_Matrix_nvals(&nvalsA, lhs_mat);
    GrB_Matrix_nvals(&nvalsB, rhs_mat);
    fprintf(stderr, "\nDEBUG: A:%lu and B:%lu in Concat after traversal\n", nvalsA, nvalsB);
    fprintf(stderr, "\nDEBUG: before mxm\n");

    GrB_Index width, height;
    GrB_Matrix_nrows(&height, rhs_mat);
    GrB_Matrix_ncols(&width, lhs_mat);
    GrB_Matrix res;
    GrB_Matrix_new(&res, GrB_BOOL, height, width);
    GRB_TRY(GrB_mxm(res, GrB_NULL, GrB_NULL,
                    sr, lhs_mat, rhs_mat, GrB_DESC_R));
    // GrB_mxm(plan->res_mat, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, lhs_mat, rhs_mat, GrB_NULL);
    plan->res_mat = res;
    GrB_Matrix_nvals(&nvalsA, plan->res_mat);
    fprintf(stderr, "\nDEBUG: A:%lu after mxm in Concat\n", nvalsA);

    return (GrB_SUCCESS);
}

static GrB_Info LAGraph_RpqMatrixKleene(RpqMatrixPlan *plan, char *msg)
{
    LG_ASSERT(plan != NULL, GrB_NULL_POINTER);
    LG_ASSERT(plan->op == RPQ_MATRIX_OP_KLEENE, GrB_INVALID_VALUE);
    LG_ASSERT(plan->res_mat == NULL, GrB_INVALID_VALUE);

    RpqMatrixPlan *lhs = plan->lhs;
    RpqMatrixPlan *rhs = plan->rhs;

    // Kleene star should have one child. Always right.
    LG_ASSERT(lhs == NULL, GrB_NULL_POINTER);
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER);

    OK(LAGraph_RPQMatrix(rhs, msg));

    GrB_Matrix B = rhs->res_mat;

    // Creating identity matrix.
    GrB_Index n;
    GRB_TRY(GrB_Matrix_nrows(&n, B));
    GrB_Matrix E;
    GRB_TRY(GrB_Matrix_new(&E, GrB_BOOL, n, n));

    GrB_Vector v;
    GRB_TRY(GrB_Vector_new(&v, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_assign_BOOL(v, NULL, NULL, true, GrB_ALL, n, NULL));

    GRB_TRY(GrB_Matrix_diag(&E, v, 0));

    GRB_TRY(GrB_Vector_free(&v));

    // B + E
    GrB_Matrix BPE;
    GRB_TRY(GrB_Matrix_new(&BPE, GrB_BOOL, n, n));
    GRB_TRY(GrB_eWiseAdd(BPE, GrB_NULL, GrB_NULL,
                         op, B, E, GrB_DESC_R));
    // S <- S x (B + E)
    GrB_Matrix S, T;
    GRB_TRY(GrB_Matrix_dup(&S, E));

    bool changed = true;
    GrB_Index nnz_S = n, nnz_T = 0;
    while (changed)
    {
        // T = S * (B + E)
        GRB_TRY(GrB_Matrix_new(&T, GrB_BOOL, n, n));
        GRB_TRY(GrB_mxm(T, GrB_NULL, GrB_NULL,
                        sr, S, BPE, GrB_DESC_R));

        GRB_TRY(GrB_Matrix_nvals(&nnz_T, T));
        if (nnz_T != nnz_S)
        {
            changed = true;
            nnz_S = nnz_T;
            GRB_TRY(GrB_Matrix_free(&S));
            S = T;
        }
        else
        {
            changed = false;
            GRB_TRY(GrB_Matrix_free(&T));
        }
    }
    plan->res_mat = S;

    GRB_TRY(GrB_Matrix_free(&E));
    GRB_TRY(GrB_Matrix_free(&BPE));
    return (GrB_SUCCESS);
}

GrB_Info LAGraph_RPQMatrix(RpqMatrixPlan *plan, char *msg)
{
    if (plan->res_mat != NULL)
    {
        GrB_Index result;
        GrB_Matrix_nvals(&result, plan->res_mat);
        fprintf(stderr, "\nDEBUG: res_mat in LAGraph_RPQMatrix: %lu", result);
        return (GrB_SUCCESS);
    }

    switch (plan->op)
    {
    case RPQ_MATRIX_OP_LABEL:
        LG_ASSERT(plan->lhs == NULL && plan->rhs == NULL, GrB_INVALID_VALUE);
        GrB_Index result;
        GrB_Matrix_nvals(&result, plan->mat);
        fprintf(stderr, "\nDEBUG: res_mat in LAGraph_RPQMatrix switch: %lu", result);
        plan->res_mat = plan->mat;
        return (GrB_SUCCESS);
    case RPQ_MATRIX_OP_LOR:
        return LAGraph_RpqMatrixLor(plan, msg);
    case RPQ_MATRIX_OP_CONCAT:
        return LAGraph_RpqMatrixConcat(plan, msg);
    case RPQ_MATRIX_OP_KLEENE:
        return LAGraph_RpqMatrixKleene(plan, msg);
    default:
        LG_ASSERT(false, GrB_INVALID_VALUE);
    }
    GrB_Index result;
    GrB_Matrix_nvals(&result, plan->res_mat);
    fprintf(stderr, "\nDEBUG: res_mat in LAGraph_RPQMatrix end: %lu", result);
    return (GrB_SUCCESS);
}

GrB_Info LAGraph_RpqMatrix_initialize()
{
    sr = LAGraph_any_one_bool;
    op = GxB_ANY_BOOL_MONOID;
}