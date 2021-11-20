
//****************************************************************************
// Unused and/or deprecated methods
//****************************************************************************

GrB_Info LAGraph_log
(
    char *caller,           // calling function
    char *message1,         // message to include (may be NULL)
    char *message2,         // message to include (may be NULL)
    int nthreads,           // # of threads used
    double t                // time taken by the test
) ;

GrB_Info LAGraph_1_to_n     // create an integer vector v = 1:n
(
    GrB_Vector *v_handle,   // vector to create
    GrB_Index n             // size of vector to create
);

GrB_Info LAGraph_ispattern  // return GrB_SUCCESS if successful
(
    bool *result,           // true if A is all one, false otherwise
    GrB_Matrix A,
    GrB_UnaryOp userop      // for A with arbitrary user-defined type.
                            // Ignored if A and B are of built-in types
);

GrB_Info LAGraph_isall      // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp op         // GrB_EQ_<type>, for the type of A and B,
                            // to check for equality.  Or use any desired
                            // operator.  The operator should return GrB_BOOL.
);

GrB_Info LAGraph_Vector_isall      // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Vector A,
    GrB_Vector B,
    GrB_BinaryOp op         // GrB_EQ_<type>, for the type of A and B,
                            // to check for equality.  Or use any desired
                            // operator.  The operator should return GrB_BOOL.
);

GrB_Info LAGraph_Vector_to_dense
(
    GrB_Vector *vdense,     // output vector
    GrB_Vector v,           // input vector
    void *id                // pointer to value to fill vdense with
);

