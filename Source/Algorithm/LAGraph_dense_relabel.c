//------------------------------------------------------------------------------
// LAGraph_dense_relabel: dense relabeling of ids to matrix indices
//------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------

// LAGraph_dense_relabel: Contributed by Marton Elekes and Gabor Szarnyas,
// Budapest University of Technology and Economics
// (with accented characters: M\'{a}rton Elekes and G\'{a}bor Sz\'{a}rnyas.

// ...

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
}

//------------------------------------------------------------------------------

GrB_Info LAGraph_dense_relabel   // compute dense relabel
(
    GrB_Matrix *MMapping_handle, // output matrix with the mapping (unfilled if NULL)
    GrB_Vector *VMapping_handle, // output vector with the mapping (unfilled if NULL)
    const GrB_Index nids,        // number of identifiers
    const GrB_Index *ids         // array of identifiers
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (MMapping_handle == NULL && VMapping_handle == NULL)
    {
        LAGRAPH_ERROR ("Both mapping arguments are NULL", GrB_NULL_POINTER) ;
    }

    GrB_Info info ;

    // ... TODO

    return (GrB_SUCCESS) ;
}

