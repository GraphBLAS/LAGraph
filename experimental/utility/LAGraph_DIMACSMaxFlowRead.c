//------------------------------------------------------------------------------
// LAGraph_DIMACSMaxFlowRead: read a DIMACS13 MaxFlow problem
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2026 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Darin Peries and Tim Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

#define BUFF_SIZE 1024

#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define LG_FREE_WORK                            \
{                                               \
    LAGraph_Free ((void **) &rows, NULL) ;      \
    LAGraph_Free ((void **) &cols, NULL) ;      \
    LAGraph_Free ((void **) &weights, NULL) ;   \
}

#define LG_FREE_ALL                             \
{                                               \
    LG_FREE_WORK ;                              \
    GrB_free (A) ;                              \
}

int LAGraph_DIMACSMaxFlowRead
(
    // output:
    GrB_Matrix* A,  // adjancency matrix, with int32 weights
    GrB_Index* s,   // source node
    GrB_Index* t,   // sink node
    // input:
    FILE* f,        // an open file containing the DIMAX MaxFlow problem
    char* msg
)
{

  LG_CLEAR_MSG ;

  int64_t n_nodes = 0, n_edges = 0 ;
  GrB_Index *rows = NULL, *cols = NULL ;
  int32_t *weights = NULL ; // For DIMACS network flow: use 32 bit integer weights
  
  char buff[BUFF_SIZE+1];

  int64_t line_count = 0 ;

  LG_ASSERT (f != NULL && A != NULL && s != NULL && t != NULL, GrB_NULL_POINTER) ;
  (*A) = NULL ;

  bool p_present = false ;
  bool s_present = false ;
  bool t_present = false ;

  while (fgets(buff, BUFF_SIZE, f))
  {
    buff [BUFF_SIZE] = '\0' ;   // ensure buff is nul-terminated

    if (buff[0] == 'c' || buff[0] == '\0' || buff[0] == '\n')
    {

        // skip comments and empty lines
        continue;

    }
    else if (buff[0] == 'p')
    {

      // problem statement: p max <n_nodes> <n_edges>
      LG_ASSERT (!p_present, GrB_INVALID_VALUE) ;
      p_present = true ;
      int result = scanf (buff, "p max %" PRIu64 " %" PRIu64, &n_nodes, &n_edges) ;
      LG_ASSERT (result == 2, GrB_INVALID_VALUE) ;
      LG_TRY (LAGraph_Malloc ((void **) &rows, n_edges, sizeof (GrB_Index), msg)) ;
      LG_TRY (LAGraph_Malloc ((void **) &cols, n_edges, sizeof (GrB_Index), msg)) ;
      LG_TRY (LAGraph_Malloc ((void **) &weights, n_edges, sizeof (int32_t), msg)) ;

    }
    else if (buff[0] == 'n')
    {

      // source node: n <source> s
      // sink node:   n <sink> t
      GrB_Index value = 0 ;
      char which = ' ';

      int result = scanf(buff, "n %" PRIu64 ", %c", &value, which) ;
      LG_ASSERT (result == 2, GrB_INVALID_VALUE) ;

      switch (which) {
      case 's':
        *s = value ;
        LG_ASSERT (!s_present, GrB_INVALID_VALUE) ;
        s_present = true ;
        break ;
      case 't':
        *t = value ;
        LG_ASSERT (!t_present, GrB_INVALID_VALUE) ;
        t_present = true ;
        break;
      default:
        LG_ASSERT (false, GrB_INVALID_VALUE) ;
        break;
      }
      LG_ASSERT (value < n_nodes, GrB_INVALID_VALUE) ;

    }
    else if (buff[0] == 'a')
    {
      // a single edge (i,j) with weight w in the graph: a <i> <j> <w>
      LG_ASSERT (line_count < n_edges, GrB_INVALID_VALUE) ;
      GrB_Index r = 0, c = 0, w = 0 ;
      int result = scanf(buff, "a %" PRIu64 " %" PRIu64 " %" PRId32, &r, &c, &w) ;
      LG_ASSERT (result == 3, GrB_INVALID_VALUE) ;
      rows[line_count] = r ;
      cols[line_count] = c ;
      weights[line_count] = w ;
      line_count++ ;
    }
  }

  // the file must contain "p max..." line, "n <src> s" and "n <sink> t"
  LG_ASSERT (p_present && s_present && t_present, GrB_INVALID_VALUE) ;
  // the file must contain exactly n_edges
  LG_ASSERT (line_count == n_edges, GrB_INVALID_VALUE) ;

  // build the adjacency matrix; no duplicates can be present
  GRB_TRY (GrB_Matrix_new (A, GrB_INT32, n_nodes, n_nodes)) ;
  GRB_TRY(GrB_Matrix_build_INT32(*A, rows, cols, weights, n_edges,
    /* NULL dup operator means no duplicates are tolerated: */ NULL)) ;

  LG_FREE_WORK ;
  return (GrB_SUCCESS) ;
}
