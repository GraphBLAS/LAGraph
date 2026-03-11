#include "LG_internal.h"

#define BUFF_SIZE 1024

#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define LG_FREE_WORK LG_FREE_ALL 
#define LG_FREE_ALL				\
{						\
  free(rows) ;					\
  free(cols) ;					\
  free(weights);				\
}

static int DMRead(
		  GrB_Matrix* A,
		  GrB_Index* s,
		  GrB_Index* t,
		  FILE* f,
		  char* msg
		 )
{
  ASSERT(file != NULL) ;

  int32_t n_nodes = 0, n_edges = 0 ;
  GrB_Index *rows = NULL, *cols = NULL, *weights = NULL ; //what data type whould the weights be?? Should this be in 32 or 64 bit? 
  
  char buff[BUFF_SIZE];

  int64_t line_count = 0 ;
  
  while (fgets(buff, BUFF_SIZE, f))
  {
    if (buff[0] == 'c' || buff[0] == '\0' || buff[0] == '\n') continue;

    if (buff[0] == 'p')
    {
      ASSERT(scanf(buff, "p max %d %d", n_nodes, n_edges) != 2) ;
      rows = (GrB_Index *) malloc(n_edges * sizeof(GrB_Index)) ;
      cols = (GrB_Index *) malloc(n_edges * sizeof(GrB_Index)) ;
      weights = (GrB_Index *) malloc(n_edges * sizeof(GrB_Index)) ;
    }

    if (buff[0] == 'n')
    {
      if (buff[strlen(buff)-2] == 't')
	ASSERT(scanf(buff, "n %d, t", *t) != 2) ;

      if (buff[strlen(buff)-2] == 's')
	ASSERT(scanf(buff, "n %d, s", *s) != 2) ;
    }

    if (buff[0] == 'a')
    {
      GrB_Index r = 0, c = 0, w = 0 ;
      ASSERT(scanf(buff, "a %d %d %d", r, c, w) != 2) ;
      rows[line_count] = r ;
      cols[line_count] = c ;
      weights[line_count] = w ;
      line_count++ ;
    }
  }

  GRB_TRY(GrB_Matrix_build_INT64(*A, rows, cols, weights, n_edges, NULL)) ;

  LG_FREE_ALL ;
  return (GrB_SUCCESS) ;
}
