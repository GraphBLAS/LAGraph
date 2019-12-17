/*
================
LAGraph
Copyright 2019 LAGraph Contributors. 
(see Contributors.txt for a full list of Contributors; see ContributionInstructions.txt for information on how you can Contribute to this project). 
All Rights Reserved.
NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED, AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.
Released under a BSD license, please see the License.txt distributed with this Software or contact permission@sei.cmu.edu for full terms.
Created, in part, with funding and support from the United States Government. (see Acknowledgments file).
This program includes and/or can make use of certain third party source code, object code, documentation and other files ("Third Party Software"). See License.txt file for more details.
DM19-0270

U. Sridhar, M. Blanco, R. Mayuranath, D. G. Spampinato, T. M. Low, and S. McMillan, “Delta-Stepping SSSP: From Vertices and Edges to GraphBLAS Implementations,” in 2019 IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW), 2019, pp. 241–250.

https://ieeexplore.ieee.org/document/8778222/references
https://arxiv.org/abs/1911.06895
*/

#include "GraphBLAS.h"
#include "read_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define LIMIT 100
#define RUNS 500

unsigned long long rdtsc()
{
  unsigned long long int x;
  unsigned a, d;

  __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));

  return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

/*
Algorithm for delta stepping sssp:
AH = A◦(A>∆)
AL = A◦(A<=∆)
t = ∞
t[σ] = 0
i = 0

while (t ≥ i∆) ̸= 0 do
  s=0
  tB =(i∆≤t<(i+1)∆)
  while tB ̸= 0 do
    tReq = ATL(t◦tB)
    s = s+tB, tB =0
    tB = (i∆ ≤ tReq < (i+1)∆)◦(tReq < t)
    t = min(t,tReq)
  od
  tReq = ATH(t◦s)
  t=min(t,tReq)
  i=i+1
od
*/


typedef double type_t;
#define T_BASE FP64
#define T GrB_FP64  
#define GRB_MIN_T GrB_MIN_FP64  
#define GRB_PLUS GrB_PLUS_FP64  
#define GRB_ID GrB_IDENTITY_FP64  
#define GRB_LT_T GrB_LT_FP64  


// Globals necessary for thresholding on delta and i:
uint64_t i;
type_t delta;
type_t lim_iters;


// comparison functions for delta (to be used in GrB unary operator:
// for  z = x <= y 
void leq_delta(void * z, const void * x){
	(*(bool*)z) = (bool)((*(type_t*)x) <= delta);
} 

void gt_delta(void * z, const void * x){
	(*(bool*)z) = (bool)((*(type_t*)x) > delta);
} 

// comparison funtions for delta_i (GrB Unary Op):
void in_delta_i_range(void * z, const void * x){
	type_t delta_i = delta * i;
	type_t delta_i_1 = delta * (i+1);
	(*(bool*)z) = (bool)(delta_i <= (*(type_t*)x)) && ((*(type_t*)x) < delta_i_1);
}

void geq_delta_i(void * z, const void * x){
	type_t delta_i = delta * i;
	(*(bool*)z) = (bool)(delta_i <= (*(type_t*)x));
}

GrB_Info sssp_delta_step(const GrB_Matrix A, type_t d, 
	GrB_Index src, GrB_Vector ** paths){
	// Global scalars:
	delta=d;
	i=0;
	// Setup operators and semirings:
	GrB_UnaryOp delta_i_range = NULL;
	GrB_UnaryOp delta_leq = NULL;
	GrB_UnaryOp delta_gt = NULL;
	GrB_UnaryOp delta_i_geq = NULL;
	GrB_Monoid min_monoid = NULL;
	GrB_Semiring min_plus_sring = NULL;
	GrB_UnaryOp_new(&delta_i_range, &in_delta_i_range, GrB_BOOL, T);
	GrB_UnaryOp_new(&delta_leq, &leq_delta, GrB_BOOL, T);
	GrB_UnaryOp_new(&delta_gt, &gt_delta, GrB_BOOL, T);
	GrB_UnaryOp_new(&delta_i_geq, &geq_delta_i, GrB_BOOL, T);
	GrB_Monoid_new(&min_monoid, GRB_MIN_T, (double)INFINITY)
	GrB_Semiring_new(&min_plus_sring, min_monoid, GrB_PLUS_FP64);

	// Setup Scalars
	GrB_Index n;
	GrB_Index m;
	GrB_Index t_size;
	GrB_Index tm_size;
	GrB_Matrix_nrows(&n, A);
	GrB_Matrix_ncols(&m, A);

	lim_iters = LIMIT * n / delta;

	// Non-boolean vectors:
	// t = ∞
	GrB_Vector t;
	GrB_Vector tReq;
	GrB_Vector tmasked;
	GrB_Vector_new(&t, T, n);	
	GrB_Vector_new(&tReq, T, n);
	GrB_Vector_new(&tmasked, T, n);

	// Boolean vectors:
	GrB_Vector tB;
	GrB_Vector tcomp;
	GrB_Vector tless;
	GrB_Vector tless1;
	GrB_Vector s;

	// Bool vectors init to false
	GrB_Vector_new(&tB, GrB_BOOL, n); 
	GrB_Vector_new(&tcomp, GrB_BOOL, n);
	GrB_Vector_new(&tless, GrB_BOOL, n);
	GrB_Vector_new(&tless1, GrB_BOOL, n);
	GrB_Vector_new(&s, GrB_BOOL, n);

	// t[σ] = 0
	GrB_Vector_setElement(t, 0, src);

	// Create high and low matrices from A based on delta:
	GrB_Matrix Ah = NULL;
	GrB_Matrix Al = NULL;
	GrB_Matrix Ab = NULL;

	GrB_Matrix_new(&Ah, T, n, m);
	GrB_Matrix_new(&Al, T, n, m);
	GrB_Matrix_new(&Ab, GrB_BOOL, n, m);

	// AL = A◦(A<=∆)
	GrB_apply(Ab, GrB_NULL, GrB_NULL, delta_leq, A, GrB_NULL );
	GrB_apply(Al, Ab, GrB_NULL, GRB_ID, A, GrB_NULL );

	// AH = A◦(A>∆)
	GrB_apply(Ab, GrB_NULL, GrB_NULL, delta_gt, A, GrB_NULL );
	GrB_apply(Ah, Ab, GrB_NULL, GRB_ID, A, GrB_NULL );

	// init i = 0
	i = 0;

	// t >= i*delta
	// Actually tcomp = t > delta*i
	GrB_Vector_apply(tless1, GrB_NULL, GrB_NULL, delta_i_geq, t, GrB_NULL);
	GrB_Vector_apply(tcomp, tless1, GrB_NULL, GrB_IDENTITY_BOOL, t, GrB_NULL);


	// while (t ≥ i∆) ̸= 0 do
	GrB_Vector_nvals(&t_size, tcomp);
	while ( t_size > 0 && i < lim_iters){
		//s=0
		GrB_Vector_clear(s);

		// Overall operation is tmasked = t.*(delta_i <= t < delta_i_1) --> tm<tB> = t
		// tB = (bool) (delta_i <= t < delta_i_1)
  		GrB_Vector_apply(tB, GrB_NULL, GrB_NULL, delta_i_range, t, GrB_NULL);
		// tmasked = t◦tB:
		// tmasked is dest. Mask and accumulator binOp are null.
		// Unary op is delta_i_range, input vector is t, and descriptor is NULL.
  		GrB_Vector_apply(tmasked, tB, GrB_NULL, GrB_IDENTITY_FP64, t, GrB_NULL);

		//   while tB ̸= 0 do
		GrB_Vector_nvals(&tm_size, tmasked);

		while (tm_size > 0) {
			// tReq = ATL*tmasked over min-plus-semiring
			GrB_vxm(tReq, GrB_NULL, GrB_NULL, min_plus_sring, tmasked, Al, GrB_NULL);
			
			// s = s+tB, tB =0
			// Note: accumulator is NOT null.
			// GrB_Vector_apply(s, GrB_NULL, GrB_LOR, GrB_IDENTITY_BOOL, tB, GrB_NULL);
			GrB_eWiseAdd(s, GrB_NULL, GrB_NULL, GrB_LOR, s, tB, GrB_NULL);
			GrB_Vector_clear(tB);

			// tB = (i∆ ≤ tReq < (i+1)∆)◦(tReq < t) in three steps:
			// tless = tReq < t
			// tB = (bool) (delta_i <= tReq < (i+1)delta) and tless
			GrB_eWiseAdd(tless, tReq, GrB_NULL, GRB_LT_T, tReq, t, GrB_NULL);
			GrB_Vector_apply(tB, tless, GrB_NULL, GrB_IDENTITY_BOOL, tReq, GrB_NULL);
			
			// t = min(t,tReq)
			GrB_eWiseAdd(t, GrB_NULL, GrB_NULL, GRB_MIN_T, t, tReq, GrB_NULL);

			GrB_Vector_nvals(&tm_size, tB);
		}
		// tmasked = (type T) t◦s
		// tReq = ATH(tmasked)
		GrB_Vector_clear(tmasked);
		GrB_Vector_apply(tmasked, s, GrB_NULL, GRB_ID, t, GrB_NULL);
		GrB_Vector_clear(tReq);
		GrB_vxm(tReq, GrB_NULL, GrB_NULL, min_plus_sring, tmasked, Ah, GrB_NULL);

		// t = min(t,tReq)
		GrB_eWiseAdd(t, GrB_NULL, GrB_NULL, GRB_MIN_T, t, tReq, GrB_NULL);
		
		// i=i+1
		i++;

		// Recompute if there are remaining unreached nodes:
		// tcomp = t > delta*i
		GrB_Vector_clear(tless1);
		GrB_Vector_clear(tcomp);
		GrB_apply(tless1, GrB_NULL, GrB_NULL, delta_i_geq, t, GrB_NULL);
		GrB_apply(tcomp, tless1, GrB_NULL, GrB_IDENTITY_BOOL, tcomp, GrB_NULL);
		GrB_Vector_nvals(&t_size, tcomp);
	}
	// Set the return paths:
	(*paths) = &t;

	// Free everything except for t, which has been 
	// passed on to paths as a returned result.
	GrB_free(&delta_i_range);
	GrB_free(&delta_leq);
	GrB_free(&delta_gt);
	GrB_free(&delta_i_geq);
	GrB_free(&min_monoid);
	GrB_free(&min_plus_sring);
	GrB_free(&s);
	GrB_free(&Ab);
	GrB_free(&Al);
	GrB_free(&Ah);
	GrB_free(&tReq);
	GrB_free(&tmasked);
	GrB_free(&tcomp);
	GrB_free(&tless);
	GrB_free(&tless1);
	GrB_free(&tB);

	return (GrB_SUCCESS);
}


int main(int argc, char** argv){
	unsigned long long t0, t1, cycles;

	if (argc != 4){
		printf("USAGE: <exec> <Matrix Filename> <src_index> <delta_step>\n");
		exit(EXIT_FAILURE);
	}
	GrB_init(GrB_BLOCKING);
	// GrB_init(GrB_NONBLOCKING);
	
	GrB_Matrix graph = NULL;
	GrB_Vector * paths = NULL;
	FILE * f = fopen(argv[1], "r");
	uint64_t src = atol(argv[2]);
	double delta = atof(argv[3]);

	read_matrix(
		&graph,          	// handle of matrix to create
		f,               	// file to read the tuples from
		false,    			// if true, return A as symmetric
		true,     			// if true, then remove self edges from A
		false,         		// if true, input matrix is 1-based
		false,           	// if true, input is GrB_BOOL, otherwise GrB_FP64
		false         		// if true, print status to stdout
		);
	fclose(f);

	GrB_Index tm_size, n, nnz;
	GrB_Matrix_nrows(&n, graph);
	GrB_Matrix_nvals(&nnz, graph);

	// printf("Finished reading. Starting SSSP from node %lu with delta %f\n", src, delta);

	for (uint32_t itr = 0; itr < RUNS; itr++){
		t0 = rdtsc();
		sssp_delta_step(graph, delta, src, &paths);
		t1 = rdtsc();
		cycles += t1-t0;
	}

	// GxB_print(*paths, GxB_SUMMARY);
	GrB_Vector_nvals(&tm_size, *paths);

	char * term_cond = (i > lim_iters) ? "ITER_LIMIT" : "NORMAL_TERMINATION";

	// Name of matrix, end condition, nodes, edges, result nnz, avg cycles
	printf("%s,%s,%u,%u,%u,%f\n", 
		argv[1], term_cond, n, nnz, tm_size,((double)cycles)/RUNS);

	GrB_free(&graph);
	// GrB_free(paths); // Causes seg fault. Why?
	GrB_finalize();
	return 0;
}
