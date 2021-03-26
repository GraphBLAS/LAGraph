Results on Intel DevCloud: 

Each row is the result from a single program that loads the urand
matrix once, and then runs the GAP TC benchmark in GraphBLAS using
6 different algorithms.  Each algorithm is run 3 times, and the
average is reported (see the qsub output for all 3 times; not listed
here).  The time to compute L=tril(A,-1) and U=triu(A,1) is included
in all run times, but not the time to load A from the input file.
Four methods make a copy of A, sorting it by ascending (+) or
descending(-) order of degree.  The 2 methods that do not sort are
fastest, since this matrix does not require sorting.  The Dot
method computes C<L>=L*U', the Dot2 method computes C<U>=U*L';
all matrices in CSR format.

The qsub jobs were submitted on as many nodes of the DevCloud as possible
(more pending).  Each row is a separate run on the same node.

The best time of the 6 algorithms was found, for all runs and all nodes:

best run times:
Dot:-  : 109.224
Dot:0  : 60.357
Dot:+  : 108.014
Dot2:- : 108.867
Dot2:0 : 62.894
Dot2:+ : 111.169

The run times below were then normalized, with t = (t/tbest-1)*100.
That is, if the result below is "3" than that run took an amount of
time that is 3% more than the best time found.

For all the runs on a single node, the best (min), average, and worst
(max) result is listed, and then for all 6 algorithms, the overall
min, avg, and max is reported on the right ("overall: ...") ;



