LAGraph/Test/Gr_Read/README.txt

As in all LAGraph/Test directories, just typing:

    make

will compile and run the test.

See make.out for the expected output.

There are two test files in the Gr_matrix folder.  Both files t1.gr and t2.gr
are unweighted, with an edge weight size of zero.  Thus, the following tests
are used:

	./build/gtest Gr_matrix/t1.gr
	./build/gtest Gr_matrix/t2.gr

t1 is the following matrix, where '.' denote no entry and 'x' is an edge.

    0 1 2 3 4 5
  0 . x . . x x
  1 x . x x . x
  2 . x . x . .
  3 . x x . x .
  4 x . . x . .
  5 x x . . . .

The printout of the matrix from GraphBLAS is the following:

GraphBLAS matrix: G 
nrows: 6 ncols: 6 max # entries: 16
format: standard CSR vlen: 6 nvec: 6 plen: 6 vdim: 6
hyper_ratio 0.0625
GraphBLAS type:  bool size: 1
number of entries: 16 
row: 0 : 3 entries [0:2]
    column 1: bool 1
    column 4: bool 1
    column 5: bool 1
row: 1 : 4 entries [3:6]
    column 0: bool 1
    column 2: bool 1
    column 3: bool 1
    column 5: bool 1
row: 2 : 2 entries [7:8]
    column 1: bool 1
    column 3: bool 1
row: 3 : 3 entries [9:11]
    column 1: bool 1
    column 2: bool 1
    column 4: bool 1
row: 4 : 2 entries [12:13]
    column 0: bool 1
    column 3: bool 1
row: 5 : 2 entries [14:15]
    column 0: bool 1
    column 1: bool 1


    0 1 2 3 4 5
  0 . x . . x x
  1 x . x x . x
  2 . x . x . .
  3 . x x . x .
  4 x . . x . .
  5 x x . . . .

The t2.gr matrix is the following matrix.  Note that it is equal to the
upper triangular part of t1.  That is, t2 = triu(t1) in MATLAB notation.

    0 1 2 3 4 5
  0 . x . . x x
  1 . . x x . x
  2 . . . x . .
  3 . . . . x .
  4 . . . . . .
  5 . . . . . .

GraphBLAS matrix: G 
nrows: 6 ncols: 6 max # entries: 8
format: standard CSR vlen: 6 nvec: 6 plen: 6 vdim: 6
hyper_ratio 0.0625
GraphBLAS type:  bool size: 1
number of entries: 8 
row: 0 : 3 entries [0:2]
    column 1: bool 1
    column 4: bool 1
    column 5: bool 1
row: 1 : 3 entries [3:5]
    column 2: bool 1
    column 3: bool 1
    column 5: bool 1
row: 2 : 1 entries [6:6]
    column 3: bool 1
row: 3 : 1 entries [7:7]
    column 4: bool 1

