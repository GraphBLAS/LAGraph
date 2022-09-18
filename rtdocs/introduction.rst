Introduction
============

A graph is a set of vertices and a set of edges between them.  This pair of sets
sets leads directly to the familiar picture of a graph as a set of dots connected
by arcs (an undirected graphs) or arrows (a directed graph). You can also represent a 
graph in terms of matrices.   Usually, this is done with an adjacency matrix where
the rows and columns correspond to the vertices and the non-empty
elements represent the edges between vertices. Since fully connected graphs (i.e., every
vertex is connected to every other vertex) are rare, matrices used for Graphs are
typically sparse (most elements are empty so it makes sense to only store the non-empty 
elements).

Representing a graph as a sparse matrix results in graph algorithms expressed in 
terms of linear algebra.    For example, if a vector represents a set of vertices,
multiplication of that vector by an adjacency matrix returns a vector of the neighbors of those vertices.
A sequence of multiplications traverses the graph in a breadth first manner.  

To cover the full range of graph algorithms, one additional ingredient is needed.  We are used to 
thinking of matrix operations over real numbers: multiple pairs of matrix elements and then
combine the resulting products through addition.  There are times, however, when those operations do
not provide the functionality needed by an algorithm. For example, it may be better to combine elements by
only keeping the minimum value.   The elements of the matrices may be Boolean values or integers or
even a user-defined type.  If the goal is to cover the full range of graph algorithms, therefore,
we need a way to generalize the type and the operators to use instead of the usual addition and multiplication.

We do this through an algebraic semiring.   This algebraic structure consists of (1) an operator
corresponding to addition, (2) the identity of that operator, (3) an operator corresponding to multiplication,
and (4) the identity of that operator.  We are all familiar with the semiring used with real numbers
consisting of $(+,0,*,1)$.  A common semiring in graph algorithms is the so-called tropical semiring 
consisting of $(min,infinity,+,0)$.  This is used in shortest path algorithms.   These semirings 
give us a mathematically rigorous way to modify the operators used in our graph algorithms.

If you work with linear algebra, you most likely know about the Basic Linear Algebra subprograms or BLAS.
Introduced in the 70's and 80's, the BLAS had a huge impact on the practice of linear algebra.  By designing
linear algebra in terms of the BLAS, an algorithm can be expressed at a high level leaving specialization to the 
low level details of a particular hardware platform to the BLAS.  So if you want to use Linear Algebra for 
Graph Algorithms, it stands to reason that you need the Basic Linear Algebra Subprograms for Graph Algorithms.
We call these the GraphBLAS (www.graphblas.org).

The GraphBLAS define opaque types for a matrix and a vector objects.  Since these objects are opaque, an implementation
has the freedom to specialize the data structures as needed to optimzie the software for a particular platform.  The
graphBLAS are great for people interested in sparse linear algebra and designing their own graph algorithms.
The GraphBLAS library, however, does not include any graph algorithms.   The GraphBLAS provide a software framework for
constructing graph algorithms, but it doesn't provide any actual Graph Algorithms.  Since most people working with
graphs use algorithms but don't develop them "from scratch", the graphBLAS are not really useful to most people.

Hence, there is a need for a library of Graph Algorithms implemented on top of the GraphBLAS.  This
library is called LAGraph.   

discuss the LAGraph idea:  graph algorithms using GraphBLAS,
with several additional data structures.  In particular the LAGraph_Graph.
Discuss its properties, and cached properties, and the differences
between Basic (``LAGraph_*``) and Advanced (``LAGr_*``).

LAGraph and GraphBLAS work together; you can have a mix of
LAGraph_Graph objects and GraphBLAS objects together in the same
program.  LAGraph methods often return GrB_Matrix or GrB_Vector
results.  An LAGraph_Graph can be modified by the user application
directly (it is not an opaque object), by including GrB_Matrix or
GrB_Vector objects into it.

