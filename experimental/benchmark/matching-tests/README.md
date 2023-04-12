This folder contains custom testing codes for LAGraph_MaximalMatching (go to file for further description):

## C++ programs:

* `gen_bipartite`: Builds a random undirected bipartite graph w/ user parameters and produces exact (maximum) or approximate (naive) matchings
* `gen_general`: Builds a random undirected general (non-bipartite) graph w/ user parameters and produces exact (maximum) or approximate (naive) matchings
* `verify_matching`: Used to verify that a GraphBLAS produced matching is correct
* `Makefile`: run `make CC=gcc CXX=g++` to build aforementioned programs

## Python programs:
* `bench.py`: Tool to auto-run list of tests