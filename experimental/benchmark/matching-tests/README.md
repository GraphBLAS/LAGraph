This folder contains custom testing codes for LAGraph_MaximalMatching (go to file for further description):

## C++ programs:
* `gen_bipartite`: Builds a random undirected bipartite graph w/ user parameters and produces exact (maximum) or approximate (naive) matchings, and optional performance data
* `gen_general`: Builds a random undirected general (non-bipartite) graph w/ user parameters and produces exact (maximum) or approximate (naive) matchings, and optional performance data
* `verify_matching`: Used to verify that a GraphBLAS produced matching is correct. For internal use by `bench.py`
* `Makefile`: run `make CC=gcc CXX=g++` to build aforementioned programs

## Python programs:
* `bench.py`: Tool to auto-run list of tests

## To use:
* Make sure SuiteSparse:GraphBLAS and LAGraph are installed and built appropriately (see the top-level `README.md` for details)
* Install Python dependencies with `pip install -r requirements.txt`
* Run `make CC=gcc CXX=g++` to build the C++ executables above
* You can now run `bench.py`, or individually run one of the C++ executables from `build/` as needed