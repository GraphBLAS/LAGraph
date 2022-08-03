LAGraph Documentation
=====================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

:ref:`genindex`

GraphBLAS Objects
-----------------
.. doxygenstruct:: LAGraph_Graph_struct
   :members:

.. doxygenclass:: LAGraph_Graph
   :members:

Example
~~~~~~~

.. code-block:: C++

   #include <iostream>
   #include <grb/grb.hpp>
   int main(int argc, char** argv) {
     // Create a new matrix, reading in from a file.
     grb::matrix<float, int> a("data/chesapeake.mtx");

     size_t m = a.shape()[0];
     size_t k = a.shape()[1];

     std::cout << "chesapeake.mtx is a " << m << " by " << k << " matrix." << std::endl;

     // Set element 12,9 (row 12, column 9) to 12.
     a[{12, 9}] = 12;

     grb::matrix<float, int> b("data/chesapeake.mtx");

     auto c = grb::multiply(a, b);

     std::cout << "Sum of elements is " << grb::sum(c) << std::endl;

     return 0;
   }

Binary Operators
----------------
Binary operators are function objects that implement binary operators, that is
operators that accept two inputs and produce a single output.  A collection of
binary operators are pre-defined by GraphBLAS.

.. doxygentypedef:: grb::plus

.. doxygentypedef:: grb::minus

.. doxygentypedef:: grb::multiplies

.. doxygentypedef:: grb::times

.. doxygentypedef:: grb::max

.. doxygentypedef:: grb::min

.. doxygentypedef:: grb::modulus

Algorithms
----------
.. doxygenfunction:: grb::multiply(AMatrixType&&, BMatrixType&&, ReduceFn&&, CombineFn&&, MaskType&&)

.. doxygenfunction:: LAGr_SingleSourceShortestPath
    :project: lagraph

.. doxygenfunction:: grb::ewise(const AMatrixType&, const BMatrixType&, const BinaryOp&)

.. doxygenfunction:: LAGraph_Matrix_IsEqual
    :project: lagraph

Utility Functions
-----------------
.. doxygenfunction:: grb::print(const grb::matrix<Args...>&, std::string)