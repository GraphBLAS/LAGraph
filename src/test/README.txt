
To run the tests:

    cd LAGraph/build
    cmake .. ; make ; make test

To run with test coverage:

    cd LAGraph/build
    cmake -DCOVERAGE=1 .. ; make ; make test_coverage

On the Mac, you should use gcc-11, from "brew install gcc" and this command:

    CC=gcc-11 CXX=g++-11 GCOV_PATH=/usr/local/lib/gcov-11 cmake -DCOVERAGE=1 .. ; make ; make test_coverage

