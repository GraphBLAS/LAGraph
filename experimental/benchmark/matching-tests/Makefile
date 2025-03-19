all:
	g++ -o build/gen_bipartite gen_bipartite.cpp -fPIC -fopenmp -I../../../include -I../../../../GraphBLAS/Include -L../../../build/src -L../../../build/experimental -L../../../../GraphBLAS/build -L../../../../GraphBLAS/build/CUDA -L../../../../GraphBLAS/build/rmm_wrap -l:liblagraph.a -l:liblagraphx.a -l:libgraphblas.a && \
	g++ -o build/gen_general gen_general.cpp -fPIC -fopenmp -I../../../include -I../../../../GraphBLAS/Include -L../../../build/src -L../../../build/experimental -L../../../../GraphBLAS/build -L../../../../GraphBLAS/build/CUDA -L../../../../GraphBLAS/build/rmm_wrap -l:liblagraph.a -l:liblagraphx.a -l:libgraphblas.a && \
	g++ -o build/verify_matching verify_matching.cpp --std=c++17 -fPIC -fopenmp -I../../../include -I../../../../GraphBLAS/Include -I../../../GraphBLAS/rmm_wrap -L../../../build/src -L../../../build/experimental -L../../../../GraphBLAS/build -l:liblagraph.a -l:liblagraphx.a -l:libgraphblas.a

clean:
	rm -rf build/*
