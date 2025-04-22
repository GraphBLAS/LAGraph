
To replicate the sync bug:

compile GraphBLAS from its cuda_bug branch of github,
ideally on a system with a V100 GPU.

compile LAGraph from its cuda_bug branch.

run the following:

    source go_10 ; tail out_0010.txt


details:


    mkdir cuda_bug
    cd cuda_bug
    git clone https://github.com/GraphBLAS/LAGraph.git
    cd LAGraph
    git checkout cuda_bug
    cd ..
    git clone https://github.com/DrTimothyAldenDavis/GraphBLAS.git
    cd GraphBLAS
    git checkout cuda_bug
    make JOBS=32
    cd ../LAGraph
    make JOBS=32
    ./go_10
    tail out_0010.txt


to run all the matrices, you first need to edit the file
GraphBLAS/CUDA/template/GB_cuda_jit_AxB_dot3_phase2.cuh, and change the
"#if 1" on line 320 to "#if 0".  Then recompile GraphBLAS as above.  Then come
back to the LAGraph folder and run this:

    make JOBS=32
    ./go156

and look at all the out*.txt files:

    ack Abort out*txt
