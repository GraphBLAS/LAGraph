
To replicate the sync. bug:

compile GraphBLAS from its cuda_bug branch of github,
ideally on a system with a V100 GPU.

compile LAGraph from its cuda_bug branch.

run the following:

    source go_10 ; tail out_0010.txt

