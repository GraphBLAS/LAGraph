################################################################################
################################################################################
##                                                                            ##
##  To benchmark LAGraph: please contact Tim Davis at davis@tamu.edu first.   ##
##                                                                            ##
################################################################################
################################################################################

LAGraph is a draft package, and its performance is not yet stable.  It includes
many draft algorithms that are sometimes posted on github in debug mode, or
with known suboptimal performance.  We ask that you not benchmark LAGraph on
your own without contacting the authors to make sure you have the right
version, and the right version of SuiteSparse:GraphBLAS to go with it.

However, assuming things are stable, the following instructions should work
for the GAP benchmarks.

Edit your .bash_profile, assuming you put GraphBLAS and LAGraph in $HOME:

    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/bin:$HOME/GraphBLAS/build:$HOME/LAGraph/build
    export LD_LIBRARY_PATH

Download GraphBLAS and LAGraph to $HOME, side-by-side (not one inside the
other).  For the latest bleeding edge versions:

    git clone -b cuda_merge_in_progress https://github.com/DrTimothyAldenDavis/GraphBLAS
    git clone https://github.com/GraphBLAS/LAGraph

But better yet, get a tagged version of each package, with the same date.
In particular, the versions of both codes with "11Dec2020" in their tag name
can be used together to give good results on BFS, BC, CC, and PR, but not
yet SSSP and TC.

You should have two folders, ~/GraphBLAS and ~/LAGraph, side-by-side.

For GraphBLAS, sometimes the latest version on the non-default MASTER
branch is non-functional; you may want to use a tagged pre-release instead:
See https://github.com/DrTimothyAldenDavis/GraphBLAS/releases .
So please ask me first (email: davis@tamu.edu), just in case.

--------------------------------------------------------------------------------
To compile and run interactively, assuming 40 threads to compile:
--------------------------------------------------------------------------------

    cd GraphBLAS
    make CC=gcc CXX=gcc JOBS=40
    cd ../LAGraph
    make
    cd Test
    make

Then, edit LAGraph/Test/BinRead/do_gap to fix the filnames, and do this once
to create the binary *.grb files:

    cd LAGraph/BinRead
    ./do_gap

Then to run the GAP benchmark do:

    cd LAGraph/Test
    ./do_gap_all > lengthy_output.txt

and watch stderr for a single line result for each benchmark.  Each line will
report method, the # of threads used ("40:" for example), the average run time,
and the matrix name, in that order.  To edit the # of threads used:

    cd LAGraph/Test
    vi */gap*.c

and edit the "#define NTHREAD_LIST" and THREAD_LIST settings.  If THREAD_LIST
is zero, the # of threads is chosen automatically, as the max number from
omp_get_num_threads, and then cutting in half for NTHREAD_LIST iterations.

--------------------------------------------------------------------------------
To compile and run on the Intel DevCloud:
--------------------------------------------------------------------------------

To compile, submit this to compile GraphBLAS and LAGraph (see contents below):

    qsub myjob_compile_plat_gcc

myjob_compile_plat_gcc contains the following:

    ----------------------------------------------------------------------
    #PBS -l nodes=1:ppn=2:plat8153
    #PBS -l walltime=24:00:00
    cd $PBS_O_WORKDIR
    cd GraphBLAS
    make CC=gcc CXX=gcc JOBS=24
    cd ../LAGraph
    make CC=gcc CXX=gcc JOBS=24
    cd Test/BinRead
    make CC=gcc CXX=gcc JOBS=24
    cd ../BetweennessCentrality
    make CC=gcc CXX=gcc JOBS=24
    cd ../BFS
    make CC=gcc CXX=gcc JOBS=24
    cd ../CC
    make CC=gcc CXX=gcc JOBS=24
    cd ../PageRank3
    make CC=gcc CXX=gcc JOBS=24
    cd ../TriangleCount
    make CC=gcc CXX=gcc JOBS=24
    cd ../SSSP
    make CC=gcc CXX=gcc JOBS=24
    ----------------------------------------------------------------------

Note that will break if you do not have the m4 command.

DO NOT USE icc; IT WILL GIVE POOR PERFORMANCE, because of how it does atomics.
I'm exploring workarounds.

You then need to convert the *mtx files to GraphBLAS binary *grb files for
faster I/O times.  My scripts assume the matrices are in ~/GAP/GAP-xxx/*.*
where xxx is kron, road, twitter, urand, or web.  Cut-and-paste the following,
just once:

    cd
    mkdir GAP
    mkdir GAP/GAP-kron
    mkdir GAP/GAP-urand
    mkdir GAP/GAP-twitter
    mkdir GAP/GAP-web
    mkdir GAP/GAP-road
    cp /data/sparse.tamu.edu/gap/GAP-kron/GAP-kron_sources.mtx ~/GAP/GAP-kron
    cp /data/sparse.tamu.edu/gap/GAP-road/GAP-road_sources.mtx ~/GAP/GAP-road
    cp /data/sparse.tamu.edu/gap/GAP-twitter/GAP-twitter_sources.mtx ~/GAP/GAP-twitter
    cp /data/sparse.tamu.edu/gap/GAP-urand/GAP-urand_sources.mtx ~/GAP/GAP-web
    cp /data/sparse.tamu.edu/gap/GAP-web/GAP-web_sources.mtx ~/GAP/GAP-road
    cd LAGraph/Test/BinRead
    qsub myjob_BinRead_gold

When that finishes, run the tests.  All the qsub scripts are already in the
LAGraph github.  You should be able to cut-and-paste these commands to run all
the scripts.  For the output, see the Test/*/myjob*.o* files that the scripts
generate.

    cd LAGraph/Test
    cd BetweennessCentrality
    qsub my1_BC_plat
    qsub my2_BC_plat
    cd ../BFS
    qsub my1_BFS_plat
    qsub my2_BFS_plat
    cd ../CC
    qsub my1_CC_plat
    qsub my2_CC_plat
    cd ../PageRank3
    qsub my1_PR_plat
    qsub my2_PR_plat
    cd ../SSSP
    qsub my1_SSSP_plat
    qsub my2_SSSP_plat
    cd ../TriangleCount
    qsub my1_TC_plat
    qsub my2_TC_plat

Those run the 2 sets: "default" (set 1, the my1* jobs), and "best effort" (the
my2* jobs).  You can also try the */my3_*_plat jobs, which use all default
settings (including # of threads, with hyperthreading enabled by default).
This will use whatever threads you have as the max, then it runs again with 1/2
that number.

Each test result will be placed in a file in each
test folder, as it completes, with the name log_HOST.txt, where
HOST is the hostname of the system that was used.

