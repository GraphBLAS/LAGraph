import os

NUM_RUNS = 100

run_graphblas_cmd = '../build/experimental/benchmark/matching_demo < data.mtx > grb_result.txt'

def do_run(size, sparse):
    run_checker_cmd = f'./gen_bipartite {size} {sparse} > my_result.txt'
    max_matching = os.system(run_checker_cmd)
    os.system(run_graphblas_cmd)

    grb_result_file = open('grb_result.txt')
    my_result_file = open('my_result.txt')

    my_result = -1
    grb_result = -1

    for ln in my_result_file:
        my_result = int(ln)

    seenEntries = False

    for ln in grb_result_file:
        if('[DBG]' in ln):
            print(f'Debug output! : {ln}')
            exit()
        if(('entries: ' in ln) and (not seenEntries)):
            grb_result = int(ln[ln.index('entries') + 9:])
            seenEntries = True
    
    return [my_result, grb_result]

sizes = [1000]#[100, 100, 100, 1000, 1000, 1000, 5000, 5000, 5000]
sparse_factors = [500]#[50, 100, 200, 500, 1000, 2000, 2500, 5000, 10000]

for i in range(len(sizes)):
    avg_ratio = 0
    runs = 0
    for _ in range(NUM_RUNS):
        run_result = do_run(sizes[i], sparse_factors[i])
        if (run_result[0] == 0):
            # empty graph
            continue
        assert(run_result[1] <= run_result[0])
        avg_ratio += run_result[1] / run_result[0]
        runs += 1
    
    print(f'(n={sizes[i]}, k={sparse_factors[i]}), optimality ratio: {avg_ratio}')

