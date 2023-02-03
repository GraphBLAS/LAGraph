import os

NUM_RUNS = 5

run_graphblas_cmd = '../build/experimental/benchmark/matching_demo < data.mtx > grb_result.txt'
run_verification_cmd = './verify < data.mtx > verif_result'

tests = [
    {
        'type': 'bipartite',
        'args': '1000 100 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive',
        'weighted': False
    },
    {
        'type': 'bipartite',
        'args': '1000 500 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive',
        'weighted': False
    },
    {
        'type': 'bipartite',
        'args': '1000 1000 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive',
        'weighted': False
    },
    {
        'type': 'general',
        'args': '1000 1000 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive',
        'weighted': False
    },
    {
        'type': 'general',
        'args': '1000 500 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive',
        'weighted': False
    },
    {
        'type': 'general',
        'args': '1000 50 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive',
        'weighted': False
    },
    {
        'type': 'general',
        'args': '1000 500 1 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive,weighted',
        'weighted': False
    },
    {
        'type': 'general',
        'args': '1000 100 1 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive,weighted',
        'weighted': False
    },
]

def build_cmd(type, args):
    return f'./gen_{type} {args} > my_result.txt'

def do_run(test):
    run_checker_cmd = build_cmd(test.get('type'), test.get('args'))

    max_matching = os.system(run_checker_cmd)
    os.system(run_graphblas_cmd)

    grb_result_file = open('grb_result.txt')
    my_result_file = open('my_result.txt')

    my_result = -1
    grb_result = -1

    for ln in grb_result_file:
        if('[DBG]' in ln):
            print(f'Debug output! : {ln}')
            exit()

    os.system(run_verification_cmd)
    verif_result_file = open('verif_result')

    for ln in verif_result_file:
        toks = list(map(int, ln.split(' ')))
        if(toks[0] == -1):
            type = test.get('type')
            args = test.get('args')
            print(f'Verification failed for {type} run w/ args {args}')
            exit()
        else:
            if(test.get('weighted')):
                grb_result = toks[1]
            else:
                grb_result = toks[0]

    os.remove('verif_result')

    for ln in my_result_file:
        my_result = int(ln)

    return [my_result, grb_result]


# sizes = [1000]#[100, 100, 100, 1000, 1000, 1000, 5000, 5000, 5000]
# sparse_factors = [500]#[50, 100, 200, 500, 1000, 2000, 2500, 5000, 10000]

for i in range(len(tests)):
    avg_ratio = 0
    runs = 0
    for _ in range(NUM_RUNS):
        run_result = do_run(tests[i])
        if (run_result[0] == 0):
            # empty graph
            continue
        avg_ratio += run_result[1] / run_result[0]
        runs += 1
    
    avg_ratio /= runs
    arg_list = tests[i].get('args').split(' ')
    arg_names_list = tests[i].get('arg_names').split(',')
    type = tests[i].get('type')
    result = f'(type: {type}, '
    for i in range(len(arg_list)):
        result += f'{arg_names_list[i]}: {arg_list[i]}'
        if(i < len(arg_list) - 1):
            result += ', '
    result += f'): performance ratio: {avg_ratio}'
    print(result)
