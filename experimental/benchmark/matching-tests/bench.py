'''
bench.py

Script that auto-runs a sequence of tests

Usage:
make
python3 bench.py

Details:

Test structure:
type [string]: <bipartite/general> (what type of graph)
performance [boolean]: <True/False> (whether to run a performance test)
args [string]: space-separated arguments for gen_bipartite/gen_general (both use same arguments, depends on which type is specified)
grb_args [string]: space-separated arguments for matching_demo (further specified in matching_demo.c)
arg_names [string]: comma-separated argument names to display for each argument in args
grb_arg_names [string]: comma-separated argument names to display for each argument in grb_args

To add new tests, simply use the following structure:
{
    'type': <type>
    'performance': <performance>
    'args': <args>
    'grb_args': <grb_args>
    'arg_names': <arg_names>
    'grb_arg_names': <grb_arg_names>
}
and add it to the tests list
'''

import os

# How many runs to do per test (results are averaged)
NUM_RUNS = 1

run_verification_cmd = './build/verify_matching < data.mtx > verif_result'


tests = [
    {
        'type': 'bipartite',
        'performance': True,
        'args': '1000000 10 1 1 0',
        'grb_args': 'data.mtx 0',
        'islight': False,
        'arg_names': 'num_nodes,sparse_factor,is_naive,perf,weighted,prefer_light',
        'grb_arg_names': 'filename,match_type,ntrials'
    },
]

# How many tests to run
NUM_TESTS = 1

'''
    {
        'type': 'bipartite',
        'performance': False,
        'args': '1000 500 0',
        'grb_args': 'stdin 0 5',
        'arg_names': 'num_nodes,sparse_factor,is_naive',
        'grb_arg_names': 'filename,match_type,ntrials',
        'weighted': False
    },
    {
        'type': 'bipartite',
        'performance': False,
        'args': '1000 1000 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive',
        'weighted': False
    },
    {
        'type': 'general',
        'performance': False,
        'args': '1000 1000 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive',
        'weighted': False
    },
    {
        'type': 'general',
        'performance': False,
        'args': '1000 500 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive',
        'weighted': False
    },
    {
        'type': 'general',
        'performance': False,
        'args': '1000 50 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive',
        'weighted': False
    },
    {
        'type': 'general',
        'performance': False,
        'args': '1000 500 1 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive,weighted',
        'weighted': False
    },
    {
        'type': 'general',
        'performance': False,
        'args': '1000 100 1 0',
        'arg_names': 'num_nodes,sparse_factor,is_naive,weighted',
        'weighted': False
    },
]
'''

def build_grb_cmd(args):
    return f'../../../build/experimental/benchmark/matching_demo {args} < data.mtx > grb_result.txt'


def build_cmd(type, args):
    return f'./build/gen_{type} {args} > my_result.txt'

def do_perf_run(test):
    run_checker_cmd = build_cmd(test.get('type'), test.get('args'))
    run_graphblas_cmd = build_grb_cmd(test.get('grb_args'))

    os.system(run_checker_cmd)
    # print('done w/ checker')
    os.system(run_graphblas_cmd)

    grb_result_file = open('grb_result.txt')
    my_result_file = open('my_result.txt')

    last_ln = ''
    for ln in grb_result_file:
        last_ln = ln

    my_time = -1
    grb_time = float(last_ln)

    for ln in my_result_file:
        my_time = float(ln)

    return [grb_time, my_time] 

def do_run(test):
    run_checker_cmd = build_cmd(test.get('type'), test.get('args'))
    run_graphblas_cmd = build_grb_cmd(test.get('grb_args'))

    os.system(run_checker_cmd)
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
            grb_result = toks[0]

    os.remove('verif_result')

    for ln in my_result_file:
        my_result = int(ln)

    return [my_result, grb_result]


for i in range(NUM_TESTS):
    avg_ratio = 0
    avg_perf_ratio = 0
    runs = 0
    perf_runs = 0
    for _ in range(NUM_RUNS):
        if(tests[i].get('performance') == True):
            # doing a perf test
            # print('passed')
            run_result = do_perf_run(tests[i])
            avg_perf_ratio += run_result[1] / run_result[0]
            perf_runs += 1
        else:
            run_result = do_run(tests[i])
            if (run_result[0] == 0):
                # empty graph
                continue
            if(tests[i].get('islight') == True):
                avg_ratio += run_result[0] / run_result[1]
            else:
                avg_ratio += run_result[1] / run_result[0]
            runs += 1
    
    if (perf_runs > 0):
        avg_perf_ratio /= perf_runs
    else:
        avg_perf_ratio = float('nan')
    
    if (runs > 0):
        avg_ratio /= runs
    else:
        avg_ratio = float('nan')

    arg_list = tests[i].get('args').split(' ')
    arg_names_list = tests[i].get('arg_names').split(',')
    grb_arg_list = tests[i].get('grb_args').split(' ')
    grb_arg_names_list = tests[i].get('grb_arg_names').split(',')
    type = tests[i].get('type')
    result = f'(type: {type}, MY ARGS: '
    for i in range(len(arg_list)):
        result += f'{arg_names_list[i]}: {arg_list[i]}'
        if(i < len(arg_list) - 1):
            result += ', '
    result += ' GRB ARGS: '
    for i in range(len(grb_arg_list)):
        result += f'{grb_arg_names_list[i]}: {grb_arg_list[i]}'
        if(i < len(grb_arg_list) - 1):
            result += ', '
    
    result += f'): quality metric: {avg_ratio}, performance metric: {avg_perf_ratio}'
    print(result)

print('All tests passed')