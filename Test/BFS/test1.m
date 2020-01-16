
clear all
system ('./build/gap_bfs_test < cover.mtx') ;

% source nodes are 0-based in s_7.mtx

A  = mread ('cover.mtx')
s  = mread ('s_7.mtx') + 1
v  = mread ('v_7.mtx')
pi = mread ('p_7.mtx')

s = full (s (end))

bfs_check (A, s, v, pi)

