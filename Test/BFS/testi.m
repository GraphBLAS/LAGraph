clear all

system ('./build/gap_bfs_test ~/indochina-2004.grb') ;

% source nodes are 0-based

Prob = ssget ('LAW/indochina-2004')
A = spones (Prob.A) ;
pi = mread ('p_7414866.mtx') ;
s = mread ('s_7414866.mtx') + 1 ;
v = mread ('v_7414866.mtx') ;

s = full (s (end))

bfs_check (A, s, v, pi)

%{
Prob = ssget ('SNAP/roadNet-CA')
A = spones (Prob.A) ;
pi = mread ('p_1971281.mtx') ;
s = mread ('s_1971281.mtx') + 1;
v = mread ('v_1971281.mtx') ;
s = full (s (end))
bfs_check (A, s, v, pi)
%}
