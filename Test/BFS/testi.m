clear all


% source nodes are 0-based

%{
% system ('./build/gap_bfs_test ~/indochina-2004.grb') ;
Prob = ssget ('LAW/indochina-2004')
A = spones (Prob.A) ;
pi = mread ('p_7414866.mtx') ;
s = mread ('s_7414866.mtx') + 1 ;
v = mread ('v_7414866.mtx') ;
s = full (s (end))
bfs_check (A, s, v, pi)
%}

%{
Prob = ssget ('SNAP/roadNet-CA')
A = spones (Prob.A) ;
pi = mread ('p_1971281.mtx') ;
s = mread ('s_1971281.mtx') + 1;
v = mread ('v_1971281.mtx') ;
s = full (s (end))
bfs_check (A, s, v, pi)
%}

% system ('./build/gap_bfs_test ~/com-Orkut.grb > o') ;
Prob = ssget ('SNAP/com-Orkut')
A = spones (Prob.A) ;
pi = mread ('p_3072441.mtx') ;
s = mread ('s_3072441.mtx') + 1 ;
v = mread ('v_3072441.mtx') ;
s = full (s (end))
bfs_check (A, s, v, pi)
%}
