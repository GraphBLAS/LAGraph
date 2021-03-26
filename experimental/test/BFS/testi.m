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

Prob = ssget ('GAP/GAP-twitter')
A = spones (Prob.A) ;
pi = mread ('p_61578415.mtx') ;
s = mread ('s_61578415.mtx') ;   % already 1-based
v = mread ('v_61578415.mtx') ;
s = full (s (end))
A = GrB (A)
bfs_check (A, s, v, pi)
%}
pi2 = mread ('ponly_61578415.mtx') ;
bfs_check (A, s, v, pi2)
clear A 
clear Prob
clear pi s v a

% system ('./build/gap_bfs_test ~/com-Orkut.grb > o') ;
Prob = ssget ('SNAP/com-Orkut')
A = spones (Prob.A) ;
pi = mread ('p_3072441.mtx') ;
s = mread ('s_3072441.mtx') ;   % already 1-based
v = mread ('v_3072441.mtx') ;
s = full (s (end))
A = GrB (A)
bfs_check (A, s, v, pi)
%}
pi2 = mread ('ponly_3072441.mtx') ;
bfs_check (A, s, v, pi2)

