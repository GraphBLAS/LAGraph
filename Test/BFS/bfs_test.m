
% graph on the cover of the book, 'Graph Algorithms in the language
% of linear algebra'.  The source node is node 4.
clear
o = 0 ;

A = [
o o o 1 o o o
1 o o o o o o
o o o 1 o 1 1
1 o o o o o 1
o 1 o o o o 1
o o 1 o 1 o o
o 1 o o o o o ] ;

A = A' ;

s = 4 ;
v = bfs_matlab (A, s) ;
n = size (A,1) ;

fprintf ('# of nodes in graph: %d\n', n) ;
fprintf ('source node: %d\n', s) ;
fprintf ('number of levels: %d\n', max (v)) ;
fprintf ('reachable nodes (incl. source): %d\n', length (find (v > 0))) ;

for level = 1:n
    q = find (v == level)' ;
    if (isempty (q))
        break ;
    end
    fprintf ('nodes in level %d: ', level) ;
    fprintf ('%d ', q) ;
    fprintf ('\n') ;
end

vok = [2 3 2 1 4 3 4]' ;
assert (isequal (v, vok)) ;

% now try LAGraph

% TODO: for now, this requires mread and mwrite from CHOLMOD:
infile = 'cover.mtx' ;
outfile = 'v' ;
mwrite (filename, sparse (A)) ;
% convert s to zero-based
system (sprintf ('./build/bfs_test %d < %s > %s', s-1, infile, outfile)) ; 
v2 = load (outfile) ;

assert (isequal (v, v2)) ;

fprintf ('bfs_test: all tests passed\n') ;
