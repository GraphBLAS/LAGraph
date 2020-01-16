function bfs_check (A, s, v, pi)
% BFS_CHECK check the result of a breadth-first-search of a graph
%
% Given the adjacency matrix A of a graph, where A(i,j) is the edge (i,j),
% and a source node s, check the result v (and optionally pi) of the BFS.
% v is the level of the node, with v(s)=1 as the first level.
% pi, if present, is the BFS tree.  pi(s)=s, otherwise pi(i) is the parent
% of node i.
%
% See also GrB/bfs, graph/bfssearch, treeplot.

if (isequal (GrB.format (A), 'by col'))
    % ensure A is by-row
    A = spones (GrB (A, 'by row')) ;
else
    A = spones (GrB (A)) ;
end

A = spones (GrB (A)) ;
n = size (A, 1) ;

nsearched = GrB.entries (v) ;
fprintf ('source: %d\nn: %d\n# nodes found: %d\n', s, n, nsearched) ;
if (nargin > 3)
    % make sure each node searched has a parent
    assert (nsearched == length (find (pi > 0))) ;
end

nlevels = max (v) ;

for k = 1:n

    % get all nodes at the current level k
    q = find (v == k) ;

    if (k == 1)
        % level 1 is just node s itself
        assert (isequal (q, s)) ;
        if (nargin > 3)
            assert (pi (s) == s) ;
        end
    end

    if (length (q) == 0)
        fprintf ('# of levels: %d\n', k-1) ;
        assert (nlevels == k-1) ;
        break ;
    end

    % get all out-neighbors of the current level k
    S = A (q, :) ;
    [i j x] = find (S) ;

    % ensure all neighbors x are in level k+1 or less
    assert (all (v (j) <= k+1)) ;

    % ensure the parents of all nodes in q are in level k-1
    if (nargin > 3 && k > 1)
        parents = full (double (pi (q))) ;
        levels = v (parents) ;
        assert (all (levels == k-1)) ;
    end
end

if (nargin > 3 && n < 10000)
    pi = full (double (pi)) ;
    pi (s) = 0 ;
    treeplot (pi) ;
end

fprintf ('BFS ok\n') ;

