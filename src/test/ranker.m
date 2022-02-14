% using MATLAB graph/centrality to compute the pagerank

clear all
version
s = dir ('../../data/*.mtx') ;

for k = 1:length(s)
    filename = ['../../data/' s(k).name] ;
    % fprintf ('\n=================== filename: %s\n', filename) ;
    try
        [A Z] = mread (filename) ;
        ok = true ;
    catch me
        ok = false ;
    end

    if (~ok)
        % fprintf ('BAD: %s bad\n', filename) ;
        continue
    end

    [m n] = size (A) ;
    if (m ~= n)
        % fprintf ('%s: rectangular\n', filename) ;
        continue ;
    end

    A = spones (spones (A) + spones (Z)) ;
    % get the out degree
    d = sum (A, 2) ;
    nsnks = length (find (d == 0)) ;
    if (nsnks > 0)
        fprintf ('%s has %d sinks\n', filename, nsnks) ;
        % continue ;
    end

    % make the graph
    G = digraph (A) ;
    pr = centrality (G, 'pagerank') ;

    fprintf ('OK: %g %6d %s\n', sum (pr), n, filename) ;

end

files = { 'karate.mtx', 'west0067.mtx', 'ldbc-directed-example.mtx' } ;

for k = 1:length (files)
    [A Z] = mread (['../../data/' files{k}]) ;
    A = spones (spones (A) + spones (Z)) ;
    G = digraph (A) ;
    pr = centrality (G, 'pagerank') ;
    fprintf ('\n%s\n', files {k}) ;
    fprintf ('%16.10f\n', pr) ;
end

