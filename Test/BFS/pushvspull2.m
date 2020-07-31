function pushvspull (name, allpush_results, allpull_results, n, nz, fig)
% pushvspull ('kron', allpush_results, allpull_results, n, nz, 1) ;

d = nz / n 
ntrials = length (allpush_results) ;
figure (fig)
hold on

for k = 1:ntrials
    
    push_result = allpush_results {k} ;
    pull_result = allpull_results {k} ;

    nq = push_result (:,1) ;
    nvisited = push_result (:,2) ;
    t_pull = pull_result (:,4) ;
    t_push = push_result (:,4) ;

    pushwork = d * nq ;
    % pushwork (f) = pushwork (f) - 4 *nvisited (f) ;

    expected = n ./ (nvisited+1) ;
    per_dot = min (d, expected) ;
    binarysearch = (3 * (1 + log2 (nq))) ;
    pullwork = (n-nvisited) .* per_dot .* binarysearch ;

    relwork = pushwork ./ pullwork ;
    reltime = t_push ./ t_pull ;

    % just = find ((nvisited-nq) < nvisited/32 & (relwork > 0.01 ))  ;
    % better = find (reltime > 1) ;

    predict = find (pullwork < pushwork)
    tt = [t_pull t_push] ;
    here= tt (predict,:)
    sum (here)

    subplot (2,2,1)
    loglog ( relwork (:), reltime (:), 'o', ...
             [1e-8 100], [1 1], 'g-', ...
             [1 1], [1e-3 10], 'g-') ;

    subplot (2,2,2)
    semilogy ( nvisited/n, reltime, 'o') ;

    subplot (2,2,3)
    x = (nvisited -nq) / n; % ./(nvisited-nq) ;
    semilogy ( x, t_push, 'o', x, t_pull, 'rx') ;

    subplot (2,2,4)
    loglog ( ...
        relwork (:), t_push (:), 'o', ...
        relwork (:), t_pull (:), 'rx', ...
             [1 1], [1e-3 10], 'g-') ;

end
title (name) ;
hold off
