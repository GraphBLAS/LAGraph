function pushvspull (name, results, n, nz, fig)
% pushvspull ('kron', results, n, nz, 1) ;

d = nz / n 
nq = results (:,1) ;
nvisited = results (:,2) ;
t_pull = results (:,3) ;
t_push = results (:,4) ;

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

figure (fig)
subplot (2,2,1)
loglog ( relwork (:), reltime (:), 'o', ...
         [1e-8 100], [1 1], 'g-', ...
         [1 1], [1e-3 10], 'g-') ;
title (name) ;

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
title (name) ;
title (name) ;

