function doplots (a2,m2,s2,q2, ap,mp,sp,q)
k = length (m2) ;

% black is push
% red is pull
subplot (1,3,1) ; semilogy (1:k, a2, 'ko-', 1:k, ap, 'ro-') ;
subplot (1,3,2) ; semilogy (1:k, m2, 'ko-', 1:k, mp, 'ro-') ;
subplot (1,3,3) ; semilogy (1:k, s2, 'ko-', 1:k, sp, 'ro-') ;

assert (isequal (q,q2)) ; 

best = min (m2, mp) ;
fprintf ('\n') ;
fprintf ('push : %12.6e\n', sum (m2)) ;
fprintf ('pull : %12.6e\n', sum (mp)) ;
fprintf ('best : %12.6e\n', sum (best)) ;

