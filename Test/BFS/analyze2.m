
% simple:
% [allpull_results, name, n, nz] = allpull_7 ;
% [allpush_results, name, n, nz] = allpush_7 ;
% pushvspull2 ('simple', allpush_results, allpull_results, n, nz, 1)

% orkut:
[pushpull_results, name, n, nz] = pushpull_3072441 ;
[ allpull_results, name, n, nz] = allpull_3072441 ;
fprintf ('\n\n=================== %s\n', name) ;
[ allpush_results, name, n, nz] = allpush_3072441 ;
pushvspull2 ('orkut', allpush_results, allpull_results, pushpull_results, n, nz, 1)

%{
% kron
[pushpull_results, name, n, nz] = pushpull_134217726 ;
[ allpull_results, name, n, nz] = allpull_134217726 ;
fprintf ('\n\n=================== %s\n', name) ;
[ allpush_results, name, n, nz] = allpush_134217726 ;
pushvspull2 ('kron', allpush_results, allpull_results, pushpull_results, n, nz, 2)
%}
