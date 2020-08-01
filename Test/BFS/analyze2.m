
% simple:
% [allpull_results, name, n, nz] = allpull_7 ;
% [allpush_results, name, n, nz] = allpush_7 ;
% pushvspull2 ('simple', allpush_results, allpull_results, n, nz, 1)

% orkut:
[pushpull_results, name, n, nz] = pushpull_3072441_slash ;
[ allpull_results, name, n, nz] = allpull_3072441_slash ;
fprintf ('\n\n=================== %s\n', name) ;
[ allpush_results, name, n, nz] = allpush_3072441_slash ;
pushvspull2 ('orkut (slash)', allpush_results, allpull_results, ...
    pushpull_results, n, nz, 1)

% orkut:
[pushpull_results, name, n, nz] = pushpull_3072441 ;
[ allpull_results, name, n, nz] = allpull_3072441 ;
fprintf ('\n\n=================== %s\n', name) ;
[ allpush_results, name, n, nz] = allpush_3072441 ;
pushvspull2 ('orkut (hypersparse)', allpush_results, allpull_results, ...
    pushpull_results, n, nz, 2)

% kron
[pushpull_results, name, n, nz] = pushpull_134217726 ;
[ allpull_results, name, n, nz] = allpull_134217726 ;
fprintf ('\n\n=================== %s\n', name) ;
[ allpush_results, name, n, nz] = allpush_134217726 ;
pushvspull2 ('kron', allpush_results, allpull_results, pushpull_results, ...
    n, nz, 3)

% urand
[pushpull_results, name, n, nz] = pushpull_134217728 ;
[ allpull_results, name, n, nz] = allpull_134217728 ;
fprintf ('\n\n=================== %s\n', name) ;
[ allpush_results, name, n, nz] = allpush_134217728 ;
pushvspull2 ('urand', allpush_results, allpull_results, pushpull_results, ...
    n, nz, 4)

