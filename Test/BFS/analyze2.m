
% simple:

[allpull_results, name, n, nz] = allpull_7 ;
[allpush_results, name, n, nz] = allpush_7 ;
pushvspull2 ('simple', allpush_results, allpull_results, n, nz, 1)
