
[a2 m2 s2 q2] = bfs2_n36000 ;
[ap mp sp  q] = bfs_pushpull_n36000 ;
figure (2) ; title ('36000') ;
doplots (a2,m2,s2,q2, ap,mp,sp,q, 1.42e7)

[ap mp sp q2] = bfs_pushpull_n9000 ;
[a2 m2 s2  q] = bfs2_n9000 ;
figure (3) ; title ('9000') ;
doplots (a2,m2,s2,q2, ap,mp,sp,q, 3.3e6)


[a2 m2 s2 q2] = bfs2_n1000000 ;
[ap mp sp  q] = bfs_pushpull_n1000000 ;
figure (1) ; title ('1000000') ;
doplots (a2,m2,s2,q2, ap,mp,sp,q, 5e6)
