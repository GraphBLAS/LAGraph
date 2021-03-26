
clear
rank3c = mread ('rank3c.mtx') ;
rank3d = mread ('rank3d.mtx') ;
rank3e = mread ('rank3e.mtx') ;

e_cd = norm (rank3c - rank3d, 1)
e_de = norm (rank3d - rank3e, 1)
e_ce = norm (rank3c - rank3e, 1)
