
clear all
matrices = [262 2852:2856]

for id = matrices
    Prob = ssget (id)
    A = GrB (Prob.A, 'logical') ;
    A = A | A' ;
    d = full (double (GrB.entries (A, 'col', 'degree'))) ;
    dfile = sprintf ('deg_%d', id) ;
    save (dfile, 'd') ;
    n = size (A,1)

    dmin = min (d)
    dmax = max (d)
    dmean = mean (d)
    dmedian = median (d)
    relative = dmean / dmedian
end
