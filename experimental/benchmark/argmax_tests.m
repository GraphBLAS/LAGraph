
% for HPEC'24 paper

if (0)
    clear all
    Prob = ssget ('GAP/GAP-twitter')
end
A = Prob.A ;
nz = nnz (A) 
n = size (A,1) ;
A = A + speye (n) ;
nz = nnz (A) 
G = GrB (A) ;

% time the GraphBLAS max and argmax methods
for thr = [1 40]
    GrB.threads (thr) ;
    for trial = 1:3
        fprintf ('\ntrial %d, threads %g\n', trial, thr) ;

        t = tic ;
        x = max (G, [ ], 1) ;
        t1 = toc (t) ;
        fprintf ('GrB colwise max: %g sec\n', t1) ;
        t = tic ;
        [x,p] = GrB.argmax (G, 1) ;
        t2 = toc (t) ;
        fprintf ('GrB colwise argmax: %g sec\n', t2) ;
        fprintf ('GrB colwise argmax time / max time: %g\n', t2/t1) ;

        t = tic ;
        x = max (G, [ ], 2) ;
        t1 = toc (t) ;
        fprintf ('GrB rowwise max: %g sec\n', t1) ;
        t = tic ;
        [x,p] = GrB.argmax (G, 2) ;
        t2 = toc (t) ;
        fprintf ('GrB rowwise argmax: %g sec\n', t2) ;
        fprintf ('GrB rowwise argmax time / max time: %g\n', t2/t1) ;

    end
end


% time the MATLAB max and argmax methods
for trial = 1:3
    fprintf ('\ntrial %d\n', trial) ;

    t = tic ;
    x = max (A, [ ], 1) ;
    t1 = toc (t) ;
    fprintf ('MATLAB colwise max: %g sec\n', t1) ;
    t = tic ;
    [x,p] = max (A, [ ], 1) ;
    t2 = toc (t) ;
    fprintf ('MATLAB colwise argmax: %g sec\n', t2) ;
    fprintf ('MATLAB colwise argmax time / max time: %g\n', t2/t1) ;

    t = tic ;
    x = max (A, [ ], 2) ;
    t1 = toc (t) ;
    fprintf ('MATLAB rowwise max: %g sec\n', t1) ;
    t = tic ;
    [x,p] = max (A, [ ], 2) ;
    t2 = toc (t) ;
    fprintf ('MATLAB rowwise argmax: %g sec\n', t2) ;
    fprintf ('MATLAB rowwise argmax time / max time: %g\n', t2/t1) ;

end

