
% Test lcc.m.  Requires ssget.

clear
index = ssget ;
f = find (index.nrows == index.ncols) ;
[~,i] = sort (index.nnz (f)) ;
f = f (i) ;

% f = [739 750 2662];
% f = [168 739 750] ; 
% f = [168] 
% f = [ 262 2459 ] ;

%   matrices = {
%   'HB/west0067'
%   'DIMACS10/road_usa'
%   'SNAP/roadNet-CA'
%   'McRae/ecology1'
%   'ND/nd3k'
%   'ND/nd12k'
%   'SNAP/soc-LiveJournal1'
%   'LAW/hollywood-2009'
%   'LAW/indochina-2004'
%   'DIMACS10/kron_g500-logn21' } ;

matrices = f ;
nmat = length (matrices)

for kk = 1:nmat

    % get a matrix
    matrix = matrices {kk} ;
    Prob = ssget (matrix, index) ;
    id = Prob.id ;
    A = Prob.A ;
    if (isfield (Prob, 'Zeros')) 
        A = A + Prob.Zeros ;
    end
    fprintf ('\n%4d %4d %-40s\n', kk, id, Prob.name) ;
    clear Prob 

    % compute the result in MATLAB
    c1 = lcc (A) ;

    % compute the result in GraphBLAS
    mwrite (A, 'A.mtx') ;
    system ('./lcc < A.mtx > lcc_results') ;
    c2 = load ('lcc_results') ;

    err = norm (c1-c2,1)
    if (err > 0)
        pause
    end
end

