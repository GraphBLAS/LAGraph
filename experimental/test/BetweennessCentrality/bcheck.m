
clear all

% 40 thread result
brand40 = mread ('40/brandes_result.mtx') ;
batch40 = mread ('40/batch_result.mtx') ;
dif = abs (brand40 - batch40) ;
err = dif / norm (brand40) ;
bad = find (err > 1e-4)
full ([brand40(bad) batch40(bad) dif(bad)]) ;

% 1 thread result
brand1 = mread ('1/brandes_result.mtx') ;
batch1 = mread ('1/batch_result.mtx') ;
dif1 = abs (brand1 - batch1) ;
err1 = dif / norm (brand1) ;
bad1 = find (err1 > 1e-4)
full ([brand1(bad1) batch1(bad1) dif(bad1)]) ;

look = [
    14267380
    38003612
    39861733
    40381015
    52064540
    53018295 ] ;


