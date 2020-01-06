
% clear all

% 40 thread result
% brand40 = mread ('40/brandes_result.mtx') ;
% batch40 = mread ('40/batch_result.mtx') ;
dif = abs (brand40 - batch40) ;
err = dif / norm (brand40) ;
bad = find (err > 1e-4)
full ([brand40(bad) batch40(bad) dif(bad)])

% 1 thread result
% brand1 = mread ('1/brandes_result.mtx') ;
% batch1 = mread ('1/batch_result.mtx') ;
dif1 = abs (brand1 - batch1) ;
err1 = dif1 / norm (brand1) ;
bad1 = find (err1 > 1e-4)
full ([brand1(bad1) batch1(bad1) dif(bad1)])

look = [
    14267380
    38003612
    39861733
    40381015
    52064540
    53018295 ] ;

% brand40b = mread ('40b/brandes_result.mtx') ;
% batch40b = mread ('40b/batch_result.mtx') ;
difb = abs (brand40b - batch40b) ;
errb = difb / norm (brand40b) ;
badb = find (errb > 1e-4)
full ([brand40b(badb) batch40b(badb) difb(badb)])



% same as brand1
e = abs (brand1-batch1) ; e (look)
e = abs (brand1-brand1) ; e (look)
e = abs (brand1-brand40) ; e (look)
e = abs (brand1-brand40b) ; e (look)

% differ from brand1
e = abs (brand1-batch40) ; e (look)
e = abs (brand1-batch40b) ; e (look)

% but these match
e = abs (batch40-batch40b) ; e (look)

e = abs (brand40b-batch40b) ; e (look)

% fp32
% brand40c = mread ('40c/brandes_result.mtx') ;
% batch40c = mread ('40c/batch_result.mtx') ;

difc = abs (brand40c - batch40c) ;
errc = difc / norm (brand40c) ;
badc = find (errc > 1e-4)
full ([brand40c(badc) batch40c(badc) difc(badc)])
