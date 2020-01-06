
% clear all

% brand40b = mread ('brandes_result.mtx') ;
% batch40b = mread ('batch_result.mtx') ;
difb = abs (brand40b - batch40b) ;
errb = dif / norm (brand40b) ;
badb = find (err > 1e-4)
full ([brand40b(badb) batch40b(badb) difb(badb)]) ;

look = [
    14267380
    38003612
    39861733
    40381015
    52064540
    53018295 ] ;

e = abs (batch1-batch40b) ; e (look)
e = abs (brand1-batch40b) ; e (look)

e = abs (brand1-batch1) ; e (look)
