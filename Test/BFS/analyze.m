% simple
%{
n = 7 ;
nz = 12 ;
results = load (sprintf ('logtime_%d_%d.txt', n, 40)) ;
pushvspull ('simple', results, n, nz, 5) ;
%}

% kron:
n = 134217726 ;
nz = 2111632322 ;
results = load (sprintf ('logtime_%d_%d.txt', n, 40)) ;
pushvspull ('kron:40', results, n, nz, 1) ;
% logtime_134217726_40.txt

results = load (sprintf ('logtime_%d_%d.txt', n, 20)) ;
pushvspull ('kron:20', results, n, nz, 2) ;

%{
% road:
% n = 23947347 ;
% nz = 28854312 ;
% d = nz / n ;
%}

% urand:
n = 134217728 ;
nz = 2147483370 ;
results = load (sprintf ('logtime_%d_%d.txt', n, 40)) ;
pushvspull ('urand', results, n, nz, 3) ;

%{
% twitter
n = 61578415 ; 
nz = 1468364884 ;
results = load (sprintf ('logtime_%d_%d.txt', n, 40)) ;
pushvspull ('urand', results, n, nz, 3) ;

% web:
n = 50636151 ;
nz = 1930292948 ;
results = load (sprintf ('logtime_%d_%d.txt', n, 40)) ;
pushvspull ('urand', results, n, nz, 4) ;

%}
