function bc_test (A, filename);

if (nargin < 2)
    filename = 'A.mtx';
end

fprintf('===== Betweenness Centrality Test - Comparing with MATLAB =====\n\n');

fprintf('Computing Betweenness Centrality with MATLAB...\n');

G = digraph(A);

% test the MATLAB version
n = size(A, 1);
v = centrality(G, 'betweenness');

% normalize
v_min = min(v);
v_max = max(v);
v = (v - v_min) / (v_max - v_min);

fprintf('# of nodes in graph: %d\n', n);

%fprintf('Betweenness Centrality:\n');
%for i = 1:length(v)
%	fprintf("v(%d) = %f\n", i, v(i));
%end

% Compare the Brandes algorithm with MATLAB
fprintf ('\nTesting LAGraph_bc:\n');
outfile = 'v_brandes';
mwrite (filename, sparse (A));

system (sprintf ('./build/bc_exe < %s > %s', filename, outfile)); 
v_brandes = load (outfile);

% normalize
v_min = min(v_brandes);
v_max = max(v_brandes);
v_brandes = (v_brandes - v_min) / (v_max - v_min);

assert (max(abs(v - v_brandes)) < 0.1);
fprintf('\nBetweenness Centrality (Brandes) - Test Passed!\n')

% Compare the batch algorithm with MATLAB
fprintf ('\nTesting LAGraph_bc_batch:\n');
outfile = 'v_batch';

system (sprintf ('./build/bc_batch_exe < %s > %s', filename, outfile)); 
v_batch = load (outfile);

% normalize
v_min = min(v_batch);
v_max = max(v_batch);
v_batch = (v_batch - v_min) / (v_max - v_min);

assert (max(abs(v - v_batch)) < 0.1);
fprintf('\nBetweenness Centrality (Batch) - Test Passed!\n');

fprintf ('bc_test: all tests passed\n') ;

