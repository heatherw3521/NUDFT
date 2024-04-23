function AA = buildcauchy(J, K, N, nodes, G, L)
%builds the Cauchy matrix given by a transformed NUDFT matrix
%build the cauchy matrix from the given indices:
J = J(:); 
K = K(:); 
w = exp(1i*pi/N); 
a = nodes(J);
b = w.^(2*K); 
A = bsxfun(@minus, a, b.'); 
A = 1./A; 


%apply hadamard product
%[~,r] =size(G); 
k = length(J);
n = length(K); 
%AA = 0; 
%for j = 1:r
AA = spdiags(G(J), 0,k,k)*A*spdiags(conj(L(K)), 0, n,n);
    %AA = AA +B; 
%end

end