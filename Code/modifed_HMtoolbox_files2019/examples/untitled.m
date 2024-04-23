%NUDFT 

%start by testing the diagonalization: 
n = 100; 
F = zeros(n,n); 
w = exp(pi*1i/n); 
for j = 1:n
    for k = 1:n
        F(j,k) = w.^(-2*(j-1)*(k-1)); 
    end
end
F = F./sqrt(n); 

Q = eye(n); 
Q = circshift(Q,1); 

D = F'*Q*F;
D(abs(D) < 1e-14) = 0;

d_tru = (diag(D)).';

d = w.^((0:2:2*(n-1))); 

max(abs(d.' - d_tru.'))


%%
% write transformation using FFT: 








 

