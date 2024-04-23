%check identities: 

%check the displacement identities in paper:

%set up parameters defining the Toeplitz matrix: 
% and original displacement equation
n = 10; 
t1 = rand(n, 1); 
t2 = rand(n,1); 
t2(1) = t1(1); 
T = toeplitz(t1, t2); 
Z = circshift(eye(n), 1); 
Zm = Z; Zm(1,end) = -1; 
%%
% 
% check transforms on Z and Zm: 
w = exp(1i*pi/n); 
D1t = w.^(2*(0:n-1)); 
Dmt = w.^(2*(1:n)-1); 
D0 = diag(w.^(0:n-1)); 

%transform on Z: 
D1 = fft((ifft(Z)).').'; 
norm(diag(D1) - D1t.')

%transform on Zm: 
Dm = D0*Zm*D0'; 
Dm = fft((ifft(Dm)).').'; 
norm(diag(Dm) - Dmt.')

 
[H, G, F] = toeplitz2cauchy(t1, t2); 

%% 
% bb version: 
F = zeros(n, n); 
for j = 1:n
    for k = 1:n
        F(j,k) = w^(j*(2*k-1)); 
    end
end
F = F/sqrt(n); 
L = w.^((2:2:2*n)); 

L2 = F*Z*F';  
L2 = L2(:); 
L2(abs(L2) < 1e-5)=0; 
L2 = reshape(L2, n,n); 
spy(L2)
Lt = diag(L2).'; 

Finv = zeros(n, n); 
for j = 1:n
    for k = 1:n
        F(j,k) = w^(-k*(2*j-1)); 
    end
end

Finv = sqrt(n)*Finv; 
