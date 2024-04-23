function y = ct_mul_tp(t0,t1, x)
% fast multiplication (C^T)x = y,
% where C is the cauchy-toeplitz matrix, 
% t0 = toeplitz row entries. 
% t1 = toeplitz col entries. 


t0 = t0(:);
t1 = t1(:);
n = length(t0); 
w = exp(pi*1i/n); 
%D0 = diag(w.^(-(0:n-1))); 
%D0 = diag(w.^((0:n-1)));
D0 = spdiags( (w.^((0:n-1))).', 0, n,n); 

%step 1: transform x
%x = sqrt(n)*ifft(x);
x = fft(x)/sqrt(n); 

%step 2: T'*x by embedded circulant
t = [t0; 0; flip(t1)]; 
t = t(1:end-1); 
t = conj(t); 
x = [x; zeros(size(x))]; 
 
x = fft(x); 
v = fft(t); 
v = spdiags(v, 0, 2*n, 2*n);
y = ifft(v*x);

%step 3: C^TX = F *D_0*(T'*x)
y = y(1:n, :); %T*x
y = D0*y; 
y = sqrt(n)*ifft(y); 

end
