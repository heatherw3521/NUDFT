function y = ct_mul(t0,t1, x)
% fast multiplication Cx = y,
% where C is the cauchy-toeplitz matrix, 
% t0 = toeplitz row entries. 
% t1 = toeplitz col entries. 


t0 = t0(:);
t1 = t1(:);
n = length(t0); 
w = exp(pi*1i/n); 
%D0 = diag(w.^(-(0:n-1))); 
D0 = spdiags((w.^(-(0:n-1))).', 0, n, n);  

%step 1: transform x
x = D0*(fft(x)/sqrt(n)); 

%step 2: T*x by embedded circulant
t = [t1; 0; flip(t0)];
t = t(1:end-1); 
x = [x; zeros(size(x))]; 
 
x = fft(x); 
v = fft(t); 
v = spdiags(v, 0, 2*n, 2*n);
y = ifft(v*x);

%step 3: C = F(Tx)
y = y(1:n, :); %T*x
y = sqrt(n)*ifft(y); %Cx = F*T*x

end
