function C = build_cauchys(n)
%
% builds a (singular) Cauchy matrix with entries $C_{jk} = 1./(w^(2*(j-1)) - w^2(k)-1,
% where w = exp(i pi/n). 
% this is related to the Fourier transform of a Toeplitz matrix. 
% diag of this matrix is undefined

w = exp(1i*pi/n); 
a = w.^(2*(1:n)); 

C = bsxfun(@minus, a.', a); 
C = 1./C; 
end