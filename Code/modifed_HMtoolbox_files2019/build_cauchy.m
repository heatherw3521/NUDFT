function C = build_cauchy(n)
%
% builds a Cauchy matrix with entries $C_{jk} = 1./(w^(2*(j-1)) - w^2(k)-1,
% where w = exp(i pi/n). 
% this is related to the Fourier transform of a Toeplitz matrix. 
% see Memo 19 for details

w = exp(1i*pi/n); 
a = w.^(2*(0:n-1)); 

C = bsxfun(@minus, a.', w*a); 
C = 1./C; 
end


