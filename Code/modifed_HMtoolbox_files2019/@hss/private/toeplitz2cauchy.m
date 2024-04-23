
function [H, G, F] = toeplitz2cauchy(t1, t2, varargin)
%given two vectors defining a toeplitz matrix, this converts
% to the HSS decomposition H of a Cauchy-like matrix using FFTs. 
% the transformed rank 2 matrices that make up the RHS of the 
% displacement eqn are also returned. 

%
n = length(t1); 
t1 = t1(:); 
t2 = t2(:); 
t0 = t1(1); 
if ~(t0 == t2(1))
    error('toeplitz2cauchy: first element of t1 does not match first element of t2.')
end

% build rank 2 matrix from displ struct, RHS is G*F'
t1 = t1(2:end); 
t2 = t2(2:end); 
e1 = [1; zeros(n-1, 1)];

G = [e1 [t0; flip(t2) + t1]]; 
F = [ [flip(t1)-t2; t0] flip(e1)]; 

% apply transformation: we use normalized fft:
D0 = exp(1i*pi/n).^(0:n); 

G = ifft(G)*sqrt(n); 
F = ifft(diag(D0)*F)*sqrt(n); 

% call cauchytoeplitz: 
H = hss_cauchytoeplitz(G, F, n, varargin{:}); 
end




