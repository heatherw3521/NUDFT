function [x,w] = clencurt(N)
% CLENCURT  nodes (Chebyshev) and weights, Clenshaw-Curtis quadrature, via FFT
%
% [x,w] = clencurt(N) returns N+1 nodes x (column vector), and N+1 weights
%   (row vector), such that for a smooth function on [-1,1] with vectorized
%   handle f, the integral of f on [-1,1] is well approximated by w*sum(f(x)).
%
%  Notes: Trefethen book, modified by Barnett to return x in increasing order.
%  FFT version using Fourier series for |sin(theta)|, by Barnett

% Copyright (C) 2008--2010, Alex Barnett and Timo Betcke

theta = pi*(N:-1:0)'/N; x = cos(theta); % note order opposite to Trefethen
W = kron(-1./((1:floor(N/2)).^2-1/4), [0 1]);   % works for even or odd
if mod(N,2)==1, W = [W 0]; end  % include extra pi-freq term if odd
w = ifft([4 W W(end-1:-1:1)]);  % 4 is the zero-freq term
w = [w(1)/2 w(2:N) w(1)/2];     % endpoints get 1/2 weight since want 1/2 circle
