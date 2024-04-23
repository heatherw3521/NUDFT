function V = densemat_nudft_vander(ucnodes,n)
% DENSEMAT_NUDFT_VANDER   dense type-2 NUDFT matrix in Vandermonde convention
%
% V = densemat_nudft_vander(ucnodes,n) returns M*n NUDFT matrix where
%  M = numel(ucnodes) and ucnodes is a list (row or col vec) of nodes gamma_j,
%  j=1,..,M, on unit circle.
%
%  The convention is Vandermonde:
%  V_{jk} = (gamma_j)^{k-1},     k=1,..,n,  j=1,..,M.
%  Thus there are "one-sided" Fourier freqs 0,1,..,n-1.
%  The nonuniform points have been mapped to unit circle (gamma_j).
%
%  Without arguments does self test

% Barnett 11/4/22
if nargin==0, test_densemat_nudft_vander; return, end
V = (ucnodes(:)).^(0:n-1);     % according to Wilber notes defn

%%%%%%%
function test_densemat_nudft_vander    % make sure conventions translate
M = 1000;             % # NU pts
N = 500;              % # freqs (even)
x = 2*pi*rand(M,1);    % NU pts: note x neq lambda's
gamma = exp(1i*x);
V = densemat_nudft_vander(gamma,N);
f = randn(N,1) + 1i*randn(N,1);    % some coeffs
cdense = V*f;
tol=1e-12;
c = do_nufft_vander(x,f,tol);
fprintf('rel l2 err of Vander dense vs FINUFFT: %.3g\n', norm(c-cdense)/norm(c))
