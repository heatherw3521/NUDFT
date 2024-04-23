function V = densemat_nudft(x,N)
% DENSEMAT_NUDFT   fill dense 1D type-2 NUDFT matrix (FINUFFT convention)
%
% V = densemat_nudft(x,N) returns M*N NUDFT matrix where
%  M = numel(x) and x is a list (row or col vec) of nonuniform points x_j,
%  j=1,..,M, on real axis. FINUFFT convention (isign=+1) is used, ie
%  the x domain is 2pi-periodic, and
%  V_{jk} = exp(i*(-N/2-1+k)*x_j)   for j=1...M and k=1...N (1-indexed matrix)
%
%  Without arguments does self test

% Barnett 11/7/22
if nargin==0, test_densemat_nudft; return, end
freqs = -floor(N/2) + (0:N-1);      % row vec
V = exp(1i*x(:)*freqs);             % outer prod  (cos + i.sin would be faster)

%%%%%%%
function test_densemat_nudft    % make sure conventions translate
M = 1000;             % # NU pts
for N = [500 501]     % test even and odd cases
  x = 2*pi*rand(M,1);   % NU pts
  V = densemat_nudft(x,N);
  f = randn(N,1) + 1i*randn(N,1);    % some coeffs
  cdense = V*f;
  tol=1e-12;
  c = finufft1d2(x,+1,tol,f);
  fprintf('rel l2 err of Vander dense vs FINUFFT: %.3g\n', norm(c-cdense)/norm(c))
end
