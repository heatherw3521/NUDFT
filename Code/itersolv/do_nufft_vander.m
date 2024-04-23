function c = do_nufft_vander(x,f,tol)
% DO_NUFFT_VANDER  Wrapper to FINUFFT 1d type-2 from Vandermonde freq convention
%
% c = do_nufft_vander(x,f) uses x as real NU pts in [0,2pi), f as complex coeffs
%  for the freqs 0,1,...,n-1, and performs a 1D type-2 fast transform.
%
% For test see: DENSEMAT_NUDFT_VANDER
if nargin<3, tol=1e-12; end
N = numel(f);
assert(mod(N,2)==0)             % same as FINUFFT "N"
c = finufft1d2(x,+1,tol,f);
c = c .* exp(1i*(N/2)*x);       % correct for the one-sided freq offset
