function w = optfrobwei(x,N,opts)
% OPTFROBWEI  optimal (in Frobenius) norm diag gridding "sinc^2" weights
%
% w = optfrobwei(x,N,opts) returns col vector of optimal "gridding" weights
%  to apply before a type 1 NUDFT to make their composition the best Frobenius
%  norm approximation to the inverse of the type 2.
%
%  That is,     w = argmin_{w in R^M} || AA^*diag{w} - I ||_F
%
%  where A is type 2 NUDFT matrix from N modes to the M=numel(x) points x.
%
%  Notes:
%  1) Leslie refers to as sinc^2 weights, but they are only approx sinc^2,
%     rather instead exactly Fejer kernel.
%  2) 1D only for now, N even only. Fixed non-Re issue 11/14/22.
%
% For test see: SOLV_ADJWEI

% Barnett 11/9/22
if nargin<3, opts=[]; end
if ~isfield(opts,'tol'), opts.tol = 1e-12; end

% get s =  "sinc^2 transform" of 1's vector... really a Fejer-kernel xform!
assert(mod(N,2)==0)       % even N for now
f = finufft1d1(x, ones(size(x)), -1, opts.tol, 2*N);   % big t1

% apply counts to turn 1D sum over j=k-k' into double sum
% over k,k' in {-N/2,N/2-1}...  (=Fejer Fourier coeffs)
f = f .* [0:N-1, N:-1:1]';      % symmetric, gives w Re
s = finufft1d2(x, +1, opts.tol, f);    % sums of squares of Fejer kernels
w = N./s;  % Choi-Munson '98 formula for optimal diag wei in Frob nrm
%w = (N-1)./s;   % found for M>N regular x grid, this is exact! why?
w = real(w);                    % since Re
