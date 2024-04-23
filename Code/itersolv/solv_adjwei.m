function c = solv_adjwei(x,N,b,opts)
% SOLV_ADJWEI  diag-weighted adjoint "gridding" approx soln of t2 NUDFT lin sys
%
% c = solv_adjwei(x,N,b) approximately solves Ac=b where A is the 1D type 2
%  (forward) NUDFT matrix with nodes x and N modes.
%  Returns c a length-N column vector of complex Fourier coeffs.
%  NU pts x are on 2pi-periodic domain.
%
%  Uses the "gridding" method of diagonal scaling (by diag matrix W) then
%  mat-vec with adjoint (A^*, fast via t1 NUFFT). Ie, capprox = A^*Wb.
%  Computes weights by fast method, but on the fly (better to prestore?)
%
% c = solv_adjwei(x,N,b,opts) controls various options:
%  opts.tol = NUFFT tolerance for A^* apply and weight computation
%  opts.w = force quad wei, useful if known
%
% For test: see test_iNUDFT

% Barnett 11/8/22
if nargin<4, opts=[]; end
if ~isfield(opts,'tol'), opts.tol = 1e-12; end
if isfield(opts,'w')            % this opt overrides optfrobwei
  w = (1/(2*pi)) * opts.w;        % Euler-Fourier formula prefac
else
  % build the vector of weights via discrete version of Leslie's sinc^2 xform:
  w = optfrobwei(x,N,opts);
end
%M = numel(b);
%w = (1/M)*ones(size(b));     % const, exact for x = unif grid over [0,2pi)

% apply weighted adjoint
c = finufft1d1(x, w.*b, -1, opts.tol, N);  % A^*Wb, note isign=-1 rel to t2
