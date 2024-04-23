function c = solv_FPadjwei(x,N,b,opts)
% SOLV_FPADJWEI  fixed point iter on diag-weighted adj, solves t2 NUDFT lin sys
%
% c = solv_FPadjwei(x,N,b) solves Ac=b where A is the 1D type 2
%  (forward) NUDFT matrix with nodes x and N modes.
%  Returns c a length-N column vector of complex Fourier coeffs.
%  NU pts x are on 2pi-periodic domain.
%
%  Uses the diagonal weights (diag matrix W) plus
%  mat-vec with adjoint (A^*, fast via t1 NUFFT), in a fixed-point iteration
%  described in Inati,...,Greengard IEEE TMI draft 2006.
%  Computes weights by fast method, once.
%
% c = solv_PFadjwei(x,N,b,opts) controls various options:
%  opts.tol = NUFFT and weight computation tolerance
%  opts.cgtol = stopping relative residual criterion
%
% For test: see test_iNUDFT

% Barnett 11/14/22
if nargin<4, opts=[]; end
if ~isfield(opts,'tol'), opts.tol = 1e-12; end
if ~isfield(opts,'cgtol'), opts.cgtol = opts.tol; end
%maxit = 2e3;   % can exceed M, can be anything...
if ~isfield(opts, 'maxit')
    maxit = 2e3;  %M, the max
else
    maxit = opts.maxit; 
end

% build the vector of weights via discrete version of Leslie's sinc^2 xform:
w = optfrobwei(x,N,opts);

nb = norm(b);
relres = 1;
r = b;
y = 0*b;
resvec = nan*b;
k=1;
while k<maxit && relres > opts.cgtol
  y = y + w.*r;
  AHy = finufft1d1(x, y, -1, opts.tol, N);  % A^* y, note isign=-1 rel to t2
  AAHy = finufft1d2(x, +1, opts.tol, AHy);  % AA^* y
  r = b - AAHy;
  relres = norm(r)/nb;
  resvec(k) = relres;
  %fprintf('k=%d:\trelres=%.3g   \t||y||_2=%.3g\n',k,relres,norm(y)) % debug
  k = k+1;
end
resvec = resvec(~isnan(resvec));            % clip

if k==maxit, fprintf('\tstopped at maxit=%d (rel resid nrm %.3g)\n',k,relres)
else, fprintf('\tdone with %d iters (rel resid nrm %.3g)\n',k,relres), end

c = finufft1d1(x, y, -1, opts.tol, N);      % A^* y, note isign=-1 rel to t2
