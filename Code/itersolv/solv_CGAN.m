function [c,data] = solv_CGAN(x,N,b,opts)
% SOLV_CGAN  conjugate gradients on adjoint normal eqns, 1D NUDFT type 2 lin sys
%
% c = solv_CGAN(x,N,b) solves Ac=b where A is the 1D type 2 (forward) NUDFT
%  matrix with nodes x and N modes. Returns c a length-N column vector of
%  complex Fourier coeffs. NU pts x are on 2pi-periodic domain.
%
%  Uses CG (or GMRES) on adjoint normal ("paranormal") equations (AA^*)y = b,
%   with NUFFT pair applying AA^*, finally returning c = A^*y.
%
% c = solv_CGAN(x,N,b,opts) controls various options:
%  opts.tol : overall tolerance (in sense of A and b in lin sys, ie, residual)
%  opts.cgtol : CG stopping criterion (relative residual)
%  opts.precond : 'sinc2' - use Leslie's optimal (in Frob nrm) diag precond
%  opts.gmres : if true, replace CG with GMRES (without precond).
%
% For test: see test_iNUDFT

% Barnett 11/12/22
if nargin<4, opts=[]; end
if ~isfield(opts,'tol'), opts.tol = 1e-12; end
if ~isfield(opts,'cgtol'), opts.cgtol = opts.tol; end

% vector v defining Toeplitz matrix A^*A...
v = finufft1d1(x, ones(size(x)), -1, opts.tol, 2*N-1);  % indices -(N-1):(N+1)
%norm(v-conj(v(end:-1:1)),inf)      % check Toep is Hermitian
vhat = fft(v);   % only needs doing once (would be Re if v were circshifted)
data = struct();

M = numel(x);
if ~isfield(opts, 'maxit')
    maxit = min(M,2e3);  %M, the max
else
    maxit = opts.maxit; 
end

if isfield(opts,'gmres') && opts.gmres  % override by GMRES
  [y,flag,relres,iter,resvec] = gmres(@(a) AAH_apply(x,a,N,opts.tol), b, [], opts.cgtol, M);
  iter = iter(end);
  
elseif ~isfield(opts,'precond')        % no precond, plain CG
  [y,flag,relres,iter,resvec] = pcg(@(a) AAH_apply(x,a,N,opts.tol), b, opts.cgtol, maxit);
  
elseif strcmp(opts.precond,'sinc2')        % Leslie optimal diag precond
  w = optfrobwei(x,N,opts);    % get 1/(sum of sinc^2) wei col vec, is inv(P)
  data.w = w;
  [y,flag,relres,iter,resvec] = pcg(@(a) AAH_apply(x,a,N,opts.tol), b, opts.cgtol, maxit, @(a) w.*a);
  
else
  error('no such preconditioner')
end
fprintf('\tflag=%d, with %d iters (rel resid nrm %.3g)\n',flag,iter,relres)

% postproc to the solution...
c = finufft1d1(x, y, -1, opts.tol, N);  % A^*y,  note isign flipped rel to t2

%%%%%%
function AAHa = AAH_apply(x,a,N,tol)
% multiplies vector a by A^* then by A, via NUFFTs. Matvec for adj normal eqns
AHa = finufft1d1(x, a, -1, tol, N);
AAHa = finufft1d2(x, +1, tol, AHa);
