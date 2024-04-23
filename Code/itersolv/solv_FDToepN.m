function [c,data] = solv_FDToepN(x,N,b,opts)
% SOLV_FDTOEPN  fast direct Toeplitz on normal eqns, 1D NUDFT type 2 lin sys
%
% c = solv_FDToepN(x,N,b) solves Ac=b where A is the 1D type 2 (forward) NUDFT
%  matrix with nodes x and N modes. Returns c a length-N column vector of
%  complex Fourier coeffs. NU pts x are on 2pi-periodic domain.
%
%  Wilber-Epperly-modified HM-Toolbox fast direct Toeplitz solver on normal
%  equations is used, with type-1 NUFFT to set up the Toeplitz vector and RHS.
%
% c = solv_FDToepN(x,N,b,opts) controls various options:
%  opts.tol : overall tolerance (in sense of A and b in lin sys, ie, residual)
%  opts.hsstol : HSS tolerance
%
% For test: see test_iNUDFT

% Barnett 11/18/22
if nargin<4, opts=[]; end
if ~isfield(opts,'tol'), opts.tol = 1e-12; end
if ~isfield(opts,'hsstol'), opts.hsstol = 10*opts.tol; end   % seems wise

rhs = finufft1d1(x, b, -1, opts.tol, N);  % A^* b,  note isign flipped rel to t2

% vector v defining Toeplitz matrix A^*A...
v = finufft1d1(x, ones(size(x)), -1, opts.tol, 2*N-1);  % indices -(N-1):(N+1)
tc = v(N:end); tr = v(N:-1:1);        % extract Matlab-form Toep col and row

try
  [c,data.H] = Toeplitz_solve_ns(tc,tr,rhs, 'tol', opts.hsstol);  % wrap mod HM-Toolbox
catch ME
  fprintf('\tToeplitz_solve failed!\n\t%s\n',ME.message);
  c = nan(N,1);
end

%fprintf('\tcheck nor eqns rel resid %.3g\n',norm(Toep_apply(c,fft([v;0]))-rhs)/norm(rhs))     % fun but optional
