function [c,data] = solv_CGN(x,N,b,opts)
% SOLV_CGN  conjugate gradients on normal eqns to solve 1D NUDFT type 2 lin sys
%
% c = solv_CGN(x,N,b) solves Ac=b where A is the 1D type 2 (forward) NUDFT
%  matrix with nodes x and N modes. Returns c a length-N column vector of
%  complex Fourier coeffs. NU pts x are on 2pi-periodic domain.
%
%  (P)CG on normal equations with FFT-apply of Toeplitz A^*A is used.
%
% c = solv_CGN(x,N,b,opts) controls various options:
%  opts.tol : overall tolerance (in sense of A and b in lin sys, ie, residual)
%  opts.precond : 'strang' for Strang's circulant preconditioner in PCG
%
% For test: see test_iNUDFT

% Barnett 11/7/22. Fixed for faster Toep_apply, 1/31/24
if nargin<4, opts=[]; end
if ~isfield(opts,'tol'), opts.tol = 1e-12; end
if ~isfield(opts,'cgtol'), opts.cgtol = opts.tol; end

rhs = finufft1d1(x, b, -1, opts.tol, N);  % A^* b,  note isign flipped rel to t2

% vector v defining Toeplitz matrix A^*A...
v = finufft1d1(x, ones(size(x)), -1, opts.tol, 2*N-1);  % indices -(N-1):(N-1)
%norm(v-conj(v(end:-1:1)),inf)      % check Toep is Hermitian
% Note padding here for new (non-bad) version of Toep_apply (swapped args too!):
vhat = fft([v;0]);   % only needs doing once (would be Re if v were circshifted)
data.vhat = vhat;

%T = toeplitz(v(N:end),v(N:-1:1)); lam=eig(T); fprintf('\t\trange(eig(T)) = [%.3g,%.3g], kappa(T)=%.3g\n',min(lam),max(lam),max(lam)/min(lam))   % debug, since SPD

%maxit = min(N,2e3);       % please edit for tests or make an opts.maxit
if ~isfield(opts, 'maxit')
    maxit = min(N,2e3);  %M, the max
else
    maxit = opts.maxit; 
end

if ~isfield(opts,'precond')        % no precond, plain CG
  c = zeros(N,size(rhs,2));
  for i = 1:size(rhs,2)
    [c(:,i),flag,relres,iter,resvec] = pcg(@(a) Toep_apply(a,vhat), rhs(:,i), opts.cgtol, maxit);
  end
  
elseif strcmp(opts.precond,'strang')
  % we don't expect Strang to help because Toep A^*A no off-diag decay...
  [Cinv,Cinvhat] = strangprecond(v);     % see tester for this func for examples
  data.Cinvhat = Cinvhat;
  [c,flag,relres,iter,resvec] = pcg(@(a) Toep_apply(a,vhat), rhs, opts.cgtol, maxit, @(a) circ_apply(Cinvhat,a));   % last arg is matvec func for P^{-1}

  if 0, Ci = circulant(Cinv);      % some dense debugging of keeping Hermitian
    %fprintf('Cinv SPD? %.3g %.3g\n',norm(Ci-Ci','fro'), min(eig(Ci)))
    T = toeplitz(v(N:end),v(N:-1:1)); 
    CinvT = Ci*T;      % debug dense
    %fprintf('\t\tkappa(C^{-1}T)=%.3g\n',cond(CinvT))
    Cis = circulant(ifft(sqrt(Cinvhat)));  % C^{-1/2}
    Tp = Cis*T*Cis';
    fprintf('\t\tkappa(C^{-1/2}TC^{*,-1/2})=%.3g\n',cond(Tp))
    figure; plot(eig(T),'.'); axis equal; hold on; plot(eig(Tp),'+');
    norm(Tp-Tp',inf)
  end
  
else
  error('no such preconditioner')
end
fprintf('\tflag=%d, with %d iters (rel resid nrm %.3g)\n',flag,iter,relres)
