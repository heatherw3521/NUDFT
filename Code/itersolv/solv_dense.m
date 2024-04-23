function c = solv_dense(x,N,b,opts)
% SOLV_DENSE  plain direct solution of 1D NUDFT (type 2) linear system
%
% c = solv_dense(x,N,b) solves Ac=b where A is the 1D type 2 (forward) NUDFT
%  matrix with nodes x and N modes. Returns c a length-N column vector of
%  complex Fourier coeffs. NU pts x are on 2pi-periodic domain.
%
% Without arguments, does self-tests for rect case.

% Barnett 11/7/22
if nargin==0, test_solv_dense; return; end
% opts unused
if numel(x)*N >= 1e7, c = nan(N,1);
  fprintf('\tdecided too big, solv_dense not attempted!\n');
  return; end

A = densemat_nudft(x,N);
c = A\b;


%%%%%
function test_solv_dense
N = 500;             % # modes (unknowns)
rng(0)
c0 = randn(N,1) + 1i*randn(N,1);    % the true coeffs
M = 800;             % # NU pts (eqns)
for nudist = 0:2  % set up various NU pt distns (easy to hard)............
  fprintf('NU dist=%d:\n',nudist);
  if nudist==0, jitter = 0.5;   % size of rand jitter from grid, wrt grid h
    x = (2*pi/M)*((0:M-1) + jitter*(2*rand(1,M)-1));
  else
    x = 2*pi*rand(1,M);    % iid unif entire domain
  end
  if nudist==2, sbwp = 8.0;   % space-bandwidth prod in half-wavelengths
    x = x*(1-sbwp/N);        % open up a gap w/ no pts in it
  end
  
  tol=1e-15;
  fprintf('fast b and resid test (tol=%.3g, no matrix A)...\n',tol)
  b = finufft1d2(x,+1,tol,c0);
  c = solv_dense(x,N,b);
  resid = finufft1d2(x,+1,tol,c) - b;
  fprintf('rel l2 soln err %.3g,   resid rel l2 nrm %3g\n', norm(c-c0)/norm(c0), norm(resid)/norm(b))
  
  fprintf('dense (exact) b and resid test...\n')
  A = densemat_nudft(x,N);
  b = A*c0;
  c = solv_dense(x,N,b);
  fprintf('rel l2 soln err %.3g,   resid rel l2 nrm %3g,   kappa(A)=%.3g\n\n', norm(c-c0)/norm(c0), norm(A*c-b)/norm(b), cond(A))
end                                         % ..............
