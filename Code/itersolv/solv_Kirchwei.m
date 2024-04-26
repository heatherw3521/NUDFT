function [c w] = solv_Kirchwei(x,N,b,opts)
% SOLV_KIRCHWEI  Kircheis-Potts'23 density weights w/ adjoint soln of t2 NUDFT lin sys
%
% c = solv_kirchwei(x,N,b) approximately solves Ac=b where A is the 1D type 2
%  (forward) NUDFT matrix with nodes x and N modes.
%  Returns c a length-N column vector of complex Fourier coeffs.
%  NU pts x are on 2pi-periodic domain.
%
%  Uses semi-direct method of Kircheis-Potts '23 FAMS paper, Sec 3.1 (Alg 3.6):
%  1) use iterative CG on adj normal eqns for the (2N-1)*M quadrature lin sys for w,
%  which are "density compensation" or quadrature weights for the Euler-F formula.
%  Then 2) c_approx = A^*Wb is directly evaluated (single NUFFT 1d2 call, v fast).
%  This does not exploit the power of the direct aspect of the method (prestoring
%  w would).
%
%  Step 1) solves Bw = e_0,   where e_0 length 2N-1, 1 in the central entry, else 0.
%  and where B_{kj} = e^{ikx_j} is the adjoint (type 1) matrix for twice-bandwidth.
%  This is done by i) solving ANE,   BB^*u = e_0,   via CG with Toeplitz apply, then
%  ii) evaluating w = B^*u.
%
%  Note: M>=2N-1 for now (overdet), so that Toeplitz matrix BB^* is full rank.
%
% [c w] = solv_kirchwei(x,N,b) also outputs the M-component w weight vector.
%
% c = solv_kirchwei(x,N,b,opts) controls various options:
%  opts.tol = NUFFT tolerance for A^* apply
%  opts.cgtol = rel resid for weight CG solve (>opts.tol)
%  opts.wei = force quad wei vec, useful if known eg from [c w] = ... output above.
%           (*** factor 2pi to sort out)
%
% For test: i) run without arguments (only passes nudist=0).
%           ii) test_iNUDFT as tested by a Kirchwei driver.
%
% Reference:
%  Kircheis M and Potts D (2023), "Fast and direct inversion methods for the
%  multivariate nonequispaced fast Fourier transform",
%  Front. Appl. Math. Stat. 9:1155484. doi: 10.3389/fams.2023.1155484
%
% To do: code the CG on (type 1) normal eqns B^*Bw = 1_M,  to handle M<2N-1 too

% Barnett 4/3/24
if nargin==0, test_solv_Kirchwei; return; end
if nargin<4, opts=[]; end
if ~isfield(opts,'tol'), opts.tol = 1e-12; end
if ~isfield(opts,'cgtol'), opts.cgtol = 1e-6; end
if isfield(opts,'wei')            % overrides step 1) below.
  w = (1/(2*pi)) * opts.wei;        % Euler-Fourier formula prefac ??
else
  % 1) build density comp weights via Kircheis-Potts (iterative, may be slow):
  v = finufft1d1(x, ones(size(x)), +1, opts.tol, 4*N-3);  % Toep vec -2N+1:2N-1
  vhat = fft([v;0]);
  e0 = zeros(2*N-1,1); e0(N)=1;    % rhs Kronecker vec for quadr wei sys
  maxit = 1e4;
  % i) iter solve...
  [u,flag,relres,iter,resvec] = pcg(@(u) Toep_apply(u,vhat), e0, opts.cgtol, maxit);
  fprintf('\tget w: CG flag=%d, with %d iters (rel resid nrm %.3g)\n',flag,iter,relres)
  w = finufft1d2(x,-1,opts.tol,u);    % ii) w=B^*u
end
%M = numel(b);
%w = (1/M)*ones(size(b));     % const, exact for x = unif grid over [0,2pi)

% 2) apply weighted adjoint (the direct stage)
c = finufft1d1(x, w.*b, -1, opts.tol, N);  % A^*Wb, note isign=-1 rel to t2


%%%%%
function test_solv_Kirchwei          % simple tester cut down from solv_dense
N = 2e3;             % # modes (unknowns), small problem
rng(0)
c0 = randn(N,1) + 1i*randn(N,1);    % the true coeffs
M = ceil(2.2*N);             % # NU pts (eqns)
for nudist = 0:2  % set up various NU pt distns (easy to hard)............
  fprintf('test_solv_Kirchwei. NU dist=%d:\n',nudist);
  if nudist==0, jitter = 0.49;   % size of rand jitter from grid, wrt grid h
    x = (2*pi/M)*((0:M-1) + jitter*(2*rand(1,M)-1));
  else
    x = 2*pi*rand(1,M);    % iid unif entire domain
  end
  if nudist==2, sbwp = 8.0;   % space-bandwidth prod in half-wavelengths
    x = x*(1-sbwp/N);        % open up a gap w/ no pts in it
  end
  
  tol=1e-14;
  fprintf('fast b and resid test (tol=%.3g, no matrix A)...\n',tol)
  b = finufft1d2(x,+1,tol,c0);
  [c w] = solv_Kirchwei(x,N,b);
  fprintf('\tKirchwei done: sum(w)=%.3g, |w|_1 = %.3g\n',sum(w),sum(abs(w)))
  resid = finufft1d2(x,+1,tol,c) - b;
  fprintf('rel l2 soln err %.3g,   resid rel l2 nrm %3g\n', norm(c-c0)/norm(c0), norm(resid)/norm(b))
end                                         % ..............
% only passes test for nudist=0 (1 and 2 fail to converge solving weights).
