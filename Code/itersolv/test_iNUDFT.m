function out = test_iNUDFT(signal,N,dataratio,solvlist,nudistlist,opts)
% test_iNUDFT   test inverting 1D type 2 NUDFT, various methods and problems
%
% test_iNUDFT(signal,N,dataratio,solvlist,nudistlist,opts) writes test timings
%  and error outputs to the MATLAB text terminal, for solving the 1D type 2
%  NUDFT MxN linear system on several cases each for several solvers. For large
%  NM the matrix is never formed (dense solver not tried) and the residual is
%  also found fast using a NUFFT.
%
% All arguments are optional and otherwise take default values:
%    signal - what 1D func to recon: 'randn' = randn F series coeffs.
%                                    'decay' = as above but decay to exp(-5)
%                                    'image' = test 1d image w/ features
%    N - how many unknowns
%    dataratio - ratio of data points (equation) to unknowns, ie, controls M/N
%    solvlist - a subset of the natural numbers listing which solvers to run
%               (see code)
%    nudistlist - a subset of {1,2,3,4} listing which NU point distributions
%               (see code)
%    opts - a struct containing maybe:
%           opts.tol (for various solvers)
%           opts.cgtol (for iterative solvers)
%           opts.verb = 0,1,.. verbosity (0=usual text, 1=also figs)
%           opts.addnoise set iid noise level added to data. (0 default)
%           opts.outputstruct - if nonzero returns 2d cell array of results
%
% A good usage it to capture this to a diary file:
%  diary('mytest')
%  test_iNUDFT(signal,N,dataratio,solvlist,nudistlist,opts)
%  diary off
%
% todo:
%  * study cons vs incons b; noiseless/noisy.
%  * empty args make take default
%  * struct of results, including residueal vecs, timings, for plotting
%    also soln vec nrm
%  * find better Toep preconds than Strang - check lit
  
% Barnett 11/7/22, onwards. Repackaged as func so can be looped, 2/2/23.

% matching list of solvers, solver-specific opts, names, all same interface...
solvers =  {@solv_dense,    @solv_CGN,        @solv_CGN,...
            @solv_adjwei,         @solv_FPadjwei, @solv_CGAN,...
            @solv_CGAN,                @solv_CGAN,           @wrapper_INUDFT,...
            @solv_FDToepN,    @solv_Kirchwei};
solvopts = {struct(),       struct(),         struct('precond','strang'),...
            struct(),             struct(),       struct(),...
            struct('precond','sinc2'), struct('gmres',1),    struct(),...
            struct(),             struct()};
solvnams = {"dense direct", "CG normal eqns", "Strang PCG nor eqns",...
            "adj matvec sinc2/quadr wei","FP adj wei",  "CG adj nor eqns",...
            "sinc2 PCG adj nor",      "GMRES adj nor eqns", "INUDFT, wrapped"...
            "FDToep nor eqns",  "Kirch adj wei (by CGAN)"};
if nargin<4, solvlist=[1:7,9:11]; end     % med size, exclude 8=GMRES (slow)

% overall problem size and target tolerance...
if nargin<2, N = 2^14; end   % # modes (unknowns) ...FDS/HSS seem to require 2^k
rng(0)
if nargin<1, signal='randn'; end      % signal type
freqinds = (-N/2:N/2-1)';           % col vec
if strcmp(signal, 'randn')          % choose signal: c0 = true Fourier coeffs
  c0 = randn(N,1) + 1i*randn(N,1);     % flat spectrum
elseif strcmp(signal, 'decay')
  c0 = exp(-10.0*abs(freqinds/N)) .* (randn(N,1) + 1i*randn(N,1));  % exp decay
elseif strcmp(signal, 'image')     % the F coeffs define this func, not fun(x)
  fun = @(x) (x<3).*sin(8*x.^2)  + (abs(x-4)<.5) + max(0,1-2*abs(x-5.5));
  nquad = 10*N;            % to be somewhat accurate for Euler-F quadrature
  xq = (2*pi/nquad)*(0:nquad-1);
  fq = fun(xq);
  tol = 1e-14; c0 = finufft1d1(xq,fq*(1/nquad),-1,tol,N);  % do the quad
  %c0 = c0 .* (1 - (freqinds/(N/2)).^2);    % Welch filter
  e=3.5; c0 = c0 .* exp(-0.5*(freqinds*(e/(N/2))).^2); % trunc Gauss filt 1e-3
end
if ~isfield(opts,'outputstruct'), opts.outputstruct=0; end
if opts.outputstruct %we want to build a struct to store results
    out = {}; 
else
    out = []; 
end
if nargin<6, opts = []; end
if ~isfield(opts,'tol'), opts.tol = 1e-12; end   % tol for FD meths
if ~isfield(opts,'cgtol'), opts.cgtol = 1e-6; end    % be nice to iter meths :)
if ~isfield(opts,'verb'), opts.verb = 0; end
if ~isfield(opts,'addnoise'), opts.addnoise = 0.0; end
if nargin<3, dataratio = 2.0; end       % how tall lin sys is (=M/N)
M = ceil(N*dataratio);     % # NU pts (eqns), very sensitive to its ratio to N
fasttest = (M*N>=1e7);        % whether run fast or dense checks
if nargin<5 | isempty(nudistlist), nudistlist = 1:4; end       % default
jj = 0; 
for nudist = nudistlist % loop over NU pt distns for x (easy to hard)..........
    jj = jj+1;
    if nudist==1, nunam='jittered grid';
    jitter = 0.49;   % size (h units) of jitter from grid; <=0.5 well-cond
    x = (2*pi/M)*((0:M-1) + jitter*(2*rand(1,M)-1));
  elseif nudist==2, nunam='quadrature pts';     % Cheby nodes on [0,2pi]
    [x w] = clencurt(M-1);
    x = pi*(1+x'); opts.w = pi*w';  % scale to domain [0,2pi]  (w=quadwei)
  elseif nudist==3, nunam='rand unif iid';
    x = 2*pi*rand(1,M);           % iid unif entire domain
  elseif nudist==4, nunam='rand iid w/ gap';
    sbwp = 8.0;   % gap space-bandwidth prod in half-wavelengths
    x = 2*pi*rand(1,M)*(1-sbwp/N);    % gap w/ no pts in it (exp bad wrt sbwp)
  end
  fprintf('\n------------- NU dist=%d (%s): -----------------\n',nudist,nunam);
  np = 6*N; xp = 2*pi*(1:np)/np;   % plot grid
  y0p = finufft1d2(xp,+1,opts.tol,c0);  % signal on plot grid
  if opts.verb, figure(nudist); clf; plot(xp,y0p,'-'); hold on; end
  
  if fasttest                % set up and test in fast way, to tol opts.tol
    b = finufft1d2(x,+1,opts.tol,c0);
    b = b + opts.addnoise*randn(M,1);
    fprintf('\nfast b and resid meas (tol=%.3g, no matrix A)...\n',opts.tol)
    if opts.verb, plot(x,real(b),'k.'); end
    for s=solvlist
      fprintf("solver %d (%s)...\n",s,solvnams{s})
      o = mergestructs(opts,solvopts{s});
      t = tic;
      c = solvers{s}(x,N,b,o);
      tim = toc(t);
      Ac = finufft1d2(x,+1,opts.tol,c);
      resid = Ac - b;
      fprintf('\t%.3g s   \trel l2 err %.3g    \tresid rel l2 nrm %.3g\n', tim, norm(c-c0)/norm(c0), norm(resid)/norm(b))
      if opts.verb, yp = finufft1d2(xp,+1,opts.tol,c);  % recon on plot grid
        plot(xp,real(yp),'-'); end
    end
  else                       % dense O(NM) set up and test, to eps_mach
    fprintf('\ndense (slow exact) b and resid meas... ')
    A = densemat_nudft(x,N);
    if M*N<1e7, fprintf('kappa(A) = %.3g\n',cond(A)); end   % a luxury
    %if ~isfield(opts,'w') && M*N<1e7, opts.w = A\eye(N)???; end   % compute weights
    b = A*c0;
    b = b + opts.addnoise*randn(M,1);
    if opts.verb, plot(x,b,'k.'); end
    kk = 0; 
    for s=solvlist
      kk = kk+1; 
      fprintf("solver %d (%s)...\n",s,solvnams{s})
      o = mergestructs(opts,solvopts{s});
      t = tic;
      c = solvers{s}(x,N,b,o);
      Ac = A*c;
      fprintf('\t%.3g s   \trel l2 err %.3g    \tresid rel l2 nrm %.3g\n', toc(t), norm(c-c0)/norm(c0), norm(Ac-b)/norm(b))
      if opts.verb, yp = finufft1d2(xp,+1,opts.tol,c);  % recon on plot grid
          plot(xp,real(yp),'-'); end
      if opts.outputstruct
        solstruc = []; 
        if opts.verb, solstruc.yp = yp; end % equi point sig vals
        solstruc.xp = xp; % equi points
        solstruc.y = Ac; % sample point sig vals
        solstruc.x = x; % sample points
        solstruc.err = norm(c-c0)/norm(c0); 
        solstruc.resid = norm(Ac-b)/norm(b);
        solstruc.rhs = b; 
        solstruc.sig = y0p; %true sig on grid
        solstruc.cfs = c0; % true coeffs
        solstruc.pts = nunam;
        solstruc.solver = solvnams{s};
        out{jj,kk} = solstruc; 
      end
    end
  end
  if opts.verb, legend(['true', 'data', solvnams{solvlist}]); end
end                                                        % ..............

%out = [];   % ***
