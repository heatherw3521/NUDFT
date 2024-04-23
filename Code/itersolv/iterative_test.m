function [relresid, relerr, times, numits, condA] = iterative_test(signal,N,dataratio,x, wts, condflag)
%  MODIFICATION OF AB CODE test_iNUDFT. This forms a figure showing the 
% performance of the iterative methods as the conditioning of the matrix
% gets increasingly worse.
%
% output = relative residual for all solvers, relativer error, and
% condition number of Vandermonde for this choice of nodes. 
%
% All arguments are optional and otherwise take default values:
%    signal - what 1D func to recon: 'randn' = randn F series coeffs.
%                                    'decay' = as above but decay to exp(-5)
%                                    'image' = test 1d image w/ features
%    N - how many unknowns
%    dataratio - ratio of data points (equation) to unknowns, ie, controls M/N
%    solvlist - a subset of the natural numbers listing which solvers to run
%               (see code)
%    
%    opts - a struct containing maybe:
%           opts.tol (for various solvers)
%           opts.cgtol (for iterative solvers)
%           opts.verb = 0,1,.. verbosity (0=usual text, 1=also figs)
%           opts.addnoise set iid noise level added to data. (0 default)
%           opts.jitter = amount of jitter

% Output is some kind of cell array of structs, TBD ***.
%
% A good usage it to capture this to a diary file:
%  diary('mytest')
%  test_iNUDFT(signal,N,dataratio,solvlist,nudistlist,opts)
%  diary off
%
% matching list of solvers, solver-specific opts, names, all same interface...
solvers =  {@solv_dense,    @solv_CGN,        @solv_CGN,...
            @solv_adjwei,         @solv_FPadjwei, @solv_CGAN,...
            @solv_CGAN,                @solv_CGAN,           @wrapper_INUDFT,...
            @wrapper_INUDFT, @solv_FDToepN};
solvopts = {struct(),       struct(),         struct('precond','strang'),...
            struct(),             struct(),       struct(),...
            struct('precond','sinc2'), struct('gmres',1),    struct(), struct(),...
            struct()};
solvnams = {"dense direct", "CG normal eqns", "Strang PCG nor eqns",...
            "adj matvec sinc2/quadr wei","FP adj wei",  "CG adj nor eqns",...
            "sinc2 PCG adj nor",      "GMRES adj nor eqns", "INUDFT, wrapped", "INUDFT, smalltol"...
            "FDToep nor eqns"};
%if nargin<4, solvlist=[1:7,9:10]; end     % med size, exclude 8=GMRES (slow)
%solvlist = [2:3,5:7,9];
solvlist = 2;
%solvlist = [9,10]; 
% overall problem size and target tolerance...
%if nargin<2, N = 2^13; end   % # modes (unknowns) ...FDS/HSS seem to require 2^k
rng(0)
%if nargin<1, signal='randn'; end      % signal type
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
opts = []; 
%if nargin<6, opts = []; end
if ~isfield(opts,'tol'), opts.tol = 1e-12; end   % tol for FD meths
if ~isfield(opts,'cgtol'), opts.cgtol = 1e-7; end    % be nice to iter meths :)
if ~isfield(opts,'verb'), opts.verb = 0; end
if ~isfield(opts,'addnoise'), opts.addnoise = 0.0; end
if nargin<3, dataratio = 2.0; end       % how tall lin sys is (=M/N)
M = ceil(N*dataratio);     % # NU pts (eqns), very sensitive to its ratio to N
if ~isempty(wts)
    opts.w = pi*wts';
end

%if nargin<5 | isempty(nudistlist), nudistlist = 1:4; end       % default 
%for nudist = nudistlist % loop over NU pt distns for x (easy to hard)..........
  
  %jitter = (nudist/4);   % size (h units) of jitter from grid; <=0.5
  %well-cond
 % x = (2*pi/M)*((0:M-1) + jitter*(2*rand(1,M)-1));
  
  %fprintf('\n------------- NU dist=%d (%s): -----------------\n',1/nudist,nunam);
  %np = 6*N; xp = 2*pi*(1:np)/np;   % plot grid
  %y0p = finufft1d2(xp,+1,opts.tol,c0);  % signal on plot grid
  %if opts.verb, figure(nudist); clf; plot(xp,y0p,'-'); hold on; end
    if ~condflag 
        b = finufft1d2(x,+1,opts.tol,c0);
        condA = 0; 
    else
        A = densemat_nudft(x,N);
        condA = cond(A);
        b = A*c0;
    end
    %if M*N<1e7, fprintf('kappa(A) = %.3g\n',cond(A)); end   % a luxury
    %if ~isfield(opts,'w') && M*N<1e7, opts.w = A\eye(N)???; end   % compute weights
    
    %b = b + opts.addnoise*randn(M,1);
    %if opts.verb, plot(x,b,'k.'); end
    kk = 0; 
    for s=solvlist
      kk = kk+1; 
      if s==9 %INUDFT with a bit more than cg tol (to get global approx cgtol)
        opts.hsstol = (1e-1)*(opts.cgtol); 
      elseif s==10 %INDUFT with higher tol 
        opts.hsstol = opts.tol; 
      end
      fprintf("solver %d (%s)...\n",s,solvnams{s})
      o = mergestructs(opts,solvopts{s});
      t = tic;
      [iter, c] = solvers{s}(x,N,b,o);
      tim = toc(t);
      Ac = finufft1d2(x,+1,opts.tol,c);
      relresid(kk) = norm(Ac-b)/norm(b); 
      relerr(kk) = norm(c-c0)/norm(c0);
      times(kk) = tim; 
      numits(kk) = iter; 
     end
  %if opts.verb, legend(['true', 'data', solvnams{solvlist}]); end
%end                                                        % ..............

