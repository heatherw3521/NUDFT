% Test script for inverting 1D type 2 NUDFT: various methods on various problems
% Barnett 11/7/22, onwards.
clear

% matching list of solvers, solver-specific opts, names, all same interface...
%solvers =  {@solv_dense};
solvers = {@wrapper_INUDFT};
solvopts = { struct() };
solvnams = {"INUDFT"};
%solvlist=[1:7,9:10];      % med size, exclude 8=GMRES, slow
%solvlist=1:10;         % all (for 2^10 small probs)
%solvlist=[2:7,10];   % list for 2^18 big probs
solvlist = [1];
% overall problem size and target tolerance...
N = 1800;             % # modes (unknowns) ...FDS/HSS seem to require 2^k
rng(0)
c0 = randn(N,1) + 1i*randn(N,1);    % true Fourier coeffs (no Fourier decay yet)
M = ceil(N*2.0);     % # NU pts (eqns), very sensitive to its ratio to N
%M = N; 
opts.tol = 1e-12;
opts.cgtol = 1e-6;    % be nice to iter for now
%fasttest = (M*N>=1e7);        % whether run fast or dense checks
fasttest = false;

for nudist = 1:4   % loop over NU pt distns for x (easy to hard)............
  if nudist==1, nunam='jittered grid';
    jitter = 0.5;   % size (h units) of jitter from grid; <=0.5 well-cond
    x = (2*pi/M)*((0:M-1) + jitter*(2*rand(1,M)-1));
  elseif nudist==2, nunam='quadrature pts';
    [x w] = clencurt(M-1);
    x = pi*(1+x'); opts.w = pi*w';    % covers entire domain [0,2pi)
  elseif nudist==3, nunam='rand unif iid';
    x = 2*pi*rand(1,M);           % iid unif entire domain
  elseif nudist==4, nunam='rand iid w/ gap';
    sbwp = 8.0;   % gap space-bandwidth prod in half-wavelengths
    x = 2*pi*rand(1,M)*(1-sbwp/N);    % gap w/ no pts in it (exp bad wrt sbwp)
  end
  fprintf('\n------------- NU dist=%d (%s): -----------------\n',nudist,nunam);
  
  if fasttest                % set up and test in fast way, to tol
    tol=1e-15;
    b = finufft1d2(x,+1,opts.tol,c0);
    fprintf('\nfast b and resid meas (tol=%.3g, no matrix A)...\n',tol)
    for s=solvlist
      fprintf("solver %d (%s)...\n",s,solvnams{s})
      o = mergestructs(opts,solvopts{s});
      t = tic;
      c = solvers{s}(x,N,b,o);
      tim = toc(t);
      resid = finufft1d2(x,+1,opts.tol,c) - b;
      fprintf('\t%.3g s   \trel l2 err %.3g    \tresid rel l2 nrm %.3g\n', tim, norm(c-c0)/norm(c0), norm(resid)/norm(b))
    end
  else                       % dense O(NM) set up and test, to eps_mach
    fprintf('\nINDUFT  b and resid meas... ')
    A = densemat_nudft(x,N);
    if M*N<1e7, fprintf('kappa(A) = %.3g\n',cond(A)); end   % a luxury
    b = A*c0;
    for s=solvlist
%      fprintf("solver %d (%s)...\n",s,solvnams{s})
      o = mergestructs(opts,solvopts{s});
      t = tic;
      c = solvers{s}(x,N,b,o);
      fprintf('\t%.3g s   \trel l2 err %.3g    \tresid rel l2 nrm %.3g\n', toc(t), norm(c-c0)/norm(c0), norm(A*c-b)/norm(b))
    end
  end
end                                                        % ..............

