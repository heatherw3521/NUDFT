%test iterative methods: 
% use AB's test suite to  compare iterative methods
% as tolerance on residual is modified

clear all; 
close all; 
N = 2^14; 
tols = [1e-2, 1e-3, 1e-4, 1e-6, 1e-8];
signal = 'randn';
dataratio = 1.8; 
M = ceil(dataratio*N); 
condflag = 0; %do not compute condition number
nudist = 3; % random nodes
solvlist = [2:3,5:7,9];

%%
%for j = 1:(length(jitterch)+2)

for j = 1:length(tols)
    if nudist==1, nunam='jittered grid';
      jitter = 0.5;   % size (h units) of jitter from grid; <=0.5 well-cond
      x = (2*pi/M)*((0:M-1) + jitter*(2*rand(1,M)-1));
      w = [];
    elseif nudist==2, nunam='quadrature pts';     % Cheby nodes on [0,2pi]
      [x,w] = clencurt(M-1);
      x = pi*(1+x');   % scale to domain [0,2pi]  (w=quadwei)
    elseif nudist==3, nunam='rand unif iid';
      x = 2*pi*rand(1,M);           % iid unif entire domain
      w = [];
    elseif nudist==4, nunam='rand iid w/ gap';
      sbwp = 8.0;   % gap space-bandwidth prod in half-wavelengths
      x = 2*pi*rand(1,M)*(1-sbwp/N);    % gap w/ no pts in it (exp bad wrt sbwp)
      w = []; 
    end
    opts = []; 
    opts.hsstol = 1e-2*tols(j); 
    opts.cgtol = tols(j); 
    opts.maxit = 10000; 
    opts.nrhs = 1; 
    for k = 1:1
     [relresid, relerror, times, condA] = iterative_test2(solvlist, signal,N,dataratio,x, w, condflag, opts);
    end
    relres(:,j) = relresid(:); 
    relerr(:,j) = relerror(:); 
    timing(:,j) = times(:); 
    condition(j) = condA; 
    j
end

%%
% now do multi-RHS for direct solvers: 
solvlist = 9;
szrhs = [10, 100, 1e3];
relresm = zeros(length(szrhs),length(tols) );
relerrm = relresm;
timingm = relresm;
for k = 1:length(szrhs)
    for j = 1:length(tols)
        if nudist==1, nunam='jittered grid';
            jitter = 0.5;   % size (h units) of jitter from grid; <=0.5 well-cond
            x = (2*pi/M)*((0:M-1) + jitter*(2*rand(1,M)-1));
            w = [];
        elseif nudist==2, nunam='quadrature pts';     % Cheby nodes on [0,2pi]
            [x,w] = clencurt(M-1);
            x = pi*(1+x');   % scale to domain [0,2pi]  (w=quadwei)
        elseif nudist==3, nunam='rand unif iid';
            x = 2*pi*rand(1,M);           % iid unif entire domain
            w = [];
        elseif nudist==4, nunam='rand iid w/ gap';
            sbwp = 8.0;   % gap space-bandwidth prod in half-wavelengths
            x = 2*pi*rand(1,M)*(1-sbwp/N);    % gap w/ no pts in it (exp bad wrt sbwp)
            w = []; 
        end
        opts = []; 
        opts.hsstol = 1e-2*tols(j); 
         %opts.cgtol = tols(j); 
         %opts.maxit = 9000; 
         %opts.nrhs = 1; 
        opts.nrhs = szrhs(k); 
        for i = 1:2
            [relresid, relerror, times, condA] = iterative_test2(solvlist, signal,N,dataratio,x, w, condflag, opts);
        end
        relresm(k,j) = relresid(:); 
        relerrm(k,j) = relerror(:); 
        timingm(k,j) = times(:)/(opts.nrhs);   
        j,k
    end
end
%%
relres = [relres; relresm];
relerr = [relerr; relerrm];
timing = [timing; timingm];

save('iter_rand_decay_214', 'relres', 'relerr', 'timing','condition', 'tols', 'szrhs')
%%


