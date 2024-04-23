%test iterative methods: 
% use AB's test suite to test direct method at lower tol than 
% the tol =1e-12 in iterative_four_node_choices.m

clear all; 
close all; 
N = 2^18; 

%jitterch = [.25, .5, 1, 10,50,100];
signal = 'randn';
dataratio = 2; 
M = dataratio*N; 
condflag = 0; %do not compute condition number


%%
%for j = 1:(length(jitterch)+2)

for j = 1:4
%for j = 1:1 
    nudist = j;  % loop over NU pt distns for x (easy to hard)..........
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
    for k = 1:2
     [relresid, relerror, times, numits, condA] = iterative_test(signal,N,dataratio,x, w, condflag);
    end
    relres(:,j) = relresid(:); 
    relerr(:,j) = relerror(:); 
    timing(:,j) = times(:); 
    numiter(:,j) = numits(:); 
    condition(j) = condA; 
    j
end


%%


