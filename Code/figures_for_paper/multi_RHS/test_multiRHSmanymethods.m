% test for how well different methods do with multiple RHS. 
% In this test we use a jittered grid, which CG likes. 


%%
clear all
close all
n = 2^13; 
m = round(2*n); m = m + mod(m,2); %KP solver needs even m?

locs_int = (-0.5:1/m:0.5-1/m) + 1/(4*m)*rand(1,m); %for KP
locs = .5 + locs_int.'; % for INUDFT
nodes = exp(-2*pi*1i*locs);
xcg = locs*2*pi; 
tol = 1e-8; % tol param for NUDFT

%%
%precomputation NUDFT + quick error check: 
%xt = rand(n,1); 
%V = nodes.^(0:n-1); 
%b = V*xt; 
b = rand(m,1); 
[L,p,x] = INUDFT(nodes,n,b, 'tol', tol);
%norm(x-xt)/norm(xt)

%%
%precomputation KP: %by default seems to do 1e-8, so we will match
% with our tol param, set CG to 1e-6. 
%xt = rand(n,10);
%V = exp(2*pi*1i*locs_int'*(-n/2:n/2-1));
%b = V*xt; 
plan = infft(locs_int, n);
%plan.f = b; 
%infft_trafo(plan)
%norm(plan.fcheck-xt)/norm(xt)

%%
% set up for CG: 
opts.cgtol = 1e-8; 
opts.tol = 1e-12; 
%[c,data] = solv_CGN(xcg,N,b,opts);


%%
tt = 8;
sz = round(logspace(0, 3, tt)); % start with small
timings = zeros(3,tt);

for j = 1:(tt)
    j
    cfs = randn(n,sz(j)) + 1i*randn(n, sz(j)); %coeffs
    B = finufft1d2(xcg,+1,opts.tol,cfs); %generate RHSs
    %INUDFT:
    for k = 1:2
        s = tic;
        X = INUDFT_solve(L, p, B);
        t = toc(s);
    end
    timings(1,j) = t;
    %KP
    for k = 1:2
        s = tic;
        [~,ncol] = size(B);
        for i = 1:ncol
           plan.f = B(:,i); 
           infft_trafo(plan);  
        end
        t = toc(s);
    end
    timings(2,j) = t; 
    %CG
    for k = 1:2
        s = tic; 
        [c,data] = solv_CGN(xcg,n,B,opts);
        t = toc(s);
    end
    timings(3,j) = t; 
end

%%
semilogy(sz, timings.')
hold on
legend('ADI-based', 'KP', 'CG')

%save('dir_multiRHS', 'sz', 'timings', 'n', 'm', 'locs')
    

