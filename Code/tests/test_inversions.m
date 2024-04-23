% test direct solvers only:

% This creates a timing figure for INUDFT, Toeplitz, and KP solvers
% (the direct solvers). 

%%

clear all
close all

ns = 2.^(10:14);
%ns = round(10^(2.5)); 
p = length(ns);

time_nudft = zeros(p,1);
time_Toep = zeros(p,1);
time_KP = zeros(p,1);
time_Gohberg = zeros(p,1);
time_bs = zeros(p,1); 
tol = 1e-10; 
warning off
for k = 1:2
    k
for j = 1:length(ns)
    j
    n = ns(j); 
    m = round(2*n); m = m + mod(m,2); %needs to be even for KP solver?
    %locs = rand(m,1); % initialize locs
    %locs = (1/m)*((0:m-1) + .5*(2*rand(1,m)-1)).'; %jittered nodes
    locs_int = (-0.5:1/m:0.5-1/m) + 1/(4*m)*rand(1,m);
    locs = .5 + locs_int.'; % for KP
    nodes = exp(-2*pi*1i*locs); %for NUDFT solver
    locs_pi = locs*2*pi; 
    rhs = rand(m,1);
    t = tic; 
    out = INUDFT(nodes, n, rhs,'tol', tol);
    tm = toc(t); 
    time_nudft(j) = tm;
    t = tic; 
    out = solv_FDToepN(locs_pi,n,rhs,[]);
    tm = toc(t); 
    time_Toep(j) = tm; 
    t = tic; 
    plan = infft(locs_int, n);
    plan.f = rhs; % Set function values
    infft_trafo(plan); % this is the fast method from KP paper
    tm = toc(t); 
    time_KP(j) = tm;
    t = tic; 
    plan2 = infft(locs_int,n,'flag_toeplitz',1);
    plan2.f = rhs; % Set function values
    infft_trafo(plan2); %this is the Gohberg solver from KP paper
    tm = toc(t);
    time_Gohberg(j) = tm;
    if ~(j == (length(ns)+1))
        t = tic; 
        V = nodes.^(0:n-1); 
        x = V\rhs; 
        tm = toc(t); 
        time_bs(j) = tm;
    end
end

end
%%
% plots
loglog(ns, time_nudft,'.-', 'linewidth', 2.5,'markersize', 20); 
hold on
loglog(ns, time_Toep,'.-', 'linewidth', 2, 'markersize', 20); 
loglog(ns, time_Gohberg, '.-', 'linewidth', 2, 'markersize', 20);
loglog(ns, time_KP, '.-', 'linewidth', 2, 'markersize', 20);
loglog(ns(1:end), time_bs(1:end), '.-', 'linewidth', 2, 'markersize', 20);
shg
legend('nudft', 'toep','gohberg','KP', 'backslash', 'location', 'northwest')

warning on
    
