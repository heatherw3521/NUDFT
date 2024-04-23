% test Toeplitz solver.
% (Barnett added complex test 11/18/22)

% simple experiment to see if Toeplitz_solve wrapper works:
clear all
n = 2^10;   % fails for 1e3, ie non-powers-of-2 ?
testcomplex = 1;     % plain solver fails for 1, "ns" version works
tr = rand(n,1);
tc = rand(n,1);
if testcomplex, tr = tr + 1i*rand(n,1); tc = tc + 1i*rand(n,1); end
tc(1) = tr(1); 
T = toeplitz(tc, tr); 
x = rand(n,1);
if testcomplex, x = x + 1i*rand(n,1); end
b = T*x; 
tol = 1e-10; 

%%
t = tic; 
xs = Toeplitz_solve(tc,tr,b, 'tol', tol);
t = toc(t); 
norm(x -xs)/norm(x)
norm(x -real(xs))/norm(x)     % why? seems to show complex tc,tr ok, not x ?

%%
% try the nonsingular version: 
xns = Toeplitz_solve_ns(tc,tr,b, 'tol', tol);

norm(x-xns)/norm(x)
norm(x-real(xns))/norm(x)


%%
% try the hadamard product version: STILL IN PROG
% NOTE: something is very wrong with this one!
xh = Toeplitz_solve_hdm(tc,tr, b, 'tol', tol);

norm(x-real(xh))/norm(x)
