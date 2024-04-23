function [err, time] = test_inudft_basics(m,n)
%test_inudft rectangular case: 
% these are just some basic tests for code functionality. 
% returns acc, time to solve for three different node sets. 
%
% Note that if m, n are chosen to be very large this code will be slow 
% because it builds and solves the full system naively to test accuracy. 

n  = 2^9; m = 2*n; cc = clencurt(m); cc= cc(1:end-1); cc = (cc +1)/2; 
nd_choice = {,exp(-2i*pi*(cc.')).',  exp(-1i*2*pi*rand(1,m)).'};
%equispaced, quadrature pts, random

for j = 1:length(nd_choice)
    nd = nd_choice{j};
    b = rand(m,1); %rhs 
    V = nd.^(0:n-1);
    vt = V\b;
    s = tic; 
    v = INUDFT(nd,n,b);
    time(j) = toc(s);
    err(j) = norm(v - vt)./norm(vt);
end

end
