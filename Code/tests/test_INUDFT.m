%test_inudft: 
clear all
close all
%%
n  = 2^11;
m = 2*n; 
% nd = linspace(0, 2*pi, m+1); 
% nd = nd(1:end-1).'; 
% nd = nd + (nd(2)-nd(1))*(.5)*rand(m, 1); 
% nd = exp(-1i*nd);
nd = exp(-2i*pi*rand(m,1));
%x = pi*(1+clencurt(M-1)');
%x = clencurt(m); x = x(1:end-1);
%nd = exp(-1i*pi*(1+x)); %gets CC pts but leaves off redundant end



%clf
%plot(nd, 'x')
%hold on
%plot(exp(1i*linspace(0, 2*pi, 100)), '.-k');

%%
b = rand(m,1); %rhs 
V = nd.^(0:n-1);
vt = V\b;

v = INUDFT(nd,n,b);
norm(v - vt)/norm(vt)
