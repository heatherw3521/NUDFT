%NUDFT 

%start by testing the diagonalization: 
n = 100; 
F = zeros(n,n); 
w = exp(pi*1i/n); 
for j = 1:n
    for k = 1:n
        F(j,k) = w.^(j*(2*k-1)); 
    end
end
F = F./sqrt(n); 

Ft 

Q = eye(n); 
Q = circshift(Q,1); 

D = F*Q*F';
D(abs(D) < 1e-14) = 0; spy(D)
d_tru = (diag(D)).';
d = w.^((2:2:2*n)); 
max(abs(d.' - d_tru.'))

%%
% using ffts: 

D2 = ft((ft(Q'))');
norm(D2-D)

%%
% using iffts: 
Dt = F'*Q*F; 
Dt2 = ift( ift(Q')');
norm(Dt - Dt2)
%%
% try the Slyvester equation out: 

%nodes: 
nd = linspace(0, 2*pi, n+1); 
nd = nd(1:end-1).'; 
nd = nd + .01*rand(n, 1); 
nd = exp(1i*nd); 


plot(nd, 'x')
hold on
plot(exp(1i*linspace(0, 2*pi, 100)), '.-k');

Dl = diag(nd); 

V = bsxfun(@(x,y) x.^y, nd, 0:n-1); 

rhs = Dl*V - V*Q;

rhs(abs(rhs)< 1e-10) = 0; clf; 
spy(rhs)
%%
%rhs directly: 
a = nd.^n - 1; 
b = zeros(n,1); 
b(end) = 1;
norm(rhs-a*b')

%% try with the transforms: 
C = ft(V')'; %cauchy matrix
L = diag(d_tru);
lft = Dl*C - C*L;
rt = a*(ft(b))';
norm(lft-rt)

%%
%explict (Fb)': 
rt1 = (ft(b))';
rte = w.^(-(1:n)*(2*n-1)); rte = rte/sqrt(n); 
norm(rte - rt1)
shg
 

%% 
% 2D vandermonde test

        
%%
%custom fft, ift that do the scaled transforms we want: 
function out = ft(x)
    % applies Fx, where
    % F_{j,k} = exp(2*pi*1i*j*(2k-1))/sqrt(length(x)). 
    [n, ~] = size(x);
    w = exp(pi*1i/n);
    J = spdiags(sqrt(n)*w.^((1:n)'), 0, n,n); W = spdiags(w.^(2*(1:n)'-2), 0, n, n);
    y = full(W*x); 
    out = full(J*ifft(y)); 
end

function out = ift(x)
    %inverse of ft(x)
    [n,~] = size(x);
    w = exp(pi*1i/n);
    J = spdiags((1/(sqrt(n)))*w.^(-(1:n)'), 0, n,n); W = spdiags(w.^(-2*(1:n)'+2), 0, n, n);
    y = full(J*x);
    out = full(W*fft(y)); 
end







 

