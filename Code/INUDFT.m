function varargout = INUDFT(nodes,n, b, varargin)
% INUDFT   Fast direct HSS solve Vandermonde form type-2 NUDFT lin sys
%
% x = INUDFT(nodes, n, b) solves the equation Vx=b, where V is a rectangular 
%  NUDFT matrix with n columns: V(j,k) = nodes(j)^(k-1).
%  Evaluation nodes are given by 'nodes', a length-m vector.
%
% x = INUDFT(nodes, n, b, 'tol', tol) solves to the tolerance tol. 
%
% [L, p, x] = INUDFT(nodes, ...) returns L, a urv factorization of a permutation
% of V, and p, the permutation information, as well as the solution x.


% This is a wrapper for @hss/private/hss_cauchy_INUFFT

%%
% we start with a square matrix
%%
% check for tolerance parameter: 
tol = 1e-12; 
if ~isempty(varargin)
    if strcmpi(varargin{1}, 'tol')
        tol = varargin{2}; 
    else
        error('could not parse input')
    end
end

%% if the matrix is small, solve directly: 
% (to do: optimize this choice)
if length(nodes) < 257 && n < 257
    V = (nodes(:)).^(0:n-1); 
    x = V\b;
    varargout = {x};
    return
end

%% Part 1: transform to diagonalized Cauchy system: 

%N = (1:n).'; 
w = exp(pi*1i/n); 
nodes = nodes(:); 

%we want to do ADI on submatrices of H, with DH - HL = ab^T. 

%permute so that nodes are ordered wrt argument:
args = mod(angle(nodes), 2*pi);
[args, p] = sort(args); 
%circshift so nodes are ordered appropriately for subdivisions. 
kk = find(args > pi/n, 1); 
p = circshift(p, -kk+1); 

nodes = nodes(p); % D = diag(nodes)

%L = w.^(2*(1:n)); %L = diag(L)

%set up RHS of Slyester equation: D*C - C*L = a*v
a = nodes.^n - 1; 
v = w.^(-(1:n)*(2*n-1)); v = v/sqrt(n);

%%
% call the hss constructor: 
warning off
H = hss('nudft', nodes, a, v,n,'tol', tol);
%figure()
%spy(H)
%title(['tol is ',num2str(tol)])

%%
%V = nodes.^(0:n-1);
%C = (ft(V'))';

%%
% solve linear system
b = b(p,:);
L = urv(H); 
x = urv_solve(L,b);
if nargout ==1
    varargout = {ift(x)};
else
    x = ift(x);
    varargout = {L,p,x};
end
end


%%%%%%%%%%%%%%%%%%%%
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




