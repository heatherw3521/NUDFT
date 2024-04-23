% test multiple RHS

clear all
%%
n  = 2^15; m = 2*n;   % n=#unknowns, m=#eqns.    Eg n=2^18 takes 100+60s 8-core
nd = exp(-1i*2*pi*rand(1,m)).';    % unit circle convention. iid rand
%%
N = 10000;
B = rand(m,N); %rhs 
tol=1e-6;

%s = tic;
%ds one RHS?)
[L,p,x] = INUDFT(nd,n, B(:,1), 'tol', tol);
%fprintf('%d-by-%d @ tol=%g, factor: %.3g s\n',m,n,tol,toc(s))
%%
s = tic;
X = INUDFT_solve(L, p, B);       % should be BLAS3 as of 9/14/23
t = toc(s);
fprintf('solve: %.3g s for %d RHS (%.3g s per RHS)\n',t,N,t/N)

if n*m<1e7
  fprintf('please wait for dense solve...\n')
  V = nd.^(0:n-1);    % do it densely
  Xt = V\B;           % "
  fprintf('err rel to dense solve: %.3g\n', norm(X- Xt)./norm(Xt))
end

% if 0 % compare sequential RHS solver for kicks...
%   tic;
%   [~,s] = size(B);
%   X = zeros(L.n, s);
%   for j=1:10
%     X(:,j) = urv_solve(L, B(:,j));
%   end
%   toc
% end
