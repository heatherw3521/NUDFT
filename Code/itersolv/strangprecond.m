function [v vhat] = strangprecond(t)
% STRANGPRECOND  Strang circulant preconditioner for square nonsymm Toeplitz
%
% [v vhat] = strangprecond(t) assumes t is vector length 2N-1 defining
%  Toeplitz (Chan convention, as Toep_apply: reversed 1st row then down the
%  rest of 1st col), and returns length-N col vec v defining circulant
%  matrix C^{-1} so that C^{-1}T is supposedly better-conditioned. The
%  DFT vhat = fft(v) is also returned so that the fast matvec by C^{-1}
%  can be done without an extra FFT.
%  Strang's rule is that C is given by circulant-izing T by taking central
%  band of width N/2.
%
% With no arguments, self-test is done.
%
% See also: CIRCULANT, CIRC_APPLY, TOEP_APPLY

% Barnett 11/9/22. Preserved Hermitian using Re of nsplit entry 11/12/22.
if nargin==0, test_strangprecond; return; end

N = floor((numel(t)+1)/2);
assert(numel(t)==2*N-1)
nsplit = floor(N/2);
t = t(:).';   % row
Cvec = [t(N:N+nsplit-1), real(t(N+nsplit)) t(N+1-nsplit:N-1)].';   % len N, col
vhat = 1./fft(Cvec);         % inv of circulant = inv its eigvals
v = ifft(vhat);

%%%%%%%
function test_strangprecond
verb = 0;
N=100;    % size, even (and may affect whether SPD below)
rng(0);

if 0
  disp('dumb direct tests when Toeplitz is circulant (should be ~eps_mach)...')
  Cvec = randn(1,N) + 1i*randn(1,N);
  v = ifft(1./fft(Cvec));     % check basic math, its a left-inv of C
  fprintf('err inverting a C: %.3g\n',norm(circulant(v)*circulant(Cvec) - eye(N),'fro'))
  t = [Cvec(2:N) Cvec];         % Toep defn vec
  [v,vhat] = strangprecond(t);  % test our func
  fprintf('Strang err with C: %.3g\n',norm(circulant(v)*circulant(Cvec) - eye(N),'fro'))
end

disp('Toeplitz in Wiener class... but smoothness controls cond#...')
%t = exp(-0.1*(-N+1:N-1).^2);      % gaussian smooth (-> ill cond), rapid decay
t = exp(-0.38*abs(-N+1:N-1));      % kink (-> well cond), decay, SS rank=1 :)
%t = (5./(5+abs(-N+1:N-1))).^5;    % rapid decay, symm
%t = (5./sqrt(5^2+(-N+1:N-1).^2)).^5;    % rapid decay, symm
%t = t.*rand(1,2*N-1); t(N)=1;      % make nonsymm and non-SS. Need SPD for CG!
t = t.*(rand(1,2*N-1) + 1i*rand(1,2*N-1)); t(N)=1;      % make nonsymm and non-SS. Need SPD for CG!
t = t + conj(t(end:-1:1));             % SPD, including Hermitian case
[v,vhat] = strangprecond(t);
%t, v, Cvec = ifft(1./vhat);
T = toeplitz(t(N:end),t(N:-1:1));
%figure; imagesc(abs(T-circulant(Cvec))); colorbar
Ci = circulant(v);
fprintf('||T-T^*||_F = %.3g, \t||C^{-1}-C^{-*}||_F = %.3g\n',norm(T-T','fro'),norm(Ci-Ci','fro'))
fprintf('kappa(T)=%.3g (min eig(T) = %.3g)\n',cond(T),min(eig(T)))
fprintf('kappa(C^{-1}T)=%.3g\n',cond(Ci*T))   % not a great precond!
fprintf('||T||_F = %.3g\n',norm(T,'fro'))
fprintf('||C^{-1}T - I||_F = %.3g\n',norm(Ci*T - eye(N),'fro'))
if verb, figure(1); clf; plot(complex(eig(T)),'+'); hold on;
  plot(complex(eig(Ci*T)),'+'); axis equal;
  legend('spec T','spec C^{-1}T'); title('test Strang precond'); end

if verb                          % show matrices
  figure; subplot(1,3,1);
  imagesc(log10(abs(T))); caxis([-15 0]); colorbar; title('T');
  axis equal tight; subplot(1,3,2);
  imagesc(log10(abs(inv(Ci)))); caxis([-15 0]); colorbar;
  title('C'); axis equal tight; subplot(1,3,3);
  imagesc(log10(abs(Ci*T-eye(N)))); caxis([-15 0]); colorbar;
  title('C^{-1}T - I'); axis equal tight
end

disp('iterative conv rate tests...'); maxit=N; tol=1e-12;
if norm(T-T','fro')>1e-14 || min(eig(T))<=0, warning('T not SPD; CG will fail!'); end
b = randn(N,1) + 1i*randn(N,1);
c0 = T\b; fprintf('dense solve soln nrm %.3g\n',norm(c0))  % exact
[c,flag,relres,iter,resvec] = pcg(T,b,tol,maxit);
fprintf('plain dense CG flag=%d: %d its\trelres=%.3g, soln err rel nrm %.3g\n',flag,iter,relres,norm(c-c0)/norm(c0))
[c,flag,relres,iter,resvec] = pcg(@(a) Toep_apply(a,fft([t,0])),b,tol,maxit);
fprintf('plain fast  CG flag=%d: %d its\trelres=%.3g, soln err rel nrm %.3g\n',flag,iter,relres,norm(c-c0)/norm(c0))
[c,flag,relres,iter,resvec] = pcg(T,b,tol,maxit,@(x) circulant(v)*x);
fprintf('Strang dens CG flag=%d: %d its\trelres=%.3g, soln err rel nrm %.3g\n',flag,iter,relres,norm(c-c0)/norm(c0))
[c,flag,relres,iter,resvec] = pcg(@(a) Toep_apply(a,fft([t,0])),b,tol,maxit,@(x) circ_apply(vhat,x));
fprintf('Strang fast CG flag=%d: %d its\trelres=%.3g, soln err rel nrm %.3g\n',flag,iter,relres,norm(c-c0)/norm(c0))
