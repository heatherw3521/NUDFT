% small example of superfast LSQ direct inversion of the type 2 NUDFT.
% Single RHS, with notation as in the paper. Needs only MATLAB. Barnett 4/25/24
addpath Code Code/modifed_HMtoolbox_files2019 Code/extras
n = 2000;       % unknowns, number of regular Fourier frequencies
m = 2*n;        % equations, number of nonuniform nodes
rng(0);          % fix seed
p = rand(m,1);   % iid random nodes in [0,1), needn't be sorted
gamma = exp(2i*pi*p);    % their Vandermonde nodes on unit circle
w = 0:n-1;          % row vec of the frequency indices
xtrue = randn(n,1) + 1i*randn(n,1);  % choose Fourier coeffs at said indices
V = gamma(:).^w;    % dense fill Vandermonde matrix (only for small tests)
b = V*xtrue;        % direct generate RHS (only for small tests)
tic;                % call our solver...
x = INUDFT(gamma,n,b);  % solve one RHS with default tolerance
fprintf('solve in %.3g s: rel l2 resid norm %.3g, rel l2 soln error norm %.3g\n', toc, norm(V*x-b)/norm(b), norm(x-xtrue)/norm(x))
tic;                % compare dense LAPACK generic LSQ solver for same task...
x = V\b;
fprintf('MATLAB dense LSQ solve in %.3g s: rel l2 resid norm %.3g, rel l2 soln error norm %.3g\n', toc, norm(V*x-b)/norm(b), norm(x-xtrue)/norm(x))
