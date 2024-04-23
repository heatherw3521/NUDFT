% check speedup of correct padding in Toep_apply. 1/31/24 Barnett.
N = 2^19;
x = randn(N,1);
t = randn(2*N-1,1);     % define nonsymm Toep: back 1st row then down 1st col
tpad = [t;0]; that = fft(tpad);     % illustrates how to pad
reps = 10;
x2=x;         % save it for 2nd method
tic; for i=1:reps, x = Toep_apply(x,that); end; t1=toc;    % power meth
that = fft(t);          % replace it; illustrates how not to pad (not)
tic; for i=1:reps, x2 = Toep_apply_badpad(that,x2); end; t2=toc;  % note arg swap, ugh
fprintf("Toep_apply speedup factor for N=%d is: %.3g\n",N,t2/t1)
fprintf("(check same answer: %.3g)\n",norm(x-x2)/norm(x))

% results: 1.5x for N=2^19  (we're lucky)
%          9x for N=3e5     (yikes!)
