function Ta = Toep_apply_badpad(vhat,a)
% TOEP_APPLY_BADPAD   fast matvec with square Toeplitz matrix (obsolete ver)
%
% Ta = Toep_apply_badpad(vhat,a) multiplies vector a by the square (generally
%  non-symmetric) Toeplitz matrix T defined by a vector v, whose DFT
%  vhat = fft(v) the user must supply. The convention for v (as in Raymond
%  Chan's book) is the 1st row of T in reverse order followed by the 2nd through
%  last elements of the 1st column in usual order. In the literature v is
%  indexed -m+1:m-1, where m is the matrix size. T*a is a discrete nonperiodic
%  convolution, and performed here by FFT & iFFT pair.
%
% Inputs: vhat = DFT of v (length 2m-1)
%         a = input column vector length m
% Output: Ta = T*a, col vec length m
%
% Without arguments does self-test
%
% NOTE: OBSOLETE; SEE TOEP_APPLY

% Barnett 11/7/22
if nargin==0, test_Toep_apply; return; end

m = numel(a);
assert(numel(vhat)==2*m-1)
apadhat = fft(a(:),2*m-1);   % zero-pads out to size of vhat
Ta = ifft(apadhat .* vhat(:));
Ta = Ta(m:end);              % extract correct indices

%%%%%%%
function test_Toep_apply
N=10;                   % compare fast against direct method
rng(0)
a = randn(N,1);
t = randn(2*N-1,1);       % define nonsymm Toep: back 1st row then down 1st col
T = toeplitz(t(N:end),t(N:-1:1));   % munge single toep vec into (C,R) format
Ta = Toep_apply_badpad(fft(t),a);
fprintf('Frob norm of diff btw fast and direct: %.3g\n',norm(T*a - Ta,'fro'))
