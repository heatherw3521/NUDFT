function A = circulant(a)
% CIRCULANT  fill dense circulant matrix given its 1st col
%
% A = circulant(a) returns square circulant matrix with first col a.
%  This matches definition in Raymond Chan book, etc.
%
% Without arguments does self test

% Barnett 2/5/08. Transposed convention to match R F Chan book, 11/20/22.
if nargin==0, test_circulant; return; end

a = a(:);
A = toeplitz(a, [a(1); a(end:-1:2)]);   % matlab cmd needs 1st col, 1st row

%%%%%%%%
function test_circulant
N=100;
x = randn(N,1)+1i*randn(N,1);
a = randn(N,1)+1i*randn(N,1);
norm(circulant(a)*x - ifft(fft(a).*fft(x)),'fro')  % fast via periodic conv
