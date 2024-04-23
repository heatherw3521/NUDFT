function Ca = circ_apply(vhat,a)
% CIRC_APPLY   Fast mat-vec by circulant mat C given DFT of 1st row of C
%
% Ca = circ_apply(vhat,a) takes vhat as DFT of 1st row of C, and returns C*a
%  fast using the FFT.

% Barnett 11/9/22
Ca = ifft(fft(a).*vhat);   % simple N-periodic convolution
