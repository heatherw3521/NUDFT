% this file tests some ideas for solving the NUDFT with 2D vandermonde matrices 
% 
%%
% product of C'*C, where C is Cauchy, will be low rank: 
n = 2^11;
m = 2*n; 
x = linspace(0, 2*pi, m).'; 
y = linspace(0, 2*pi, n); y = y+.2;  
C = 1./(x-y); 
N = C.'*C; 
%Nbul = N(1:250, 251:500); %low rank
%NBlr = N(251:500, 1:250);
H = hm(N); %compress off diag blocks;
spy(H)
%%
% structure of the 2D vandermonde matrix



