%test cauchytoeplitz

s = 11; 
n = 2^s; %size of matrix

C = build_cauchy(n); 

t = 4; %partition level
pp = n/2^t; %partition index
k = 25;

%%