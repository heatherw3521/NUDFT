function [T, Tinv, gam, M] = mobiusT(I)
%given I = [a b c d] where [a b] [c d] are two disjoint intervals 
% on the real line (or on circle), T(I) maps to the four points [-gamma, -1, 1, gamma]. 
% M is the cross-ratio. 

a = I(1); 
b = I(2); 
c = I(3); 
d = I(4); 

%parameters
M = abs((c-a)*(d-b)/((c-b)*abs(d-a))); 
gam = -1+2*M+2*sqrt(M^2-M); 
A = -gam*a*(1-gam)+gam*(c-gam*d)+c*gam-gam*d; 
B = -gam*a*(c*gam-d)-a*(c*gam-gam*d)-gam*(c*d-gam*d*c); 
C = a*(1-gam)+gam*(c-d)+c*gam-d; 
D = -gam*a*(c-d)-a*(c-gam*d)+c*d-gam*d*c; 

T = @(z) (A*z+B)./(C*z+D);
Tinv = @(z) (D*z-B)./(-C*z+A);
end