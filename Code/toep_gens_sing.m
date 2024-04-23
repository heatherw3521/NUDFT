function [G, H] = toep_gens_sing(tr, tc) 
%given a toeplitz matrix with 
% rowvec tr, colvec tc, 
% this constructs generators GH^*
% so that Q_1T - TQ_1 = GH^*
tr = tr(:); 
tc = tc(:); 
n = length(tr); 

r1 = flip(tc(2:end)) - tr(2:end);
r1 = [r1; 0];

 
r2 = flip(tr(2:end)) - tc(2:end); 
r2 = [0; r2]; 
eye = [1; zeros(n-1, 1)];
G = [ eye r2]; 
H = [r1 flip(eye)]; 
H = conj(H); % so that GH^* gets the right answer
end