function [U, J] = rand_interpcol(Phi,p)
%creates a row interpolative decomposition 

[~, U, P] = qr(Phi');  %col-piv QR
[J, ~] = find(P); %vector of indices for LHS
[z, ~] = find(P'); %vector of indices for RHS   
U = U(1:p, :); 
U = [eye(p) U(1:p, 1:p)\U(:, p+1:end)]; 
U = U(:, z); %permute rows 
J = J(1:p); %indices for ID row selection
U = U';
end

