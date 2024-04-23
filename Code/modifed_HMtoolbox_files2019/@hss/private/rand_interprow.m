function [U, J] = rand_interprow(Phi,p)
%creates a row interpolative decomposition 
%Phi = Y(ridx, :) - D*X(ridx,:); 
[~, U, P] = qr((Phi)'); %col-piv QR
[v, ~] = find(P'); %vector of indices for LHS
[J, ~] = find(P); %vector of indices for RHS 
U = U';
%d = diag(U); 
%k = length(d(abs(d) > 10*eps)); 
U = U(:, 1:p); 
U = [eye(p); U(p+1:end, :)/U(1:p, 1:p)];  
U = U(v, :); %permute rows 
J = J(1:p);
end
