function [U, J] = fADI_row(sz, I, ridx, cidx, n, tol, G)
% performs ADI on the block C(ridx, cidx) of the square matrix C. 
% using shift parameters p and q. 
% returns U, an approx to the colspace, and indices J, where C approx = U*C(J,:). 
%%

% set up Sylvester operators and get shift parameters: 
a = sz; 
w = exp(pi*1i/n);  
%Dp = spdiags(w.^(2*ridx), 0, a, a); 
Dp = w.^(2*ridx);
 
[p, q] = getshifts_adi(I, 'tol', tol);
pp = -p; 
qq = -q; 
%do fADI on cols. 
[~,r] =size(G); %rank of RHS
k = length(p);  
%Im = speye(size(Dp)); 
%Z(:, 1:r) = (Dp+Im*qq(1))\G;
Z(:,1:r) = G./(Dp+qq(1))*(q(1)-p(1));
%DD =  (q(1)-p(1)).*ones(r,1);

    
for i = 1:k-1 
    %Z(:, i*r+(1:r)) = (Dp+pp(i)*Im)*((Dp+qq(i+1)*Im)\Z(:,(i-1)*r+(1:r)));
    Z(:,i*r+(1:r)) = (Z(:,(i-1)*r+(1:r))./(Dp + qq(i+1))).*(Dp + pp(i))*(q(i+1)-p(i+1)); 
    %DD = [DD; (q(i+1)-p(i+1)).*ones(r,1)];
end

%% col-pivoted qr on ZZ to get ID: 
k = k*r;  
%ZZ = Z; 
%[~, U, P] = qr((ZZ*diag(DD))'); %col-piv QR
[~, U, P] = qr(Z'); %col-piv QR
[v, ~] = find(P'); %vector of indices for LHS
[J, ~] = find(P); %vector of indices for RHS 
U = U';
U = U(:, 1:k); 
U = [eye(k); U(k+1:end, :)/U(1:k, 1:k)];  
U = U(v, :); %permute rows 
J = J(1:k); %indices for ID row selection
end






