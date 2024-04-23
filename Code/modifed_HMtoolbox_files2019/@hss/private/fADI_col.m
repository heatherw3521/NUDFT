function [U, J] = fADI_col(b, I, kk, ridx, cidx, n, tol, H)
% performs ADI on the block C(ridx, cidx) of the square matrix C. 
% using shift parameters p and q. 
% returns U, an approx to the rowspace, and indices J, where C approx = C(,J)*U. 
%%

% set up Sylvester operators and get ADI shifts: 
 
w = exp(pi*1i/n); 
%Dn = spdiags(w.^(2*cidx), 0, b,b); 
Dn = w.^(2*cidx); 
[p, q] = getshifts_adi(I, 'tol', tol);
cp = -conj(p); 
cq = -conj(q); 

%do fADI on cols. 
%Dn = Dn'; 
Dn = conj(Dn);
%[~,r] =size(H);
r = 2; 
k = length(p);  
%In = speye(size(Dn));   
%Y(:,1:r) = (Dn+cp(1)*In)\H; 
Y(:,1:r) = H./(Dn + cp(1)); 
DD =  (q(1)-p(1)).*ones(r,1); 
    
for j = 1:k-1 
    %Y(:, j*r+(1:r)) = (Dn + cq(j)*In)*((Dn + cp(j+1)*In)\Y(:,(j-1)*r+(1:r)));
    Y(:,j*r+(1:r)) = ( ( Y(:,(j-1)*r+(1:r)) )./ (Dn + cp(j+1)) ).*(Dn + cq(j));
    DD = [DD; (q(j+1)-p(j+1)).*ones(r,1)];
end
%% col-pivoted qr on YY to get ID: 
 % TO DO: add option to use strong RRQR:
 
YY = Y; 
%[~, U, P] = qr(diag(DD)*YY'); %col-piv QR
[~, U, P] = qr(spdiags(DD,0, r*k, r*k)*YY');
[J, ~] = find(P); %vector of indices for LHS
[z, ~] = find(P'); %vector of indices for RHS   
k = k*r; 
U = U(1:k, :); 
U = [eye(k) U(1:k, 1:k)\U(:, k+1:end)]; 
U = U(:, z); %permute rows 
J = J(1:k); %indices for ID row selection
U = U'; 
end

%d = diag(U); 
%k = length(d(abs(d) > tol*abs(d(1))));
%k = kk;





