function [U, J] = vfADI_col(sz, I, cidx,n, tot,tol, H)
% performs ADI on the block C(ridx, cidx) of the matrix C = VF'. 
% returns U, an approx to the rowspace, and indices J, where C approx = C(,J)*U. 
%%

% set up Sylvester operators and get ADI shifts: 
w = exp(pi*1i/n); 
%Dn = spdiags(w.^(2*cidx), 0, sz,sz); 
[p, q] = getshifts_adi(I, tot);
cp = -conj(p); 
cq = -conj(q); 

%do fADI on cols. 
Dn = conj(w.^(2*cidx)); 
k = length(p);  
%In = speye(size(Dn));   
%Y(:,1) = (Dn+cp(1)*In)\H; 
Y(:,1) = H./( Dn + cp(1)); 
DD = (q(1)-p(1)) ; 
    
for j = 1:k-1 
    %Y(:, j+1) = (Dn + cq(j)*In)*((Dn + cp(j+1)*In)\Y(:,j));
    Y(:, j+1) = ( (Dn + cq(j)).*(Y(:,j)./(Dn + cp(j+1))) );
    DD = [DD; (q(j+1)-p(j+1))];
end
%% col-pivoted qr on YY to get ID: 
 % TO DO: add option to use strong RRQR:
 
%YY = Y; 
[~, U, P] = qr(spdiags(DD,0,k,k)*Y'); %col-piv QR
%[~,U,P] = qr(Y'); 
%check for more decay in terms:
dd = abs(diag(U))> tol*1e-2; k = min(sum(dd),k); 
[J, ~] = find(P); %vector of indices for LHS
[z, ~] = find(P'); %vector of indices for RHS    
U = U(1:k, :); 
U = [eye(k) U(1:k, 1:k)\U(:, k+1:end)]; 
U = U(:, z); %permute rows 
J = J(1:k); %indices for ID row selection
U = U'; 
end

%d = diag(U); 
%k = length(d(abs(d) > tol*abs(d(1))));
%k = kk;





