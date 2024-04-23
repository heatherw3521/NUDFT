function [U, J] = fADI_col1(C, ridx, cidx, n, tol)
% performs ADI on the block C(ridx, cidx) of the square matrix C. 
% using shift parameters p and q. 
% returns U, an approx to the rowspace, and indices J, where C approx = C(,J)*U. 
%%

% set up Sylvester operators: 
[~, b] = size(C); 

w = exp(pi*1i/n); 

Dn = diag(w.^(2*cidx - 1));
V = ones(b,1);  
I = [w.^(2*(ridx(1)-1)), w.^(2*(ridx(end)-1)), w.^(2*min(cidx)-1), w.^(2*max(cidx)-1)]; 
[p, q] = getshifts_adi(I, 'tol', tol);
cp = -conj(p);
cq = -conj(q);
%do fADI on rows. 
Dn = Dn'; 
r = 1; 
k = length(p);  
In = speye(size(Dn)); 
Y(:,1:r) = (Dn-cp(1)*In)\V;   
%YY = Y; 
DD =  (q(1)-p(1)).*ones(r,1); 
    %for i = 1:k-1 
       % Y = Y + (Dn-cp(i+1)*In)\((cp(i+1)-cq(i))*Y);
       % YY = [YY Y];
       % DD = [DD; (q(i+1)-p(i+1)).*ones(r,1)];
    %end
    for j = 1:k-1 
        Y(:, j*r+(1:r)) = (Dn + cq(j)*In)*((Dn + cp(j+1)*In)\Y(:,(j-1)*r+(1:r)));
        DD = [DD; (q(j+1)-p(j+1)).*ones(r,1)];
    end
%% col-pivoted qr on YY to get ID: 
 % TO DO: replace with strong RRQR
YY = Y; 
[~, U, P] = qr(diag(DD)*YY'); %col-piv QR
[J, ~] = find(P); %vector of indices for LHS
[z, ~] = find(P'); %vector of indices for RHS   
U = [eye(k) U(1:k, 1:k)\U(:, k+1:end)]; 
U = U(:, z); %permute rows 
J = J(1:k); %indices for ID row selection
U = U'; 
end


% [Z, D, Y, ~, ~] = HSSblock_fADI(C, ridx, cidx, k);
% Y = D*Y';
% [Q, RR, P] = qr(Y); 
% [v, ~] = find(P); %left side idx
% [z, ~] = find(P'); % right side idx
% R1 = RR(1:k, 1:k); 
% 
% M = Z*Q*R1; 
% R = [eye(k) R1\RR(:, k+1:end)]; 
% R = R(:,z); 
% v = v(1:k); 





