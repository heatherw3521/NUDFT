function [U, J] = fADI_row1(C, ridx, cidx, n, tol)
% performs ADI on the block C(ridx, cidx) of the square matrix C. 
% using shift parameters p and q. 
% returns U, an approx to the colspace, and indices J, where C approx = U*C(J,:). 
%%

% set up Sylvester operators: 
[a, ~] = size(C); 

w = exp(pi*1i/n);  
Dp = diag(w.^(2*(ridx-1))); 
U = ones(a, 1); 

I = [w^(2*(min(ridx)-1)), w^(2*(max(ridx)-1)), w.^(2*cidx(1) - 1), w.^(2*cidx(end) - 1)]; 
[p, q] = getshifts_adi(I, 'tol', tol);
pp = -p; 
qq = -q;
%do fADI on cols. 
r =1; %rank of RHS
k = length(p);  
Im = speye(size(Dp)); 
Z(:, 1:r) = (Dp-Im*q(1))\U;
%ZZ = Z; 
DD =  (q(1)-p(1)).*ones(r,1);
    %for i = 1:k-1 
       % Z = Z + (Dp - q(i+1)*Im)\((q(i+1)-p(i))*Z);
       % ZZ = [ZZ Z];
       % DD = [DD; (q(i+1)-p(i+1)).*ones(r,1)];
    %end
    for i = 1:k-1 
        %Z = Z + (Dp - q(i+1)*Im)\((q(i+1)-p(i))*Z);
        %Z = Z + (Dp + qq(i+1)*Im)\((q(i+1)-p(i))*Z);
        %ZZ = [ZZ Z];
        Z(:, i*r+(1:r)) = (Dp+pp(i)*Im)*((Dp+qq(i+1)*Im)\Z(:,(i-1)*r+(1:r))); 
        DD = [DD; (q(i+1)-p(i+1)).*ones(r,1)];
    end

%% col-pivoted qr on ZZ to get ID: 
ZZ = Z;  
[~, U, P] = qr((ZZ*diag(DD))'); %col-piv QR
[v, ~] = find(P'); %vector of indices for LHS
[J, ~] = find(P); %vector of indices for RHS 
U = U';  
U = [eye(k); U(k+1:end, :)/U(1:k, 1:k)];  
U = U(v, :); %permute rows 
J = J(1:k); %indices for ID row selection
end






