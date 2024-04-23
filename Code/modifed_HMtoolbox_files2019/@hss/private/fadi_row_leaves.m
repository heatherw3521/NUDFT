function x = fadi_row_leaves(G,H,m, N, tol)

%m = block size
%N = matrix size
%GH^* is RHS

w = exp(1i*pi/N); 
[~, r] = size(G); 
I = [w.^(2), w.^(2*m), w.^(2*(m+1)), 1];
%find shift parameters:
[p,q]= getshifts_adi(I1, 'tol', tol);

%set up RHS array: 
G = reshape(G, m, r*N/m); 
%set up colscale: 
s = [w.^(2*m) w.^(2*m)];
j = 0:N/m-1; 
ss = bsxfun(@(s,j) s.^j, s.',j);
S = spdiags(ss, 0, N/m*r); 

%set up LHS ops: 


%do fadi on all RHSs at once: 

Im = speye(size(Dp)); 
Z(:, 1:r) = (Dp-Im*q(1))\(G*S);
ZZ = Z; 
DD =  (q(1)-p(1)).*ones(r,1);
    for i = 1:k-1 
        Z = Z + (Dp - q(i+1)*Im)\((q(i+1)-p(i))*Z*S);
        ZZ = [ZZ Z];
        DD = [DD; (q(i+1)-p(i+1)).*ones(r,1)];
    end


