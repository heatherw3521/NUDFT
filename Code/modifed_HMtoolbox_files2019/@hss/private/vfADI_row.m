function [U, J] = vfADI_row(sz, I, ridx, nodes,n, kk,tol, G)
% performs ADI on the block C(ridx, cidx) of the matrix C = VF'. 
% returns U, an approx to the colspace, and indices J, where C approx = U*C(J,:). 
%%

% set up Sylvester operators and get shift parameters:   
%Dp = spdiags(nodes(ridx), 0, sz, sz); 
Dp  = nodes(ridx); 
%if tol < 1
    %[p, q] = getshifts_adi(I, 'tol', tol);
    %if length(p) > length(ridx)
       % warning('suboptimal block sizes')
       % [p, q] = getshifts_adi(I, length(ridx));
        %if this happens, it means that there was no compression
        %possible on this pass. It shouldn't happen. The above code
        % saves us from breakdown, but the min block-size picker in the
        % tree setup should be adjusted.
    %end
%else
[p, q] = getshifts_adi(I, kk); %kk = numits
%end
pp = -p; 
qq = -q; 
%do fADI on cols. 
%[~,r] =size(G); %rank of RHS
k = length(p);  
%Im = speye(size(Dp)); 
%Z(:, 1) = (Dp+Im*qq(1))\G;
Z(:,1) = (G./(Dp + qq(1)) )*(q(1)-p(1)); 
%DD =  (q(1)-p(1));

    
for i = 1:k-1 
    %Z(:, i+1) = (Dp+pp(i)*Im)*((Dp+qq(i+1)*Im)\Z(:,i)); 
    Z(:,i+1) = ( (Dp+pp(i)).*( Z(:,i)./(Dp + qq(i+1)) ) )*(q(i+1)-p(i+1));                
    %DD = [DD; (q(i+1)-p(i+1))];
end

%% col-pivoted qr on ZZ to get ID:   
%ZZ = Z; 
%[~,U,P] = qr( (ZZ*spdiags(DD,0,k,k))');
[~,U,P] = qr(Z');
%dd = abs(diag(U))> tol*1e-2; k = min(sum(dd),k); 
% DO NOT CHECK diag decay because we want k to match col. choice.

[v, ~] = find(P'); %vector of indices for LHS
[J, ~] = find(P); %vector of indices for RHS 
U = U';
U = U(:, 1:k); 
U = [eye(k); U(k+1:end, :)/U(1:k, 1:k)];  
U = U(v, :); %permute rows 
J = J(1:k); %indices for ID row selection
end






