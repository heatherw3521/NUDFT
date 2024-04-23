function H = hss_cauchytoeplitz2(N,G, L, dg, varargin)
% this code constructs an HSS factorization
% for the Cauchy matrix that comes from 
% applying the DFT to a Toeplitz matrix.
%
% It preserves an interpolative structure so that
% the ULV-like solver in [1] can be used. 
%
% This version of the code uses the singular disp eqn. 
%
 
block_size = 128; %smallest blocks. %need to adjust this to tol; 
tol = hssoption('threshold'); %rel tol. 
if ~isempty(varargin)
    tolcheck = find(strcmp('tol', varargin));
    if ~isempty(tolcheck)
        tol = varargin{tolcheck+1};  
    end
end


%% For the leaves, we only need one U, V pair: 
%find u:
ridx = 1:block_size; 
cidx = (block_size+1):N; 
w = exp(1i*pi/N);
I = [w^(2*(ridx(1))), w^(2*(ridx(end))), w.^(2*(cidx(1))), w.^(2*(cidx(end)))];
a = block_size; 
[Ua, ~] = fADI_row(a, I,ridx.', cidx.', N,tol, ones(block_size, 1));
I = [w^(2*(cidx(1))), w^(2*(cidx(end))), w.^(2*(ridx(1))), w.^(2*(ridx(end)))];
[Va, ~] = fADI_col(a, I,1, cidx.', ridx.', N, tol, ones(block_size, 1)); 
%% PART 2: Build tree


% Prepare the tree for the HSS structure -- leaving all the blocks empty
H = hss_build_hss_tree(N, N, block_size);

H = BuildHSS_iter(H, tol, 0, 0, N, G, L, Ua, Va, dg);
end



%% SUBROUTINES
% 
function H = BuildHSS_iter(H, tol, mh, nh, N, G, L, Ua, Va, dg)

m = size(H, 1);
n = size(H, 2);

if H.leafnode    
    % Let's do it for the HSS block row first.
    ridx = (mh+1:mh+m).'; 
    %cidx = [(1:mh) (mh+m+1:N)].'; 
    %w = exp(pi*1i/N); 
    %I = [w^(2*(ridx(1)-1)), w^(2*(ridx(end)-1)), w.^(2*(mh+m+1) - 1), w.^(2*mh - 1)];
    %I = [w^(2*(ridx(1))), w^(2*(ridx(end))), w.^(2*(mh+m+1)), w.^(2*mh)];
    %cidx = (nh+n+1:N).'; 

    %% use fadi to build an ID:
    %a = length(ridx);  
    %[H.U, J] = fADI_row(a, I, ridx, cidx, N,tol, G(ridx, :)); %N is global size 
    
     
    [~,k] = size(Ua); 
     %U = [spdiags(G(ridx, 1), 0, m, m)*Ua, spdiags(G(ridx, 2),0, m, m)*Ua];
    U = [diag(G(ridx, 1))*Ua, diag(G(ridx, 2))*Ua];
     [H.U, J] = rand_interprow(U,2*k);  
    H.ridx = J + mh; %keep track of row indices.
    % And now do the columns
    %ridx = (mh+m+1:N).'; 
    %ridx = [(1:nh) (nh+n+1:N)].'; 
    cidx = (nh+1:nh+n).'; 
    %V = [spdiags(L(cidx, 1), 0, n, n)*Va, spdiags(L(cidx, 2), 0, n, n)*Va]; 
    V = [diag(L(cidx, 1))*Va, diag(L(cidx, 2))*Va];
    [H.V, K] = rand_interpcol(V, 2*k);  
    %b = length(cidx); 
    %I = [w.^(2*(nh+n+1-1)), w.^(2*(nh-1)), w.^(2*min(cidx)-1), w.^(2*max(cidx)-1)]; 
    %I = [w.^(2*(nh+n+1)), w.^(2*(nh)), w.^(2*min(cidx)), w.^(2*max(cidx))]; 
    %[H.V, J] = fADI_col(b, I, kk,ridx, cidx, N,tol, L(cidx, :)); %N is global size
    H.cidx = K +nh; %keep track of col indices. 
    
    H.D = buildcauchy(mh+1:mh+m, nh+1:nh+n, N, G, L, dg); 
else
    % Call the constructor recursively on the left and right children.
      H.A11 = BuildHSS_iter(H.A11, tol, mh, nh, N, G, L, Ua, Va, dg);
    
      [mm, nn] = size(H.A11); 
      %H.A22 = shift_build(H.A11,mh,nh,N); 
      mh = mh + mm; 
      nh = nh + nn; 
      H.A22 = BuildHSS_iter(H.A22,tol, mh,nh,N, G, L, Ua, Va, dg);
      %now build glue matrices for level below:
      H.B12 = buildcauchy(H.A11.ridx, H.A22.cidx, N, G, L, []); 
      H.B21 = buildcauchy(H.A22.ridx, H.A11.cidx, N, G, L,[] ); 
      
      if H.topnode
        return;
      end
      
      %row translation matrices:
      ridx = [H.A11.ridx; H.A22.ridx]; %subset of rows;
      %cidx = 2*nh+1:N; 
      idx1 = mh - mm +1; %square start
      idx2 = mh+mm; %square end
      cidx = [(1:idx1-1) (idx2+1:N)].'; 
      a = length(ridx); 
      %B = buildcauchy(ridx, cidx, N); 
      %I = [w^(2*(min(ridx)-1)), w^(2*(max(ridx)-1)), w.^(2*(mh+m+1) - 1), w.^(2*mh - 1)];
      w = exp(1i*pi/N); 
      %I = [w^(2*(min(ridx)-1)), w^(2*(max(ridx)-1)), w.^(2*(idx2+1) - 1), w.^(2*(idx1-1) - 1)];
       I = [w^(2*(min(ridx))), w^(2*(max(ridx))), w.^(2*(idx2+1)), w.^(2*(idx1-1))];  
      [R, J] = fADI_row(a, I, ridx, cidx, N,tol, G(ridx, :));
      a = length(H.A11.ridx); 
      H.Rl = R(1:a, :); 
      H.Rr = R(a+1:end, :); 
      H.ridx = ridx(J); 
      kk = length(J);
     
      %now col translation matrices: 
      %ridx = 2*mh+1:N; 
      ridx = [(1:idx1-1) (idx2+1:N)].';
      cidx = [H.A11.cidx; H.A22.cidx]; 
      b = length(cidx); 
      %B = buildcauchy(ridx, cidx, N); 
      %I = [w.^(2*(idx2+1-1)), w.^(2*(idx1-1-1)), w.^(2*min(cidx)-1), w.^(2*max(cidx)-1)]; 
      I = [w.^(2*(idx2+1)), w.^(2*(idx1-1)), w.^(2*min(cidx)), w.^(2*max(cidx))]; 
      [R, K] = fADI_col(b,I,kk, ridx, cidx, N, tol, L(cidx, :)); 
      H.Wl = R(1:a, :); 
      H.Wr = R(a+1:end, :); 
      H.cidx = cidx(K); 
      
        
end
    

    
 
end


%%
function AA = buildcauchy(J, K, N, G, L, dg)
%build the cauchy matrix from the given indices:
J = J(:); 
K = K(:); 
w = exp(1i*pi/N); 

a = w.^(2*J); 
b = w.^(2*K); 
%a = w.^(2*(J-1));
%b = w.^(2*K-1); 
A = bsxfun(@minus, a, b.'); 
A = 1./A; 
A(isinf(A)) = 0; 
  

%apply hadamard product
[~,r] =size(G); 
k = length(J);
n = length(K); 
AA = 0; 
for j = 1:r
    B = spdiags(G(J, j), 0,k,k)*A*spdiags(conj(L(K,j)), 0, n,n);
    AA = AA +B; 
end

if ~isempty(dg) %should be J == K if this is called.
    AA = AA + diag(dg(J));
end
  

end






