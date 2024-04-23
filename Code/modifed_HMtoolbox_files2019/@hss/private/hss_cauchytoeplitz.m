function H = hss_cauchytoeplitz(N,G, L, varargin)
% this code constructs an HSS factorization
% for the Cauchy-like matrix that comes from 
% applying the (modified) DFT to a Toeplitz matrix.
%
% This code uses the nonsingular displacement formula
% shown in [1]. 
%
% It preserves an interpolative structure so that
% the ULV-like solver in [1] can be used if desired. 
%
% Code is meant to build the blocks for the HSS matrix
% so that it works with the HSS toolbox. 


% TO DO (to further optimize): 
% 1) add ULV-like solver, 
% 2) replace fadi subroutine with mex to code fadi in C
% 3) create option for using Hadamard product + cauchy (cauchy HSS factorization 
% is a precomputation that works for all T). 


%%
% References: [1] Xia, Xi, Gu. "A Superfast Structured Solver
% for Toeplitz Linear Systems via Randomized Sampling."
% SIMAX, Vol. 33 No 3, p.837-858, 2012.
%
% [2] Wilber, Heather. "Computing numerically with rational functions." 
% see Ch. 4. 

%%%
 
block_size = 128; %smallest blocks. %TO DO: adjust this to tol; 
tol = hssoption('threshold'); %rel tol. 
if ~isempty(varargin)
    tolcheck = find(strcmp('tol', varargin));
    if ~isempty(tolcheck)
        tol = varargin{tolcheck+1};  
    end
end


%% Build tree

% Prepare the tree for the HSS structure -- leaving all the blocks empty:
H = hss_build_hss_tree(N, N, block_size);

%build the HSS representation recursively:
H = BuildHSS_iter(H, tol, 0, 0, N, G, L);
end



%% SUBROUTINES

function H = BuildHSS_iter(H, tol, mh, nh, N, G, L)

m = size(H, 1);
n = size(H, 2);

if H.leafnode    
    % Let's do it for the HSS block row first.
    ridx = (mh+1:mh+m).'; 
    cidx = [(1:mh) (mh+m+1:N)].'; 
    w = exp(pi*1i/N); 
    I = [w^(2*(ridx(1)-1)), w^(2*(ridx(end)-1)), w.^(2*(mh+m+1) - 1), w.^(2*mh - 1)]; 

    % use fadi to build an ID:
    a = length(ridx);  
    [H.U, J] = fADI_row(a, I, ridx, cidx, N,tol, G(ridx, :)); %N is global size 
    H.ridx = J + mh; %keep track of row indices. 
   
    % And now do the columns 
    ridx = [(1:nh) (nh+n+1:N)].'; 
    cidx = (nh+1:nh+n).'; 
    %kk = length(J); %heuristic  here for possible future use. 
    b = length(cidx); 
    I = [w.^(2*(nh+n+1-1)), w.^(2*(nh-1)), w.^(2*min(cidx)-1), w.^(2*max(cidx)-1)]; 
    [H.V, J] = fADI_col(b, I, 1,ridx, cidx, N,tol, L(cidx, :)); %N is global size
    H.cidx = J +nh; %keep track of col indices. 
    
    H.D = buildcauchy(mh+1:mh+m, nh+1:nh+n, N, G, L); 
else
    % Call the constructor recursively on the left and right children.
      H.A11 = BuildHSS_iter(H.A11, tol, mh, nh, N, G, L);
    
      [mm, nn] = size(H.A11); 
      mh = mh + mm; 
      nh = nh + nn; 
      H.A22 = BuildHSS_iter(H.A22,tol, mh,nh,N, G, L);
      
      %now build glue matrices for level below:
      H.B12 = buildcauchy(H.A11.ridx, H.A22.cidx, N, G, L); 
      H.B21 = buildcauchy(H.A22.ridx, H.A11.cidx, N, G, L); 
      
      if H.topnode
        return;
      end
      
      %row translation matrices:
      ridx = [H.A11.ridx; H.A22.ridx]; %subset of rows;
       
      idx1 = mh - mm +1; %square start
      idx2 = mh+mm; %square end
      cidx = [(1:idx1-1) (idx2+1:N)].'; 
      a = length(ridx); 
      
      w = exp(1i*pi/N); 
      I = [w^(2*(min(ridx)-1)), w^(2*(max(ridx)-1)), w.^(2*(idx2+1) - 1), w.^(2*(idx1-1) - 1)];
        
      [R, J] = fADI_row(a, I, ridx, cidx, N,tol, G(ridx, :));
      a = length(H.A11.ridx); 
      H.Rl = R(1:a, :); 
      H.Rr = R(a+1:end, :); 
      H.ridx = ridx(J); 
      kk = length(J);
     
      %now col translation matrices: 
      
      ridx = [(1:idx1-1) (idx2+1:N)].';
      cidx = [H.A11.cidx; H.A22.cidx]; 
      b = length(cidx); 
       
      I = [w.^(2*(idx2+1-1)), w.^(2*(idx1-1-1)), w.^(2*min(cidx)-1), w.^(2*max(cidx)-1)]; 
       
      [R, K] = fADI_col(b,I,kk, ridx, cidx, N, tol, L(cidx, :)); 
      H.Wl = R(1:a, :); 
      H.Wr = R(a+1:end, :); 
      H.cidx = cidx(K); 
end
    
end


%%
function AA = buildcauchy(J, K, N, G, L)
%build the cauchy matrix from the given indices:
J = J(:); 
K = K(:); 
w = exp(1i*pi/N); 
a = w.^(2*(J-1));
b = w.^(2*K-1); 
A = bsxfun(@minus, a, b.'); 
A = 1./A; 

%apply hadamard product
[~,r] =size(G); 
k = length(J);
n = length(K); 
AA = 0; 
for j = 1:r
    B = spdiags(G(J, j), 0,k,k)*A*spdiags(conj(L(K,j)), 0, n,n);
    AA = AA +B; 
end


end






