function H = hss_cauchytoeplitzc(N,varargin)
% this code constructs an HSS factorization
% for the Cauchy matrix that comes from 
% applying the DFT to a Toeplitz matrix.
%
% 
%
% note: this uses the nonsingular disp structure and
% builds an HSS approximation for just the Cauchy matrix. 
% for a generic Toeplitz, this is called by hss_cauchy_toep_hadamard
% as a precomputation. Then, a hadamard product is applied to get 
% the transformed Toeplitz matrix in HSS form. 

 
block_size = 128; %smallest blocks. %need to adjust this to tol; 
tol = hssoption('threshold'); %rel tol. 
if ~isempty(varargin)
    tolcheck = find(strcmp('tol', varargin));
    if ~isempty(tolcheck)
        tol = varargin{tolcheck+1};  
    end
end


%% PART 2: Build tree


% Prepare the tree for the HSS structure -- leaving all the blocks empty
H = hss_build_hss_tree(N, N, block_size);

H = BuildHSS_iter(H, tol, 0, 0, N);
end



%% SUBROUTINES
% 
function H = BuildHSS_iter(H, tol, mh, nh, N)

m = size(H, 1);
n = size(H, 2);

if H.leafnode    
    % Let's do it for the HSS block row first.
    % call the first block:
    ridx = mh+1:mh+m; 
    cidx = nh+n+1:N; 
    B = buildcauchy(ridx, cidx, N); 
    %% use fadi to build an ID:
  
    [H.U, J] = fADI_row1(B, ridx, cidx, N,tol); %N is global size 
    H.ridx = J; %keep track of row indices. 
   
    % And now do the columns
    ridx = mh+m+1:N; 
    cidx = nh+1:nh+n; 
    B = buildcauchy(ridx, cidx, N); 
    [H.V, J] = fADI_col1(B, ridx, cidx, N,tol); %N is global size
    H.cidx = J; %keep track of col indices. 
    
    H.D = buildcauchy(mh+1:mh+m, nh+1:nh+n, N); 
else
    % Call the constructor recursively on the left child then build the
    % right child from the left one.  
      H.A11 = BuildHSS_iter(H.A11, tol, mh, nh, N);
    
      [mh, nh] = size(H.A11); 
      H.A22 = shift_build(H.A11,mh,nh,N); 
      
      %now build glue matrices for level below:
      H.B12 = buildcauchy(H.A11.ridx, H.A22.cidx, N); 
      H.B21 = buildcauchy(H.A22.ridx, H.A11.cidx, N); 
      
      if H.topnode
        return;
      end
      
      %row transition matrices:
      ridx = [H.A11.ridx; H.A22.ridx]; %subset of rows;
      cidx = 2*nh+1:N; 
      B = buildcauchy(ridx, cidx, N); 
      [R, J] = fADI_row1(B, ridx, cidx, N,tol);
      a = length(H.A11.ridx); 
      H.Rl = R(1:a, :); 
      H.Rr = R(a+1:end, :); 
      H.ridx = ridx(J); 
     
      %now col transition matrices: 
      ridx = 2*mh+1:N; 
      cidx = [H.A11.cidx; H.A22.cidx]; 
      B = buildcauchy(ridx, cidx, N); 
      [R, K] = fADI_col1(B, ridx, cidx, N, tol); 
      H.Wl = R(1:a, :); 
      H.Wr = R(a+1:end, :); 
      H.cidx = cidx(K); 
        
end
    

    
 
end


%%
 function Y = shift_build(X, mh, nh, N)
      %build a right child by shifting indices. 
      
      Y = X; 
      [m, n] = size(X); 
      if Y.leafnode %assign elements
          Y.ridx = Y.ridx + m; 
          Y.cidx = Y.cidx + n; 
          %Y.D = A(mh+1:mh+m, nh+1:nh+n); 
          Y.D = buildcauchy(mh+1:mh+m, nh+1:nh+n, N); 
      else
          %transition matrices match X; nothing to do.
          
          %shift IDX
          Y.ridx = Y.ridx + m; 
          Y.cidx = Y.cidx + n; 
          
          %drill down recursively and build children
          Y.A11 = shift_build(X.A22, mh, nh, N);
          Y.A22 = shift_build(Y.A11, mh+size(Y.A11, 1), nh+size(Y.A11, 2), N); 
          
          %build glue for level below: 
          Y.B12 = buildcauchy(Y.A11.ridx, Y.A22.cidx, N); 
          Y.B21 = buildcauchy(Y.A22.ridx, Y.A11.cidx, N); 
          
      end
 end

%%
function A = buildcauchy(J, K, N)
%build the cauchy matrix from the given indices:
J = J(:); 
K = K(:); 
w = exp(1i*pi/N); 
a = w.^(2*(J-1));
b = w.^(2*K-1); 
A = bsxfun(@minus, a, b.'); 
A = 1./A; 
end






