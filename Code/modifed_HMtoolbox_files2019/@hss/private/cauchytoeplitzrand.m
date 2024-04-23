function H = cauchytoeplitzrand(N,G, L,t0, t1, varargin)
% this code constructs an HSS factorization
% for the Cauchy matrix that comes from 
% applying the DFT to a Toeplitz matrix. It uses
% randomized LA to do the compressions. 
% t0 = toeplitz row vec
% t1 = toeplitz col vec
%
% This method is based on [1]. 
%
% It preserves an interpolative structure so that
% the ULV-like solver in [1] can be used. 
%
 
block_size = 128; %smallest blocks. %need to adjust this to tol; 
tol = hssoption('threshold'); %rel tol. 
if ~isempty(varargin)
    tolcheck = find(strcmp('tol', varargin));
    if ~isempty(tolcheck)
        tol = varargin{tolcheck+1};  
    end
end

%% PART 1: Get approx colspace/rowspace
%determine rank p and oversampling parameter q: 
%to do this we just take the finest block: 
% and use Zolotarev bounds:

mu = 8; %oversampling parameter
w = exp(pi*1i/N);  
ridx = 1:block_size; 
cidx = block_size+1:N; 
[~,r] = size(G); %rank of RHS
I = [w^(2*(min(ridx)-1)), w^(2*(max(ridx)-1)) , w.^(2*cidx(1) - 1), w.^(2*cidx(end) - 1)];

% compute cross-ratio:
cr = get_cr(I); 

%bound on rank:
p = r*ceil(1/pi^2*log(4/tol)*log(16*cr));
p = p + mu; 

% now compute approx colspan/rowspan: 
X = rand(N, p); %draw random matrix

%colspan = CX
Y = ct_mul(t0,t1, X); 
%rowspan = C^TX 
Z = ct_mul_tp(t0,t1, X);

%% PART 2: Build tree
% Prepare the tree for the HSS structure -- leaving all the blocks empty
H = hss_build_hss_tree(N, N, block_size);

H = BuildHSS_iter(H, tol, 0, 0, N, G, L, Y, Z, X);
end



%% SUBROUTINES
% 
function H = BuildHSS_iter(H, tol, mh, nh, N, G, L, Y, Z,X)

m = size(H, 1);
n = size(H, 2);

if H.leafnode    
    % Let's do it for the HSS block row first.
    ridx = (mh+1:mh+m).';  
    %generate the diagonal block: 
    D = buildcauchy(mh+1:mh+m, nh+1:nh+n, N, G, L);
    
    % check rank: (for now this is for timings)
    rridx = (mh+1:mh+m).'; 
    ccidx = [(1:mh).'; (mh+m+1:N).'];
    w = exp(pi/N*1i); 
    I = [w^(2*(min(rridx)-1)), w^(2*(max(rridx)-1)) , w.^(2*ccidx(1) - 1), w.^(2*ccidx(end) - 1)];
    % compute cross-ratio:
    cr = get_cr(I); 
    [~,r] = size(G); 
    %bound on rank:
    pn = r*ceil(1/pi^2*log(4/tol)*log(16*cr));

    
    %%
    [~,p] = size(X); 
    Phi = Y(ridx, :) - D*X(ridx,:); 
    [H.U, J] = rand_interprow(Phi,p); 
    H.ridx = J + mh; 
    H.Phi = Phi; 
    H.J = J; 
    
    % And now do the columns
    cidx = (nh+1:nh+n).'; 
    Psi = Z(cidx, :) - D'*X(cidx, :); 
    [H.V, J] = rand_interpcol(Psi, p); 
    H.cidx = J +nh; %keep track of col indices. 
    H.Psi = Psi; 
    H.K = J; 
    H.D = D;
    H.hY = (H.V)'*X(ridx,:); 
    H.hZ = (H.U)'*X(cidx,:); 
    
else
    % Call the constructor recursively on the left and right children.
      H.A11 = BuildHSS_iter(H.A11, tol, mh, nh, N, G, L,Y, Z,X); 
      [mm, nn] = size(H.A11); 
      mh = mh + mm; 
      nh = nh + nn; 
      
    %check rank (for now this is for timings)
    rridx = (mh - mm +1: mh+mm).' ; 
    ccidx = [(1:mh-mm+1).'; (mh+mm+1:N).'];
    w = exp(pi/N*1i);
    I = [w^(2*(min(rridx)-1)), w^(2*(max(rridx)-1)) , w.^(2*ccidx(1) - 1), w.^(2*ccidx(end) - 1)];
    % compute cross-ratio:
    cr = get_cr(I); 
    [~,r] = size(G); 
    
    %bound on rank:
    pn = r*ceil(1/pi^2*log(4/tol)*log(16*cr));
      
    
    %%
      H.A22 = BuildHSS_iter(H.A22,tol, mh,nh,N, G, L, Y, Z,X);
      %now build glue matrices for level below:
      H.B12 = buildcauchy(H.A11.ridx, H.A22.cidx, N, G, L); 
      H.B21 = buildcauchy(H.A22.ridx, H.A11.cidx, N, G, L); 
      
      if H.topnode
        return;
      end
      
      %row translation matrices:
      idx1 = mh - mm +1; %square start %useful for debugging.
      idx2 = mh+mm; %square end 
      %update row span:   
      Phi1 = H.A11.Phi(H.A11.J,:) - H.B12*H.A22.hY; 
      Phi2 = H.A22.Phi(H.A22.J, :) - H.B21*H.A11.hY; 
      Phi = [Phi1; Phi2]; 
      p = length(H.A11.J); 
      [R, J] = rand_interprow(Phi, p);  
      a = length(H.A11.ridx); 
      H.Rl = R(1:a, :); 
      H.Rr = R(a+1:end, :); 
      rr = [H.A11.ridx; H.A22.ridx];
      H.ridx = rr(J); 
      H.J = J; 
      H.Phi = Phi; 
      
      %now col translation matrices: 
      idx1 = nh - nn +1; %square start
      idx2 = nh+nn; %square end
      %cidx = idx1:idx2;
      Psi1 = H.A11.Psi(H.A11.K,:) - H.B21'*H.A22.hZ; 
      Psi2 = H.A22.Psi(H.A22.K, :) - H.B12'*H.A11.hZ; 
      Psi = [Psi1; Psi2]; 
      p = length(H.A11.K); 
      [W, K] = rand_interpcol(Psi, p); 
      a = length(K); 
      H.Wl = W(1:a, :); 
      H.Wr = W(a+1:end, :); 
      rr = [H.A11.cidx; H.A22.cidx];
      H.cidx = rr(K); 
      H.K = K; 
      H.Psi = Psi; 
      
      %update randomized product matrices:
      H.hY = W'*[H.A11.hY; H.A22.hY]; 
      H.hZ = R'*[H.A11.hZ; H.A22.hZ]; 
              
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



function M = get_cr(I)
%given I = [a b c d] where [a b] [c d] are two disjoint intervals 
% on the real line, 
% M is the cross-ratio. 

a = I(1); 
b = I(2); 
c = I(3); 
d = I(4); 

%parameters
M = abs((c-a)*(d-b)/((c-b)*abs(d-a))); 
end


