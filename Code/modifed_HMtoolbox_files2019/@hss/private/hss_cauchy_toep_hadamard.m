function H = hss_cauchytoep_hadamard(varargin)
%
% builds an HSS factorization for a cauchy-like matrix
% using the hadamard product. Calls cauchytoeplitzc. 
% this is a non-interpolative structured version of the
% factorization
%
% hss_cauchytoeplitz(t1,t2, varargin)
% c1 and c2 are N by 2 matrices and
% C \circ (c1*c2.') is the Cauchy-like matrix
% that we want to approximate. 
%
% hss_cauchytoeplitz(N, varargin), 
% where N is the size of C, returns just the
% Cauchy matrix (N by N). (i.e., c1 = c2 = ones(N,1)$. 
%
% hss_cauchytoeplitz(..., 'tol', tol) sets tolerance. 
 

%% Part 1: parse input

if  all((1 == size(varargin{1})))
    %first argument is size of cauchy. 
    H = hss_cauchytoeplitzc(varargin{:}); 
    return; 
end
%otherwise, we build cauchy matrix first, then 
% build the rest from it:  

%set tol:
tol = hssoption('threshold'); %rel tol. 
if ~isempty(varargin)
    tolcheck = find(strcmp('tol', varargin));
    if ~isempty(tolcheck)
        tol = varargin{tolcheck+1};  
    end
end

%get vectors and set size: should be rhs = FG'
G = varargin{1}; [a, b] = size(G); 
if (a < b) 
    G = G.'; 
end
F = varargin{2}; [a, b] = size(F); 
if (a < b) 
    F = F.';
end
N = max(size(G)); 

%% part 2: cauchy matrix
%get Cauchy matrix: 
H1 = hss_cauchytoeplitzc(N, 'tol', tol); 

%% part 3: apply hadamard product factors
% build the Cauchy-like factorization on top of the Cauchy HSS structure:
mh = 0; nh=0; 
[H, ~, ~] = hss_hadamardcauchy(H1, G, F,mh, nh, N); 


end

%% subroutines
function [Y, mh, nh] = hss_hadamardcauchy(H, G, F,mh, nh, N)
Y = H; 
if H.leafnode 
    [m, n] = size(H); 
    %pick out vectors
    dridx = mh+1:mh+m; dG1 = G(dridx, 1); dG2 = G(dridx, 2);
    dcidx = nh+1:nh+n; dF1 = F(dcidx,1); dF2 = F(dcidx,2);
 
    %new factors
    Y.U = [diag(dG1)*H.U, diag(dG2)*H.U]; % U = [D_1U D_2U] 
    Y.V = [diag(dF1)*H.V, diag(dF2)*H.V];  % V' = [V'D_1' V'D_2'] 
    
    %new diagonal
    Y.D = diag(dG1)*H.D*diag(dF1') + diag(dG2)*H.D*diag(dF2'); 
    mh = mh+m; 
    nh = nh+n; 
else
    
    %recurse 
    [Y.A11, mh, nh] = hss_hadamardcauchy(H.A11, G, F, mh, nh, N);
    [Y.A22, mh, nh] = hss_hadamardcauchy(H.A22, G, F, mh, nh, N); 
    

    %set glue blocks for one level below: 
    Y.B12 = blkdiag(H.B12, H.B12); 
    Y.B21 = blkdiag(H.B21, H.B21); 
    
    if H.topnode
        return;
    end
    
    %set transition matrices 
    Y.Rl = blkdiag(H.Rl, H.Rl);   
    Y.Rr = blkdiag(H.Rr, H.Rr); 
 
    
    Y.Wl = blkdiag(H.Wl, H.Wl);   
    Y.Wr = blkdiag(H.Wr, H.Wr);     
end




end





