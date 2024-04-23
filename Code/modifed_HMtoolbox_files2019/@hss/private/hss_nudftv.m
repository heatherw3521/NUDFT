function H = hss_nudftv(nodes, a, b,N, varargin)
% this code constructs an HSS factorization
% for the Cauchy-like matrix VF^*, where V is the nudft matrix
% and F is a scaled DFT. V is size M by N, where M = length(nodes)
%
% It preserves an interpolative structure so that
% the ULV-like solver in [1] can be used. 
% 
% We use ADI on the disp eqn D(VF^*) - (VF^*)L = a*b, where
% D = diag(nodes), L = diag(w^(2*(1:n))), w = exp(pi*1i/n).
 
block_size = 128; %smallest blocks in col dim.  
tol = hssoption('threshold'); %rel tol. 
if ~isempty(varargin)
    tolcheck = find(strcmp('tol', varargin));
    if ~isempty(tolcheck)
        tol = varargin{tolcheck+1};  
    end
end 

%% Build tree
M = length(nodes); 

nda = angle(nodes); jj = find(nda < 0,1); 
nda(jj:end) = nda(jj:end) + 2*pi; 
rou = 2*pi*(1:N)/N;
% Prepare the tree for the HSS structure -- leaving all the blocks empty:
%H2 = hss_build_hss_tree(M, N, block_size);
H = hss_build_hss_tree_nudft(M, N, nda, rou, block_size);

%build the HSS representation recursively:
H = BuildHSS_iter(H, tol, 0, 0,N, nodes, a, b');
end



%% SUBROUTINES
% 
function H = BuildHSS_iter(H, tol, mh, nh,N,nodes, a, b)

m = size(H, 1);
n = size(H, 2);

if H.leafnode

    %STEP 1: determine the rank: 
    w = exp(pi*1i/N); 

    % rank of HSS block row :
    ridx1 = (mh+1:mh+m).'; 
    cidx1 = [(1:nh) (nh+n+1:N)].'; 
    cidxl = circshift(cidx1, -nh);
    Ir = [w.^(2*(cidxl(1))), w.^(2*cidxl(end))];
    Ir = [flip(Ir.*[w.^(-1), w.^(1)]), Ir];
    [~, ~, ~, cr] = mobiusT(Ir);
    kr = ceil(1/pi^2*log(4/tol)*log(16*cr));

    % rank of HSS block column:
    M = length(nodes);
    %ridx2 = [(1:mh) (mh+m+1:M)].'; 
    %ridxl = circshift(ridx2, -mh); 
    cidx2 = (nh+1:nh+n).'; 
    Ic = [w.^(2*(cidx2(1))), w.^(2*cidx2(end))];
    Ic = [flip(Ic.*[w.^(-1), w.^(1)]), Ic];
    [~, ~, ~, cr] = mobiusT(Ic);
    kc = ceil(1/pi^2*log(4/tol)*log(16*cr));

    %[p, ~] = getshifts_adi(I, 'tol', tol);
    % I = interval on unit circle used for determining ADI shift params. 
    %I = [nodes(ridx(1)),nodes(
    % ridx(end)), w.^(2*(cidxl(1))), w.^(2*cidxl(end))];
     
    kk = min([kr,kc,m, n]); %smallest number is the rank
    

    % use fadi to build an ID for each block:

    %COL BLOCK
    bb = length(cidx2); 
    [H.V, J] = vfADI_col(bb, Ic,cidx2, N,kk, tol, b(cidx2)); %N is global size
    H.cidx = J +nh; %keep track of col indices. 
    kk = length(J); 
    %ROW BLOCK
    aa = length(ridx1);  
    [H.U, J] = vfADI_row(aa, Ir, ridx1, nodes,N, kk,tol, a(ridx1)); %N is global size 
    H.ridx = J + mh; %keep track of row indices. 
    
    %H.D = buildcauchy(mh+1:mh+m, nh+1:nh+n, N, nodes, a, b); 
    H.D = vbuildcauchydiags(nodes,mh+1:mh+m, nh+1:nh+n, N);
else
    % Call the constructor recursively on the left and right children.
      H.A11 = BuildHSS_iter(H.A11, tol, mh, nh, N, nodes, a, b);
    
      [mm1, nn] = size(H.A11); 
      [mm2, ~] = size(H.A22);
      mh = mh + mm1; 
      nh = nh + nn; 
      H.A22 = BuildHSS_iter(H.A22,tol, mh,nh,N, nodes, a, b);
      
      %now build glue matrices for level below:
      %H.B12 = buildcauchy(H.A11.ridx, H.A22.cidx, N, nodes, a, b); 
      %H.B21 = buildcauchy(H.A22.ridx, H.A11.cidx, N, nodes, a, b); 
      H.B12 = vbuildcauchydiags(nodes, H.A11.ridx, H.A22.cidx, N); 
      H.B21 = vbuildcauchydiags(nodes, H.A22.ridx, H.A11.cidx, N); 

      %try something: 
      %B12 = vbuildcauchydiags(nodes, H.A11.ridx, H.A22.cidx, N);

      if H.topnode
        return;
      end
      % translation matrices:

      %estimate row rank
      w = exp(1i*pi/N); 
      ridx1 = [H.A11.ridx; H.A22.ridx]; %subset of rows;
      idx1 = nh - nn +1; %col left
      idx2 = nh+nn; %col right
      cidx1 = [(1:idx1-1) (idx2+1:N)].';
      cidxl = circshift(cidx1, -(nh-nn));
      Ir = [w.^(2*(cidxl(1))), w.^(2*cidxl(end))];
      Ir = [flip(Ir.*[w.^(-1), w.^(1)]), Ir];
      [~, ~, ~, cr] = mobiusT(Ir);
      kr = ceil(1/pi^2*log(4/tol)*log(16*cr));
      aa = length(ridx1);

      %estimate col rank
      %M = length(nodes);
      %idx1 = mh - mm1 +1; %row top
      %idx2 = mh+mm2; %row bottom
      %ridx2 = [(1:idx1-1) (idx2+1:M)].';
      cidx2 = [H.A11.cidx; H.A22.cidx]; 
      bb = length(cidx2); 
      %ridxl = circshift(ridx, -(mh-mm)); 
      Ic = [w.^(2*min(cidx2)), w.^(2*max(cidx2))];
      Ic = [flip(Ic.*[w.^(-1), w.^(1)]), Ic];
      [~, ~, ~, cr] = mobiusT(Ic);
      kc = ceil(1/pi^2*log(4/tol)*log(16*cr));
      %bb = length(cidx2);
      kk = min([kr,kc,m, n]); %smallest number is the rank

      

      %build col translation matrices: 
      [R, K] = vfADI_col(bb,Ic,cidx2, N, kk,tol, b(cidx2)); 
      bb = length(H.A11.cidx);
      H.Wl = R(1:bb, :); 
      H.Wr = R(bb+1:end, :); 
      H.cidx = cidx2(K);
      kk = length(K); 
       
      
      %build row translation matrices: 
      [R, J] = vfADI_row(aa, Ir, ridx1, nodes, N,kk,tol, a(ridx1));
      aa = length(H.A11.ridx); 
      H.Rl = R(1:aa, :); 
      H.Rr = R(aa+1:end, :); 
      H.ridx = ridx1(J); 
      

    %     %DELETE WHEN DONE begin error check: 
%     testrows = buildcauchy(ridx, cidx, N, nodes, a, b);
%     buildcheck = R*testrows(J,:);
%     if norm(buildcheck-testrows) > 1e-5
%         error('error line 136 (transl. block row)')
%     end
%     %end error check


     
       


    %     %DELETE WHEN DONE begin error check: 
%     testcols = buildcauchy(ridx, cidx, N, nodes, a, b);
%     buildcheck = testcols(:,K)*R';
%     if norm(buildcheck-testcols) > 1e-5
%         error('error line 162 (transl. block col)')
%     end
%     %end error check
end
    
end


%%

function [T, Tinv, gam, M] = mobiusT(I)
%given I = [a b c d] where [a b] [c d] are two disjoint intervals 
% on the real line, T(I) maps to the four points [-gamma, -1, 1, gamma]. 
% M is the cross-ratio. 

a = I(1); 
b = I(2); 
c = I(3); 
d = I(4); 

%parameters
M = abs((c-a)*(d-b)/((c-b)*abs(d-a))); 
gam = -1+2*M+2*sqrt(M^2-M); 
A = -gam*a*(1-gam)+gam*(c-gam*d)+c*gam-gam*d; 
B = -gam*a*(c*gam-d)-a*(c*gam-gam*d)-gam*(c*d-gam*d*c); 
C = a*(1-gam)+gam*(c-d)+c*gam-d; 
D = -gam*a*(c-d)-a*(c-gam*d)+c*d-gam*d*c; 

T = @(z) (A*z+B)./(C*z+D);
Tinv = @(z) (D*z-B)./(-C*z+A);
end







