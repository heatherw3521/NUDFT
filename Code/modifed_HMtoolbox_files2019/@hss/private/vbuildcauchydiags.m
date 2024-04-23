function C = vbuildcauchydiags(nodes, ridx, cidx, n)
%build diagonal block sumatrices C(ridx,cidx) where C = VF'. 
%this constructs the submatrix entries directly using a formula based 
% on geometric sums rather than using the Sylvester equation, which
% can be numerically unstable

nodes = nodes(:);
ridx = ridx(:); cidx = cidx(:)';
%V = nodes.^(0:n-1);
%out = (ft(V'))';
%m = length(ridx);
%nn = length(cidx);
%lam = angle(nodes);
lam = angle(nodes(ridx));
%C = zeros(m,nn); 
ww = exp(-1i*pi/n);

%geometric sum is CC(j,k) = (nodes(j).^(0:n-1))*( (ww^(k)).^(2*(1:n)-1).' )/sqrt(n);
% for j = 1:m
%     jj = ridx(j);
%     for k = 1:nn
%         kk = cidx(k); 
%         if abs(1-nodes(jj)*ww^(2*kk)) < 1e-9 %singularity
%             C(j,k) = ww^(kk)/sqrt(n)*(1 - ( nodes(jj)*ww^(2*kk) )^n)./...
%             (1-nodes(jj)*ww^(2*kk));
%         else
%             C(j,k) = ww^(kk)/sqrt(n)*(expm1(1i*lam(jj)*n))./...
%                 expm1(1i*lam(jj)-2*pi*1i*kk/n);
%         end
%     end
% end
 

C = ( ww.^(cidx)./(expm1(1i*lam-2*pi*1i*cidx/n)) ).*( expm1(1i*lam*n))/sqrt(n);

 
%check for singularity: 
filter = abs(1-nodes(ridx).*(ww.^(2*cidx))) < 1e-8;
if any(filter, 'all')
    Calt = (ww.^(cidx).*(1-(nodes(ridx).*(ww.^(2*cidx))).^n))...
        ./(1-nodes(ridx).*(ww.^(2*cidx)))/sqrt(n); 
C(filter) = Calt(filter); 
end
%C = Calt;
end
%%

%function out = ft(x)
    % applies Fx, where
    % F_{j,k} = exp(2*pi*1i*j*(2k-1))/sqrt(length(x)). 
   % [n, ~] = size(x);
%    w = exp(pi*1i/n);
 %   J = spdiags(sqrt(n)*w.^((1:n)'), 0, n,n); W = spdiags(w.^(2*(1:n)'-2), 0, n, n);
  %  y = full(W*x); 
  %  out = full(J*ifft(y)); 
%end