 function C = toeplitz2cauchy(tr, tc, type, varargin)
%transforms a Toeplitz matrix to a Cauchy-like 
% matrix C. 
% If type = 'sing', we use the transform
% associated with the singular displacement
% equation QT-TQ, where Q = circshift(eye(n)) 
% is the circumshift matrix.
%
% If type = 'nonsing', we use the transform
% associated with a nonsingular displacement equation, 
% QT - TU, where U = Q except that U(1,end) = -1.
%
% if the keyword 'full' is used, then this code
% actually constructs the Toeplitz matrix and takes
% the transform. If no keyword is used, the code constructs
% C using generators. 


n = length(tc); 
w = exp(pi*1i/n);
N = (1:n).';

gens = 1; %full version or generator version
if~isempty(varargin)
    v1 = varargin{1}; 
    if strcmpi(v1, 'full')
        gens = 0;
    end
end

if gens==0 
T = toeplitz(tc, tr); %make Toeplitz
    if strcmpi(type, 'sing') || strcmpi(type, 'singular')
        Dx = spdiags((w.^N), 0, n,n); %diag row scaling
        Dy = spdiags(w.^(-2+2*N), 0, n,n); %diag col scaling
        Dxi = spdiags((w.^(-N)), 0, n,n);
        Dyi = spdiags(w.^(-(-2+2*N)), 0, n,n);
        
        C = Dx*ifft(Dy*T); 
        C = (Dxi*fft(Dyi*C.')).';
        return;
    else
        D0 = spdiags(w.^(-(N-1)), 0, n,n); 
        C = ifft(T*D0); 
        C = (fft(C.')).';
        return;
    end
else %use generators to make C: 
    if strcmpi(type, 'sing') || strcmpi(type, 'singular')
        [Gs, Ls] = toep_gens_sing(tr, tc); 
        
        
        %now write using FFTs: 
        Dx = spdiags(w.^(N), 0, n,n); %diag row scaling
        Dy = spdiags(w.^(-2+2*(N)), 0, n,n); %diag col scaling

        G = sqrt(n)*Dx*ifft(Dy*Gs);
        L = sqrt(n)*Dx*ifft(Dy*Ls);
        
        CC = build_cauchys(n); 
        CC(isinf(CC))=0; 
        
        C = 0; 
        for j = 1:2
            B = spdiags(G(:, j), 0,n,n)*CC*spdiags(conj(L(:,j)), 0, n,n);
            C = C +B; 
        end
        
        cdiag = fast_toep2cauchydiags(tr, tc); 
        C = C + diag(cdiag); 
        return
    else
        
        [GG, LL] = toep_gens(tr, tc); 
        CC = build_cauchy(n); 
        
        D0 = spdiags(w.^(-(N-1)), 0, n,n);
        G = sqrt(n)*ifft(GG);
        L = sqrt(n)*ifft(D0'*LL);
        
        %C = CC.*(G*L'); 
        C = 0; 
        for j = 1:2
            B = spdiags(G(:, j), 0,n,n)*CC*spdiags(conj(L(:,j)), 0, n,n);
            C = C+B; 
        end
    end
end

        
        

