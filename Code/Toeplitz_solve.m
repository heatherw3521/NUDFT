function x = Toeplitz_solve(tc,tr, b, varargin)
% Solves the equation Tx =b, where T = toeplitz(tc, tr).
% 
% Toeplitz(tc, tr, b, 'tol', tol) sets a tolerance parameter.
% 

% This is a wrapper for applying @hss/private/hss_cauchy_toeplitz2

%%

%%
% check for tolerance parameter: 
tol = 1e-12; 
if ~isempty(varargin)
    if strcmpi(varargin{1}, 'tol')
        tol = varargin{2}; 
    else
        error('could not parse input')
    end
end
    


%% Part 1: transform to Cauchy system: 

%get transformed diagonal elements:


%dg = fast_toep2cauchydiags(tr, tc);
dg = fftToepDiag(tc,tr); 
 
n = length(tr); 
N = (1:n).'; 
w = exp(pi*1i/n); 

%get Toeplitz generators:
[GG, LL] = toep_gens_sing(tr, tc); 

%transform 2 Cauchy-like generators:
Dx = spdiags(w.^(N), 0, n,n); %diag row scaling
Dy = spdiags(w.^(-2+2*(N)), 0, n,n); %diag col scaling

G = sqrt(n)*Dx*ifft(Dy*GG);
L = sqrt(n)*Dx*ifft(Dy*LL);

%transform RHS: 
b = sqrt(n)*Dx*ifft(Dy*b); 

%% Part 2: HSS 
% build HSS approx to C: 
warning off
H = hss('cauchytoeplitz2', n, G, L, dg,'tol', tol);
warning on
%% Part 3: solve Cauchy-like system
% solve H x = b:
x = H\b; 

%% Part 4:
% transform x back: 
Dxi = spdiags(w.^(-N), 0, n,n); 
Dyi = spdiags(w.^((-1)*(-2+2*(N))),0,n,n);

x = Dyi*fft(Dxi*x)/sqrt(n); 
 
end

