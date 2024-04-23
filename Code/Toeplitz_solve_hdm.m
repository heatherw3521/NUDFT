
function x = Toeplitz_solve_hdm(tc,tr, b, varargin)
% Solves the equation Tx =b, where T = toeplitz(tc, tr).
% using the nonsingular version of the disp struct
% Toeplitz(tc, tr, b, 'tol', tol) sets a tolerance parameter.
% 

% This is a wrapper for @hss/private/hss_cauchy_toeplitz2

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


%C = toeplitz2cauchy(tr, tc,'nonsing');

n = length(tr); 
N = (1:n).'; 
w = exp(pi*1i/n); 

%get Toeplitz generators:
        [GG, LL] = toep_gens(tr, tc); 
        %CC = build_cauchy(n); 

        %transform 2 Cauchy-like generators:
        D0 = spdiags(w.^(-(N-1)), 0, n,n);
        G = sqrt(n)*ifft(GG);
        L = sqrt(n)*ifft(D0'*LL);

%transform RHS: 
b = sqrt(n)*ifft(b); 

%% Part 2: HSS 
% build HSS approx to C: 
warning off
H = hss('cauchytoeplitzunstr', n, G, L,'tol', tol);
warning on
%% Part 3: solve Cauchy-like system
% solve H x = b:
x = H\b; 

%% Part 4:
% transform x back: 

x = D0*fft(x)/sqrt(n);  
 
end

