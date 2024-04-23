function [x,H] = Toeplitz_solve_ns(tc,tr, b, varargin)
% Solves the equation Tx =b, where T = toeplitz(tc, tr).
% using the nonsingular version of the disp struct
% Toeplitz_solve_ns(tc, tr, b, 'tol', tol) sets a tolerance parameter.

% This is a wrapper for @hss/private/hss_cauchy_toeplitz

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

n = length(tr); 
block_size = 128;

%% Part 1: transform to Cauchy system: 

if ~(log2(n)==floor(log2(n)))
    error('For now, size of Toeplitz matrix must be a power of 2')
end
N = (1:n).'; 
w = exp(pi*1i/n); 

%get Toeplitz generators:
[GG, LL] = toep_gens(tr, tc); 
        
%transform to Cauchy-like generators:
D0 = spdiags(w.^(-(N-1)), 0, n,n);
G = sqrt(n)*ifft(GG);
L = sqrt(n)*ifft(D0'*LL);

%transform RHS: 
b = sqrt(n)*ifft(b); 

%% Part 2: HSS 
% build HSS approx to C:

if n < block_size %if block size is very small, we want to skip fadi (not low rank!)
    aa = w.^(2*(0:n-1)); 
    C = bsxfun(@minus, aa.', w*aa); 
    C = 1./C; 
    H = hss((G*L').*C);
else
warning off
H = hss('cauchytoeplitz', n, G, L,'tol', tol);
warning on
end
%% Part 3: solve Cauchy-like system

x = H\b; 

%% Part 4:
% transform x back: 

x = D0*fft(x)/sqrt(n);  
 
end

