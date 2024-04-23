function d = fftToepDiag(c,r) 

c = c(:); 
r = r(:); r= r.'; 
n = length(c); 
%c = rand(n,1); r = rand(1,n); r(1) = c(1);
%T = toeplitz(c,r);

% Compute orthogonal projection of Toeplitz matrix onto the subspace of cyclic matrices
rnew = ( r(2:end).*[n-1:-1:1] + c(end:-1:2)'.*[1:n-1] ) / n;

%C = toeplitz([c(1),rnew(end:-1:1)]',[r(1),rnew]);
cc = [c(1),flip(rnew(1:end))]';
% test orthogonality
%trace((T-C)'*C)

% Compute diagonal of transformed null space projection
% the DFT we apply is unconventional: 

%F = fft(eye(n)) / sqrt(n);
%td = diag( F*C*F' );
%fd = fft(C(:,1));
%norm( td - fd )
w = exp(pi*1i/n); 
N = (1:n).'; 
Dx = w.^(N); %diag row scaling
Dy = w.^(-2+2*(N)); %diag col scaling


%d = sqrt(n)*Dx.*(ifft(Dy.*cc));
d = flip(fft(cc)); 
end
