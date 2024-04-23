function out = my_fft(x)
    % applies Fx, where
    % F_{j,k} = exp(2*pi*1i*j*(2k-1))/sqrt(length(x)). 
    % this is a nonconventional DFT
    [n, ~] = size(x);
    w = exp(pi*1i/n);
    J = spdiags(sqrt(n)*w.^((1:n)'), 0, n,n); W = spdiags(w.^(2*(1:n)'-2), 0, n, n);
    y = full(W*x); 
    out = full(J*ifft(y)); 
end