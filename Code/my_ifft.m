function out = my_ifft(x)
    %inverse of my_fft(x)
    [n,~] = size(x);
    w = exp(pi*1i/n);
    J = spdiags((1/(sqrt(n)))*w.^(-(1:n)'), 0, n,n); W = spdiags(w.^(-2*(1:n)'+2), 0, n, n);
    y = full(J*x);
    out = full(W*fft(y)); 
end