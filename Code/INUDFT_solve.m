function X = INUDFT_solve(L,p,B)
%solves VX = B where L stores permuted urv factorization for the transformed 
%matrix V. p is the permutation info.  B is multiple RHS. 
%
% Example: First call [L, ~] = INUDFT(nodes, n, b1); to produce L, where 
% nodes define A and A is m by n. 
% Now solve system AX = B with X = INUDFT_solve(L,p, B); 

B = B(p,:);
X = urv_solve(L, B);
X = ift(X);

end

function out = ift(x)
    %inverse of ft(x)
    [n,~] = size(x);
    w = exp(pi*1i/n);
    J = spdiags((1/(sqrt(n)))*w.^(-(1:n)'), 0, n,n); W = spdiags(w.^(-2*(1:n)'+2), 0, n, n);
    y = full(J*x);
    out = full(W*fft(y)); 
end