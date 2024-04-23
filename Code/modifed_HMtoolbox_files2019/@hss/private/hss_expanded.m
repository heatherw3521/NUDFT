function [expanded, perm] = hss_expanded(A,alpha,beta)
%HSS_EXPANDED     form the perfect shuffled expanded matrix from A
%
%	       hssA = HSS_EXPANDED(hssA, alpha, beta) forms the perfect
%	       shuffle of the expanded matrix [alpha*eye hssA;hssA' beta*eye]

expanded = A;

if A.leafnode
    expanded.D = [alpha*eye(size(A.D,1)) A.D;
                      A.D' beta*eye(size(A.D,2))];
    expanded.U = blkdiag(A.U, A.V);
    expanded.V = expanded.U;
    perm = 1:size(expanded.D,1);
    return
end

expanded.ml = A.ml+A.nl;
expanded.nl = A.ml+A.nl;
expanded.mr = A.mr+A.nr;
expanded.nr = A.mr+A.nr;

expanded.B12 = [zeros(size(A.B12,1),size(A.B21,1)) A.B12;
                    A.B21' zeros(size(A.B21,2),size(A.B12,2))];
expanded.B21 = expanded.B12';
expanded.Rl = blkdiag(A.Rl, A.Wl);
expanded.Wl = expanded.Rl;
expanded.Rr = blkdiag(A.Rr, A.Wr);
expanded.Wr = expanded.Rr;

[expanded.A11, perm_l] = hss_expanded(A.A11,alpha,beta);
[expanded.A22, perm_r] = hss_expanded(A.A22,alpha,beta);

m1 = size(A.A11,1);
m2 = m1 + size(A.A22,1);
m3 = m2 + size(A.A11,2);
m4 = m3 + size(A.A22,2);
perm = [1:m1 (m2+1):m3 (m1+1):m2 (m3+1):m4];
n = length([1:m1 (m2+1):m3]);

perm1 = perm(1:n);
perm1 = perm1(perm_l);
perm2 = perm((n+1):end);
perm2 = perm2(perm_r);

perm = [perm1 perm2];

end

