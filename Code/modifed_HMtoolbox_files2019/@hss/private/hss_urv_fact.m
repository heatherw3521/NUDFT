function F = hss_urv_fact(A)
% HSS_URV_FACT computes the  URV factorization
%	       of A
%

F = struct();
F.topnode = A.topnode;
F.leafnode = A.leafnode;
F.B12 = A.B12;
F.B21 = A.B21;
F.Rl = A.Rl;
F.Rr = A.Rr;
if A.leafnode
    F.left = [];
    F.right = [];
    F.m = size(A.D, 1);
    F.n = size(A.D, 2);
    D = A.D;
    U = A.U;
    V = A.V;
else
    F.left = hss_urv_fact(A.A11);
    F.right = hss_urv_fact(A.A22);
    F.m = F.left.m + F.right.m;
    F.n = F.left.n + F.right.n;
    D = [F.left.D22 F.left.U2 * A.B12 * F.right.V';
         F.right.U2 * A.B21 * F.left.V' F.right.D22];
    if isempty(A.Rl) && isempty(A.Rr) && isempty(A.Wl) && isempty(A.Wr)
        Rl = zeros(size(F.left.U2, 2), 0);
        Rr = zeros(size(F.right.U2, 2), 0);
        Wl = zeros(size(F.left.V, 2), 0);
        Wr = zeros(size(F.right.V, 2), 0);
    else
        Rl = A.Rl; Rr = A.Rr; Wl = A.Wl; Wr = A.Wr;
    end
    U = [F.left.U2 * Rl;F.right.U2 * Rr];
    V = [F.left.V * Wl;F.right.V * Wr];
end

if A.topnode
    [F.Q, F.D11] = qr(D);
    return
end

if size(U, 1) >= size([U D], 2)
    %[F.Om,tmp] = qr([U D]);
    [F.Om,tmp] = qr([U D],"econ");
    k = min(size([U D])); l = size(U,2);
    U = tmp(1:k,1:l);
    D = tmp(1:k,l+1:end);
else
    F.Om = [];
end

r = size(V, 2);
[F.P,tmp] = qr(V); F.P = F.P(:,end:-1:1);
F.V = tmp(r:-1:1,:);

D = D * F.P;
[F.Q,~] = qr(D(:,1:end-r));
D = F.Q' * D; U = F.Q' * U;

k = size(D,2) - r;
F.D11 = D(1:k,1:k); F.D12 = D(1:k,k+1:end);
F.D22 = D(k+1:end,k+1:end);
F.U1 = U(1:k,:); F.U2 = U(k+1:end,:);

end
