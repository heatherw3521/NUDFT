function C = hss_hadamard_mul(A,B);
tol = hssoption('threshold');
C = hss_hadamard_mul_ric(A,B);
C = hss_compress(C,tol);
end

function C = hss_hadamard_mul_ric(A,B)
C = hss();
if  A.topnode ~= B.topnode || A.leafnode ~= B.leafnode
    error('HSS_HADAMARD_MUL:: the two hss matrices have not compatible partitioning')
else
    C.topnode = A.topnode;
    C.leafnode = A.leafnode;
    if C.leafnode == 0
        if A.ml ~= B.ml || A.nl ~= B.nl ||A.mr ~= B.mr ||A.nr ~= B.nr
            error('HSS_HADAMARD_MUL:: the two hss matrices have not compatible partitioning')
        end
        C.ml = A.ml;
        C.nl = A.nl;
        C.mr = A.mr;
        C.nr = A.nr;
    end
end

if C.leafnode == 1
    C.D = A.D .* B.D;
    C.U = [];
    C.V = [];
    for j = 1:size(A.U,2)  % Is there a better way?
        C.U = [C.U, diag(A.U(:,j)) * B.U];
    end
    for j = 1:size(A.V,2)  % Is there a better way?
        C.V = [C.V, diag(A.V(:,j)) * B.V];
    end
else
    C.B12 = kron(A.B12, B.B12);
    C.B21 = kron(A.B21, B.B21);
    
    if C.topnode == 0
        C.Rl = kron(A.Rl, B.Rl);
        C.Rr = kron(A.Rr, B.Rr);
        C.Wl = kron(A.Wl, B.Wl);
        C.Wr = kron(A.Wr, B.Wr);
    end
    C.A11 = hss_hadamard_mul_ric(A.A11, B.A11);
    C.A22 = hss_hadamard_mul_ric(A.A22, B.A22);
end
end
