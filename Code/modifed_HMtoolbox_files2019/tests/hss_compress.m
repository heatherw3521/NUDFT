% function that performs the compression of an HSS matrix
function B = hss_compress(A,tol)
B = hss_proper(A);
B = backward_stage(B,tol);
end

function A = backward_stage(A,tol,S,T)
if(A.leafnode==1)
    return
end
if (A.topnode == 1)
    [U,Su,V] = tsvd(A.Bu, tol);
    A.Bu = Su; Tl = Su.';
    if A.hssl.leafnode ==0
        A.hssl.Rl = A.hssl.Rl * U;
        A.hssl.Rr = A.hssl.Rr * U;
        A.hssr.Wl = A.hssr.Wl * V;
        A.hssr.Wr = A.hssr.Wr * V;
    else
        A.hssl.U = A.hssl.U * U;
        A.hssr.V = A.hssr.V * V;
    end
    [U,Sl,V] = tsvd(A.Bl, tol);
    A.Bl = Sl; Tu = Sl.';
    if A.hssl.leafnode ==0
        A.hssr.Rl = A.hssr.Rl * U;
        A.hssr.Rr = A.hssr.Rr * U;
        A.hssl.Wl = A.hssl.Wl * V;
        A.hssl.Wr = A.hssl.Wr * V;
        A.hssl = backward_stage(A.hssl, tol, Su, Tu);
        A.hssr = backward_stage(A.hssr, tol, Sl, Tl);
    else
        A.hssr.U = A.hssr.U * U;
        A.hssl.V = A.hssl.V * V;
    end
else
    Su = [A.Bu, A.Rl * S];
    Tl = [A.Bu.', A.Wr * T]; % possibile inghippo
    [Us,Su,Vs] = tsvd(Su, tol);
    [Ut,Tl,Vt] = tsvd(Tl, tol);
    k = size(A.Bu,2);
    A.Bu = Su * Vs(1:k, :)' * Ut;
    A.Rl = Us' * A.Rl;
    A.Wr = Ut' * A.Wr;
    if A.hssl.leafnode == 0
        A.hssl.Rl = A.hssl.Rl * Us;
        A.hssl.Rr = A.hssl.Rr * Us;
        A.hssr.Wl = A.hssr.Wl * Ut;
        A.hssr.Wr = A.hssr.Wr * Ut;
    else
        A.hssl.U = A.hssl.U * Us;
        A.hssr.V = A.hssr.V * Ut;
    end
    
    Sl = [A.Bl, A.Rr * S];
    Tu = [A.Bl.', A.Wl * T]; % possibile inghippo
    [Us,Sl,Vs] = tsvd(Sl, tol);
    [Ut,Tu,Vt] = tsvd(Tu, tol);
    k = size(A.Bl,2);
    A.Bl = Sl * Vs(1:k, :)' * Ut;
    A.Rr = Us' * A.Rr;
    A.Wl = Ut' * A.Wl;
    if A.hssr.leafnode == 0
        A.hssr.Rl = A.hssr.Rl * Us;
        A.hssr.Rr = A.hssr.Rr * Us;
        A.hssl.Wl = A.hssl.Wl * Ut;
        A.hssl.Wr = A.hssl.Wr * Ut;
        A.hssl = backward_stage(A.hssl, tol, Su, Tu);
        A.hssr = backward_stage(A.hssr, tol, Sl, Tl);
    else
        
        A.hssr.U = A.hssr.U * Us;
        A.hssl.V = A.hssl.V * Ut;
    end
end
end
