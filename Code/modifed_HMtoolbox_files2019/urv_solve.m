function x = urv_solve(F, b)
%URV_SOLVE     solve the system A X = B where F contains the URV factorization of A
%
%	       X = URV_SOLVE(F, B) computes A\B with F = URV(A);
if (~isstruct(F))
    error('F is not of the correct type');
end
x = hss_urv_fact_solve(F, b);
end
