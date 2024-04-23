function F = urv(A)
%ULV   Computes the URV factorization of A.
%
%      F = URV(A) returns a structure containing a
%      parametrization of the URV of A.
%      This can be used to compute X = A\B with the command
%      X = URV_SOLVE(F, B);
F = hss_urv_fact(A);
end
