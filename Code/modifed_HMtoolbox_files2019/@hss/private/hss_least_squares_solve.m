function x = hss_least_squares_solve(A, b)
% HSS_LEAST_SQUARES_SOLVE computes the least squares solution to the linear
%          system A x = b by means of the expanded equations
%

F = urv(A);
x = urv_solve(F, b);

end