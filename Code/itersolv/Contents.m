% ITERSOLV
%
% Files
%   circ_apply            - Fast mat-vec by circulant mat C given DFT of 1st row of C
%   circulant             - fill dense circulant matrix given its 1st col
%   clencurt              - nodes (Chebyshev) and weights, Clenshaw-Curtis quadrature, via FFT
%   densemat_nudft        - fill dense 1D type-2 NUDFT matrix (FINUFFT convention)
%   densemat_nudft_vander - dense type-2 NUDFT matrix in Vandermonde convention
%   do_nufft_vander       - Wrapper to FINUFFT 1d type-2 from Vandermonde freq convention
%   mergestructs          - combine two structs with distinct fields
%   optfrobwei            - optimal (in Frobenius) norm diag gridding "sinc^2" weights
%   solv_adjwei           - diag-weighted adjoint "gridding" approx soln of t2 NUDFT lin sys
%   solv_CGAN             - conjugate gradients on adjoint normal eqns, 1D NUDFT type 2 lin sys
%   solv_CGN              - conjugate gradients on normal eqns to solve 1D NUDFT type 2 lin sys
%   solv_dense            - plain direct solution of 1D NUDFT (type 2) linear system
%   solv_FPadjwei         - fixed point iter on diag-weighted adj, solves t2 NUDFT lin sys
%   strangprecond         - Strang circulant preconditioner for square nonsymm Toeplitz
%   test_iNUDFT           - test script for inverting 1D type 2 NUDFT: various methods on various problems
%   Toep_apply            - fast matvec with square Toeplitz matrix (good pad)
%   Toep_apply_badpad     - fast matvec with square Toeplitz matrix (obsolete)
