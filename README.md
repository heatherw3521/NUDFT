# NUDFT
MATLAB code for developing a fast least-squares solver for type-II NUDFT

This repository contains code for a superfast solver for the inverse type-II nonuniform discrete Fourier transform problem. Code includes all experiments described in the 
paper "A superfast direct inversion method for the nonuniform discrete Fourier transform."

The NUDFT solver itself has no package dependencies, but the implementation of iterative methods, other direct methods, and related tests requires at various points the FINUFFT package (and the FFTW package), the NFFT package. Less critical packages that are called at various points include the Chebfun package (only to produce Clenshaw-Curtis points in certain files), the Altmany/export-fig package (to produce PDFs of figures found in the paper), along the linspecer package (for plot colors).  
