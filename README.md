# NUDFT
MATLAB code for developing a fast least-squares solver for type-II NUDFT

*Authors: Heather Wilber, Ethan N. Epperly, and Alex H. Barnett*

This repository contains code for a superfast solver for the inverse type-II nonuniform discrete Fourier transform problem, in one dimension. There are also implementations of various iterative solution methods for the same problem. This code includes all experiments described in the preprint:

"A superfast direct inversion method for the nonuniform discrete Fourier transform," H. Wilber, E. N. Epperley, and A. H. Barnett, 26 pages,
submitted, SIAM J. Sci. Comput., 2024. https://arxiv.org/abs/2404.13223

It you find this code useful, please cite that paper.

# Installation

The NUDFT solver itself has no package dependencies, but the implementation of iterative methods, other direct methods, and related tests requires at various points the FINUFFT package (and the FFTW package), the NFFT package. Less critical packages that are called at various points include the Chebfun package (only to produce Clenshaw-Curtis points in certain files), the Altmany/export-fig package (to produce PDFs of figures found in the paper), along the linspecer package (for plot colors).  
