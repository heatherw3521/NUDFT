# NUDFT
A fast least-squares solver for type-II NUDFT linear systems in MATLAB

*Authors: Heather Wilber, Ethan N. Epperly, and Alex H. Barnett*

This repository contains code for a superfast solver for the inverse type-II nonuniform discrete Fourier transform problem, in one dimension. There are also implementations of various iterative and direct solution methods for the same problem. This code includes all experiments described in the preprint:

"A superfast direct inversion method for the nonuniform discrete Fourier transform," H. Wilber, E. N. Epperley, and A. H. Barnett, 26 pages,
submitted, SIAM J. Sci. Comput., 2024. https://arxiv.org/abs/2404.13223

It you find this code useful, please cite that paper.

*Warning*: this is currently research code that will gradually get cleaned up over time (at which point, interfaces may break). Use at your own risk.


### Example use

After cloning this repo, `cd NUDFT`, open MATLAB, then try this basic
small example which solves the linear system $Vx = b$
for a single right-hand side vector, using the
conventions of the above preprint:
```matlab
addpath Code Code/modifed_HMtoolbox_files2019 Code/extras
n = 2000;                % unknowns, number of regular Fourier frequencies
m = 2*n;                 % equations, number of nonuniform nodes
p = rand(m,1);           % iid random nodes in [0,1), needn't be sorted
gamma = exp(2i*pi*p);    % their Vandermonde nodes on unit circle
w = 0:n-1;               % row vec of the frequency indices
xtrue = randn(n,1) + 1i*randn(n,1);    % choose Fourier coeffs at said indices
V = gamma(:).^w;         % dense fill Vandermonde matrix (only for small tests)
b = V*xtrue;             % direct generate RHS (only for small tests)
tic;                     % time our solver
x = INUDFT(gamma,n,b);   % solve one RHS with default tolerance
fprintf('solve in %.3g s: rel l2 resid norm %.3g, rel l2 soln error norm %.3g\n', toc, norm(V*x-b)/norm(b), norm(x-xtrue)/norm(x))
```
This completes in less than a second with the output:
```
done in 0.23 s: rel l2 resid norm 3.81e-11, rel l2 soln error norm 2.96e-10
```
For comparison, MATLAB's dense direct solver takes 7 seconds for this problem.
The condition number is $\kappa_2(V) \approx 1.03 \times 10^3$ (for random seed `rng(0)`).
For complete code see `small_example.m` in the top directory.

Next, look in `Code/itersolv` and its README for documented implementations of all iterative solution methods, `Code/wrapper_INUDFT`,
and the comparison code `test_iNUDFT.m`. They all use a different convention
with nodes in $[0,2\pi)$ and symmetrized frequencies, as in FINUFFT.
Also see [FINUFFT iterative inversion tutorial](https://finufft.readthedocs.io/en/latest/tutorial/inv1d2.html) which uses the latter notation for this same task.


### Installation

The inverse NUDFT solver itself has no package dependencies and should
run on any recent MATLAB version (eg R2021 or later).

To generate data for larger tests in linear time, or to run the iterative solvers, one needs to install:

   * [FINUFFT](https://finufft.readthedocs.io) (in particular, its MATLAB interface). Please compile its MEX file, insure FINUFFT MATLAB tests pass, then `addpath your-path-to-finufft/matlab`.

Some of our comparions tests, figure-generating codes, and benchmarking need the additional packages:

     * [NFFT3.5](https://github.com/NFFT/nfft)  (in particular, its MATLAB interfaces in `infft1d`). Please compile its MEX files, run its self-tests, then insure its file `infft.m` is in your path. This is needed for comparison tests to the Kircheis--Potts sparse direct algorithm of 2022.
     * [chebfun](https://www.chebfun.org/) package (only to produce Clenshaw--Curtis points in certain files)
     * [export-fig](https://github.com/altmany/export_fig) package (to produce PDFs of figures found in the paper)
     * [linspecer](https://github.com/davidkun/linspecer) package (for plot colors)
     * [memorygraph](https://github.com/ahbarnett/memorygraph) for measuring RAM and CPU usage vs time

Please contact us if you have problems running the basic example above, but not for problems running the comparison tests, figure-generating codes, or integration with other packages listed above.


