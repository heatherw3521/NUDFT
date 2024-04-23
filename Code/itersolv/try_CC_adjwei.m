% driver script for inv NUDFT 1D type 2: Clenshaw-Curtis grid w/ known quad wei.
% In particular, test solv=11 Kircheis-P'23 direct adj+wei (w solved by CGAN[B])
% Barnett 4/3/24

clear
signal = 'randn';  % image has decay coeffs
N = 2^10; % 2^10 small; or 2^16 = bigger test.  2^11 nudist=3 gives kappa~1e4
% for Kirchwei: following must be >=2.0:
dataratio = 2.2; %2.2; %1.53;   % pi/2 critical nud=2; Smaller kappa gro for nud>=3.
solvlist = [2 3 4 6 7 9 10 11];
opts.verb = 0;     % 1 = make plots
%opts.tol=1e-6;
%opts.addnoise = 1e-2;
nudistlist = 3;         % 2 (CC quad) here fails (kappa->inf) for dataratio<pi/2
out = test_iNUDFT(signal,N,dataratio,solvlist,nudistlist,opts);
% note solv=11 fails for nudist=3 when other iter methods succeed.  :(


% nud=1, N=2^14, dataratio=2.2 :
% $$$ ------------- NU dist=1 (jittered grid): -----------------
% $$$ 
% $$$ fast b and resid meas (tol=1e-12, no matrix A)...
% $$$ solver 2 (CG normal eqns)...
% $$$ 	flag=0, with 11 iters (rel resid nrm 3.34e-07)
% $$$ 	0.022 s   	rel l2 err 4.59e-07    	resid rel l2 nrm 3.79e-07
% $$$ solver 3 (Strang PCG nor eqns)...
% $$$ 	flag=0, with 11 iters (rel resid nrm 7.49e-07)
% $$$ 	0.0267 s   	rel l2 err 9.93e-07    	resid rel l2 nrm 8.41e-07
% $$$ solver 4 (adj matvec sinc2/quadr wei)...
% $$$ 	0.0114 s   	rel l2 err 0.144    	resid rel l2 nrm 0.135
% $$$ solver 6 (CG adj nor eqns)...
% $$$ 	flag=0, with 11 iters (rel resid nrm 3.96e-07)
% $$$ 	0.0794 s   	rel l2 err 4.41e-07    	resid rel l2 nrm 3.96e-07
% $$$ solver 7 (sinc2 PCG adj nor)...
% $$$ 	flag=0, with 10 iters (rel resid nrm 4.65e-07)
% $$$ 	0.0773 s   	rel l2 err 5.25e-07    	resid rel l2 nrm 4.65e-07
% $$$ solver 9 (INUDFT, wrapped)...
% $$$ 	1.87 s   	rel l2 err 5.59e-10    	resid rel l2 nrm 5.66e-10
% $$$ solver 10 (FDToep nor eqns)...
% $$$ 	2.83 s   	rel l2 err 1.1e-11    	resid rel l2 nrm 1.03e-11
% $$$ solver 11 (Kirch adj wei (by CGAN))...
% $$$ 	get w: CG flag=0, with 80 iters (rel resid nrm 8.88e-07)
% $$$ 	0.426 s   	rel l2 err 5.94e-07    	resid rel l2 nrm 5.97e-07

% we see Kircheis-Potts solv=11 not a threat, for M <= 2.2N