% driver script for inv NUDFT 1D type 2 with plots, showing signal & recons.
% Barnett 2/2/23

clear
signal = 'image';  % image has decay coeffs
N = 2^14; % 2^10 small; or 2^16 = bigger test.  2^11 nudist=3 gives kappa~1e4
dataratio = 1.8;    % 1.8 for run 11/18/22. Smaller gives kappa gro for nud>=3
%solvlist = [2 5 6 9 10];
solvlist = [2 6 9 10];
opts.verb = 1;     % 1 = make plots
nudistlist = 3;         % 2 (CC quad) here fails for dataratio<pi/2
out = test_iNUDFT(signal,N,dataratio,solvlist,nudistlist,opts);
