% crappy script to test iNUDFT on big probs (1D).
% Needs: memorygraph.  Barnett 9/14/23
addpath ..

clear
signal = 'image';
Ns=[2^19 2^20] % both used >10GB.    2^21 2^22]; 2^21 not on laptop.
dataratio = 1.8;    % as for run 11/18/22
solvlist = [2,6,9,10];   % only include 1 if N < few k.
nudist = 3;
opts.tol = 1e-6;
opts.cgtol = opts.tol;
fnam = sprintf('results/%s_NU%d_big_dataratio%g_tol%g.txt',signal,nudist,dataratio,opts.tol)
system(['rm -f ',fnam]);       % discard prev
diary(fnam);                   % start logging screen text
omg.dt=0.1; memorygraph('start',omg);
for i=1:numel(Ns), N=Ns(i);
  nam = sprintf('N=%d, M=%d, tol=%g',N,ceil(N*dataratio),opts.tol);
  fprintf(['\n' nam '......................\n'])
  memorygraph('label',nam);
  outs{i} = test_iNUDFT(signal,N,dataratio,solvlist,nudist,opts);  % dummy out for now
  diary off                      % appends to file
  [bytes est_times cpu_times cpu_usages labelstrings labeltimes] = memorygraph('plot'); drawnow;
end
memorygraph('done');
