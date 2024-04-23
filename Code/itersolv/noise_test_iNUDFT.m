% Script to run multiple noise-levels test_iNUDFT
% Barnett 9/14/23

clear
signal = 'image';
N=2^11;
dataratio = 1.8;    % as for run 11/18/22
solvlist = [1:2, 6, 9:10];   % only include 1 if N < few k.
nudists = [1 4];
noises = [0 1e-6 1e-4 1e-2];  % small set for now
fnam = sprintf('results/%s_N%d_dataratio%g_noises.txt',signal,N,dataratio)
system(['rm -f ',fnam]);       % discard prev
diary(fnam);                   % start logging screen text
opts = [];
for i=1:numel(noises), noise=noises(i);
  fprintf('\nN=%d, dataratio(M/N)=%g, noise=%g......................................\n',N,dataratio,noise);
  opts.addnoise = noise;
  outs{i} = test_iNUDFT(signal,N,dataratio,solvlist,nudists,opts);
  diary off                      % appends to file
end
