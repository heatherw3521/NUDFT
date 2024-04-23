%image figs and noise

clear all
close all
signal = 'image';
%signal = 'decay';
N=2^11;
dataratio =1.8;    % as for run 11/18/22
solvlist = [1:2, 6, 9:10];   % only include 1 if N < few k.
nudists = [1 4];
noises = [0 1e-6 1e-4 1e-2];  % small set for now
%fnam = sprintf('results/%s_N%d_dataratio%g_noisesT.txt',signal,N,dataratio)
%system(['rm -f ',fnam]);       % discard prev
%diary(fnam);                   % start logging screen text
opts = [];
opts.verb = 1; 
opts.outputstruct = 1; 
opts.hsstol = 1e-10; 
for i=1:numel(noises), noise=noises(i);
  fprintf('\nN=%d, dataratio(M/N)=%g, noise=%g......................................\n',N,dataratio,noise);
  opts.addnoise = noise;
  outputs{i} = test_iNUDFT(signal,N,dataratio,solvlist,nudists,opts);
  %diary off                      % appends to file
end

%save('decay_test_noise', 'outputs', 'noises', 'dataratio', 'N')
%save('image_test_noise', 'outputs', 'noises', 'dataratio', 'N')