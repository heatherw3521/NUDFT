% Script to run multiple test_iNUDFT runs and save outputs as text files for now
% Barnett 2/2/23

clear
signal = 'decay';
Ns = 2.^[11,14,18];
dataratio = 1.8;    % as for run 11/18/22
for i=1:numel(Ns), N=Ns(i);
  solvlist = [2:7,9:10];
  if N<=2^11, solvlist = [1,solvlist]; end  % small case only check dense cond #s
  fnam = sprintf('results/%s_N%d_dataratio%g.txt',signal,N,dataratio)
  system(['rm -f ',fnam]);       % discard prev
  diary(fnam);                   % start logging screen text
  fprintf('\nN=%d, dataratio(M/N)=%g......................................\n',N,dataratio);
  outs{i} = test_iNUDFT(signal,N,dataratio,solvlist);
  diary off                      % appends to file
end
