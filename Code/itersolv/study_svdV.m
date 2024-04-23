% Study sing vals of type-II system mat V, to get at CG nor conv type.
% Barnett 4/18/24

clear
signal = 'randn';
N = 2^11;
dataratio = 1.8;

solvlist = [1 2];
opts.verb = 0;     % 1 = make plots
opts.cgtol=1e-12;
opts.outputstruct=1;
%opts.addnoise = 1e-2;
nudistlist = [1 2 3 4];
out = test_iNUDFT(signal,N,dataratio,solvlist,nudistlist,opts);
for j = 1:numel(nudistlist)
  nudist = nudistlist(j)
  x = out{j,1}.x;      % get NU pts on [0,2pi)
  A = densemat_nudft(x,N);   % rebuild A (=V)
  [U,S,V] = svd(A,0);
  sv = diag(S);
  figure(nudist);
  subplot(4,1,1); semilogy(1:N,sv,'+'); title(sprintf("nudist %d: sigma_j(V)",nudist))
  subplot(4,1,2); plot(sv, 0*sv, '+'); title("sing vals V")
  v = axis; axis([0 v(2:4)]);    % enclose origin
  subplot(4,1,3); hist(sv,20); title("histogram sing vals V")
  v = axis; axis([0 v(2:4)]);
  subplot(4,1,4); hist(log10(sv),20); title("histogram log10(sing vals V)")
  drawnow
end
% saved group of 4 figs as results/singval_spec_*.png
% results: nudist=3,4 roughly flat sing val dist over [0,100]

