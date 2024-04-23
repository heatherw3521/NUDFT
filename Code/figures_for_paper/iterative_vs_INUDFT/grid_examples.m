% depends on external pkgs: chebfun, export_fig
% Wilber.
% Barnett fixed n=m/2 issue 4/3/24.

%grid styles: 
%%
clear all
close all
%%
m = 150; 
nds1 = linspace(0, 1, m).' + .5/m*(2*rand(m,1)-1); 
clf
plot(nds1, ones(m,1),'.k', 'markersize', 25)
hold on
plot([0, 1], [1, 1], '--k')
axis off
%%
nds2 = chebpts(m, [0, 1]); 
%nds = nds(1:end-1); 
clf
plot(nds2, ones(m,1),'.k', 'markersize', 25)
hold on
plot([0, 1], [1, 1], '--k')
axis off
%%
nds3 = rand(m,1); 
clf
plot(nds3, ones(m,1),'.k', 'markersize', 25)
hold on
plot([0, 1], [1, 1], '--k')
axis off
%%
sbwp = 8.0;   % gap space-bandwidth prod in half-wavelengths
n = m/2;         % alex added
nds4 = rand(m,1)*(1-sbwp/n); %assume num modes = m/2 ...but was paren typo sbwp/m/2 not sbwp/(m/2) :(
clf
plot(nds4, ones(m,1),'.k', 'markersize', 25)
hold on
plot([0, 1], [1, 1], '--k')
axis off

%%
% ALL TOGETHER: 
clf
plot(nds4, 0*ones(m,1),'.k', 'markersize', 25)
hold on
plot(nds3, .2*ones(m,1),'.k', 'markersize', 25)
plot(nds2, .4*ones(m,1),'.k', 'markersize', 25)
plot(nds1, .6*ones(m,1),'.k', 'markersize', 25)
plot([0, 1], 0*[1, 1], '--k')
plot([0, 1], .2*[1, 1], '--k')
plot([0, 1], .4*[1, 1], '--k')
plot([0, 1], .6*[1, 1], '--k')
axis equal
axis off
set(gcf,'color','w')
%export_fig example_samples.pdf
%%
% use circles: to deal with ugly spacing, we split into
% four images: 
sz = 20; 
clf
circ = exp(1i*2*pi*linspace(0, 1, 300));
gry = [.4 .4 .4]; 
plot(circ, '-k')
hold on
plot(exp(1i*2*pi*nds1),'.', 'color', gry, 'markersize', sz, 'markerfacecolor', gry)
%title('1. jitttered','Position', [0, -1.4, 0])
axis equal
axis off
set(gcf,'color','w')
%export_fig example_samples_circ1.pdf

%%
clf
plot(circ, '-k')
hold on
plot(exp(1i*2*pi*nds2),'.', 'color', gry, 'markersize', sz, 'markerfacecolor', gry)
%title('2. Clenshaw-Curtis','Position', [0, -1.4, 0])
axis equal
axis off
set(gcf,'color','w')
%export_fig example_samples_circ2.pdf
%%

clf
plot(circ, '-k')
hold on
plot(exp(1i*2*pi*nds3),'.', 'color', gry, 'markersize', sz, 'markerfacecolor', gry)
%title('3. Random','Position', [0, -1.4, 0])
axis equal
axis off
set(gcf,'color','w')
%export_fig example_samples_circ3.pdf

%%
clf
plot(circ, '-k')
hold on
plot(exp(1i*2*pi*nds4),'.', 'color', gry, 'markersize', sz, 'markerfacecolor', gry)
%title('4. Random with gaps','Position', [0, -1.4, 0])
axis equal
axis off
set(gcf,'color','w')
export_fig example_samples_circ4.pdf

%export_fig example_samples_circ.pdf
