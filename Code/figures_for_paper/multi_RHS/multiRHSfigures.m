load("dir_multiRHS.mat")
%%
% we want to match colors used in the dir solver plots: 

colors = get(gca, 'ColorOrder'); 


clf
loglog(sz, timings(1,:), '.-', 'markersize', 25, 'linewidth', 2, 'color', colors(1,:))
hold on
loglog(sz, timings(2,:), '.-', 'markersize', 25, 'linewidth', 2, 'color', colors(4,:))
loglog(sz, timings(3,:), '.--', 'markersize', 25, 'linewidth', 2, 'color', 'k')
axis tight
set(gcf,'color','w')
set(gca, 'fontsize', 14)
%%