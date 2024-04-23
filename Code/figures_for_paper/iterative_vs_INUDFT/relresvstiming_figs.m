% The relres vs timing plot: 

% last run 3/5/24: n = 2^14, m  = ceil(1.8*n), maxits = 10000
%%
%solvers included are 
% row 1 = cgn normal
% row 2 = cgn normal with strang precond
% row 3 = FP adj weights sinc quad (GG, Inati)
% row 4 = CG on adj eqns
% row 5 = CG on adj with sinc2 precond
% row 6 = inudft
% row 7 = inudft with 10 rhs (timings are per RHS)
% row 8 = inudft with 100 rhs
% row 9 = inudft with 1000 rhs
%%
load('iter_rand_decay_214.mat')
%%
% plot with separate markers
C = linspecer(6); %colors
CO = get(gca, 'ColorOrder');
relres = flip(relres, 2); %we want bigger residual to the right
timing = flip(timing, 2);
%%
clf
loglog(relres(1,:), timing(1,:),'o-', 'color', CO(4,:), 'markersize', 6, 'markerfacecolor', CO(4,:));
hold on
%loglog(relres(2,:), timing(2,:), '*-', 'color', CO(2,:), 'markersize', 8);
loglog( relres(3,:),timing(3,:),'+-', 'color', CO(6,:), 'markersize', 8);
loglog( relres(4,:), timing(4,:),'s-', 'color', CO(1,:), 'markersize', 8, 'markerfacecolor', CO(1,:));
loglog( relres(5,:), timing(5,:),'d-', 'color', CO(3,:), 'markersize', 8, 'markerfacecolor', CO(3,:));
loglog( relres(6,:), timing(6,:),'p-', 'color', 'k', 'markersize', 8, 'markerfacecolor', 'k');
loglog( relres(7,:), timing(7,:),'^-', 'color', 'k', 'markersize', 8, 'markerfacecolor', 'k');
loglog( relres(8,:), timing(8,:),'<-', 'color', 'k', 'markersize', 8, 'markerfacecolor', 'k');
loglog( relres(9,:), timing(9,:),'>-', 'color', 'k', 'markersize', 8, 'markerfacecolor', 'k');
%axis tight
set(gca, 'fontsize', 14)
legend('CG nor','FP adj-sinc', 'CG adj','PCG adj sinc', '', '', '', '', 'location', 'northwest')
set(gcf,'color','w')
%export_fig iterative_1RHS_relresvstimings.pdf
xlim([1e-10, .02439])
grid on
%export_fig relresvstimingsplots.pdf