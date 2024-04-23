% plots for iterative methods with different node sets

%solvers included are 
% row 1 = cgn normal
% row 2 = cgn normal with strang precond
% row 4 = FP adj weights sinc quad (GG, Inati)
% row 5 = CG on adj eqns
% row 6 = inudft
% (check tols in the itersolve/iterative_test.m)


% vars are 
%numiter = # of iterations (when relevant)
%relres = relative residual
%relerr = relative error
%timing = wall clock time for solve

% 1st col = test with jittered nodes
% 2nd col = test with CC quadrature nodes
% 3rd col = test with rand unif. iid nodes
% 4th col = rand unif. iid nodes with gaps. 

%N = 2^18
% M = 2*N
%%
%%
clear all; close all
load('iter_4nodesets_decay_BIG.mat') 
dirtime6 = mean(timing(end,:));
dirres6 = mean(relres(end,:));

%%
% plot residual vs time

%basic vis:
loglog(timing.', relres.','o-', 'linewidth', 1, 'markersize', 5)
legend('cgnormal', 'cgnor strang','FP', 'CGad','CGad_sinc', 'inudft' )

%% 
% plot with separate markers
C = linspecer(10); %colors
clf
loglog(timing(1,:),relres(1,:), 'o-', 'color', 'k', 'markersize', 6, 'markerfacecolor', 'k');
hold on
loglog(timing(2,:),relres(2,:), '*-', 'color', C(2,:), 'markersize', 8);
loglog(timing(3,:), relres(3,:),'+-', 'color', C(8,:), 'markersize', 8);
loglog(timing(4,:), relres(4,:),'s-', 'color', C(1,:), 'markersize', 8, 'markerfacecolor', C(1,:));
loglog(timing(5,:), relres(5,:),'d-', 'color', C(9,:), 'markersize', 8, 'markerfacecolor', C(9,:));
loglog(dirtime6, dirres6,'p', 'color', 'r', 'markersize', 10, 'markerfacecolor','r');
axis tight
set(gca, 'fontsize', 14)
legend('cgnor', 'pcgnor (strang)','FP+sinc', 'cgadnor','pcgadnor (sinc)', 'inudft', 'location', 'northwest')
set(gcf,'color','w')
%export_fig iterative_1RHS_relresvstimings.pdf
%%
% BAR PLOTS timing: 
clf
bar(timing(1:end-1,:).')
set(gca,'Yscale','log')
hold on
semilogy([0, 5], [dirtime6, dirtime6], '--k', 'linewidth', 2)
%hold on
%semilogy([0, 5], [dirtime11, dirtime11], '--k', 'linewidth', 2)
set(gcf,'color','w')
set(gca, 'fontsize', 14)
xlim([.5, 4.5])
legend('CG nor', 'PCG nor strang', 'FP adj-sinc', 'CG adj', 'PCG adj sinc', 'location', 'northwest')
%legend off
export_fig iterative_solvers_4_nodes_decay218_timingsBAR.pdf
%%
% BAR PLOTS resid: 
clf
bar(relres(1:end-1,:).')
set(gca,'Yscale','log')
hold on
semilogy([0, 5], [dirres6, dirres6], '--k', 'linewidth', 2)
%hold on
%semilogy([0, 5], [dirtime11, dirtime11], '--k', 'linewidth', 2)
set(gcf,'color','w')
set(gca, 'fontsize', 14)
xlim([.5, 4.5])
legend('CG nor', 'PCG nor strang', 'FP adj-sinc', 'CG adj', 'PCG adj sinc', 'location', 'northwest')
%legend off
export_fig iterative_solvers_4_nodes_decay218_residsBAR.pdf
