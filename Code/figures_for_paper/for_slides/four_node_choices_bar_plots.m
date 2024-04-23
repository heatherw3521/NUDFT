% bar graph figures for iterative tests
% see Code/iterative_four_node_choices.m

%solvers included are 
% row 1 = cgn normal
% row 2 = cgn normal with strang precond
% row 3 = adj matvec sinc quad weights
% row 4 = FP adj weights sinc quad (GG, Inati)
% row 5 = CG on adj eqns
% row 6 = sinc2 PCG on adj eqns
% row 7 = inudft with hsstol = cgtol (1e-6)
% row 8 = inudft with hsstol = 1e-11

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
% leave out solver 3, which is using wei. adj AWA^*, not 
% an iterative method. 
%%
load('iter_4nodesets_decay_BIG.mat')
numiter(3,:) = []; 
relres(3,:) =[]; 
relerr(3,:) = []; 
timing(3,:) = []; 
dirtime6 = mean(timing(end-1,:));
dirtime11 = mean(timing(end, :));

[timingsorted,srtidx] = sort(timing, 2);
relressorted = zeros(size(relres));
for j = 1:7
  relressorted(j,:) = relres(j,srtidx(j));
end

%%
clf
loglog(timingsorted(1:end-2,:).' ,relres(1:end-2,:).', '.-', 'linewidth', 2, 'markersize', 30)
legend('cg nor', 'pcg nor strang', 'FP adj. sinc wei', 'cg adj', 'sinc2 pcg adj', 'location', 'northwest')



%%
bar(timing(1:end-2,:).')
set(gca,'Yscale','log')
hold on
%semilogy([0, 5], [dirtime6, dirtime6], '--k', 'linewidth', 2)
%hold on
semilogy([0, 5], [dirtime11, dirtime11], '--k', 'linewidth', 2)
set(gcf,'color','w')
set(gca, 'fontsize', 14)
xlim([.5, 4.5])
legend('cg nor', 'pcg nor strang', 'FP adj. sinc wei', 'cg adj', 'sinc2 pcg adj')
legend off
%export_fig iterative_solvers_4_nodes_decay218.pdf
%%
clf
bar(relres(1:end-2,:).')
set(gca,'Yscale','log')
hold on
semilogy([1, 2, 3, 4], relres(end,:), '.--k', 'linewidth', 2, 'markersize', 30)
%semilogy([1, 2, 3, 4], relres(end,:), '.--b', 'linewidth', 2, 'markersize', 30)
set(gcf,'color','w')
set(gca, 'fontsize', 14)
xlim([.5, 4.5])
%legend('cg nor', 'pcg nor strang', 'FP adj. sinc wei', 'cg adj', 'sinc2 pcg adj')
legend off


