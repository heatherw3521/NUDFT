n = 2^16; % 2 ^ 17;
m = 2*n;

num_rhs = 10;
opts.tol = 1e-15;
opts.cgtol = 1e-6;
opts.hsstol = 1e-6;
opts.addnoise = 0;
        
% Generate problem
% [nodes, ~] = clencurt(m-1);
nodes = sort(2*pi*linspace(0,1-1/m,m) + 0.01*randn(1,m));
freqinds = (-n/2:n/2-1)';
C0 = exp(-10.0*abs(freqinds/n)) .* (randn(n,num_rhs)+1i*randn(n,num_rhs));
B = finufft1d2(nodes,+1,opts.tol,C0);
B = B + opts.addnoise*randn(m,num_rhs);

methods = {@wrapper_INUDFT,@solv_FDToepN};
methodnames = {'Direct', 'Toeplitz'};
tols = logspace(-14, -2, 13);
times = zeros(length(tols), length(methods));
mems = zeros(length(tols), length(methods));
accs = zeros(length(tols), length(methods));
addpath('../itersolv/')

for method_idx = 1:length(methods)
    method = methods{method_idx};
    fprintf('%s\n',methodnames{method_idx});
    for tol_idx = 1:length(tols)
        opts.hsstol = tols(tol_idx);
        fprintf('\t%e\n',opts.hsstol);
        tic;
        [X,data] = method(nodes, n, B, opts);
        times(tol_idx, method_idx) = toc;
        mems(tol_idx, method_idx) = getMemSize(data);
        accs(tol_idx, method_idx) = norm(X - C0,'fro') / norm(C0, 'fro');
    end
end

%% Plot
close all
figure(1)
loglog(tols, times(:,1), '--', 'Color', "#D95319"); hold on
loglog(tols, times(:,2),':', 'Color', "#EDB120")
xlabel('HSS compression tolerance')
ylabel('Time (sec)')
legend('Direct (Ours)', 'ToeplitzNE')
saveas(gcf, '../figs/fig2_times.png')
saveas(gcf, '../figs/fig2_times.fig')

figure(2)
loglog(tols, mems(:,1)/1e9, '--', 'Color', "#D95319"); hold on
loglog(tols, mems(:,2)/1e9,':','Color',"#EDB120");
xlabel('HSS compression tolerance')
ylabel('Memory (GB)')
saveas(gcf, '../figs/fig2_mems.png')
saveas(gcf, '../figs/fig2_mems.fig')

figure(3)
loglog(tols, accs(:,1), '--', 'Color', "#D95319"); hold on
loglog(tols, accs(:,2),':','Color',"#EDB120");
xlabel('HSS compression tolerance')
ylabel('Relative error')
saveas(gcf, '../figs/fig2_accs.png')
saveas(gcf, '../figs/fig2_accs.fig')