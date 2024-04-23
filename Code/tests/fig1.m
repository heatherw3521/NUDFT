ns = 2 .^ (9:20);
num_rhs = 10;
opts.tol = 1e-12;
opts.cgtol = 1e-6;
opts.hsstol = 1e-6;
opts.addnoise = 1e-8;
methods = {@solv_CGN,@wrapper_INUDFT,@solv_FDToepN};
methodnames = {'CGN', 'Direct', 'Toeplitz'};
methodmaxs = [50000 Inf Inf];
times = zeros(length(ns), length(methods));
mems = zeros(length(ns), length(methods));
addpath('../itersolv/')

for method_idx = 1:length(methods)
    method = methods{method_idx};
    fprintf('%s\n',methodnames{method_idx});
    for n_idx = 1:length(ns)
        n = ns(n_idx);
        if n > methodmaxs(method_idx)
            break
        end
        fprintf('\t%d\n',n)
        m = 2*n;
        
        % Generate problem
        % [nodes, ~] = clencurt(m-1);
        nodes = sort(2*pi*linspace(0,1-1/m,m) + 0.01*randn(1,m));
        freqinds = (-n/2:n/2-1)';
        C0 = exp(-10.0*abs(freqinds/n)) .* (randn(n,num_rhs)+1i*randn(n,num_rhs));
        B = finufft1d2(nodes,+1,opts.tol,C0);
        B = B + opts.addnoise*randn(m,num_rhs);

        tic;
        [x,data] = method(nodes, n, B, opts);
        times(n_idx, method_idx) = toc;
        mems(n_idx, method_idx) = getMemSize(data);
    end
end

%% Plot
close all
figure(1)
loglog(ns, times(:,1)); hold on
loglog(ns, times(:,2),'--')
loglog(ns, times(:,3),':')
xlabel('Number of columns $n$')
ylabel('Time (sec)')
legend('CGNE', 'Direct (Ours)', 'ToeplitzNE')

figure(2)
loglog(ns, mems(:,1)/1e9); hold on
loglog(ns, mems(:,2)/1e9,'--');
loglog(ns, mems(:,3)/1e9,':');
xlabel('Number of columns $n$')
ylabel('Memory (GB)')