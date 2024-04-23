%test iterative methods: 
% use AB's test suite to create an image comparing iterative methods
% vs conditioning of matrix

clear all; 
close all; 
n = 2^13; 

jitterch = [.25, .5, 1, 10,50,100];
signal = 'randn';
dataratio = 2; 
M = dataratio*n; 
w = []; condflag = 1; % compute condition numbers (takes a long time)
%%
%for j = 1:(length(jitterch)+2)

for j = 1:(length(jitterch))
%for j = 1:1 
    if j <= length(jitterch)
        jitter = jitterch(j); 
        x = (2*pi/M)*((0:M-1) + jitter*(2*rand(1,M)-1));
    elseif j == (length(jitterch)+1)
        x = 2*pi*rand(1,M);
    else
        sbwp = 8.0;   % gap space-bandwidth prod in half-wavelengths
        x = 2*pi*rand(1,M)*(1-sbwp/n);
    end
    for k = 1:2
    [relresid, relerror, times, numits, condA] = iterative_test(signal,n,dataratio,x, condflag);
    end
    relres(:,j) = relresid(:); 
    relerr(:,j) = relerror(:); 
    timing(:,j) = times(:); 
    numiter(:,j) = numits(:); 
    condition(j) = condA; 
    j
end

save('iter_fig_cond', 'relres', 'relerr', 'timing', 'numiter', 'condition')
%%
% now create a figure with timing
% leave out solver 4. 
[~,ord] = sort(condition);
dirtime = mean(timing(end,:));
loglog(condition(ord), timing([1:2, 4:end-1], ord).', '.-', 'linewidth', 2.5, 'markersize', 20);
hold on
loglog([condition(1), condition(end)], [dirtime, dirtime], '--k', 'linewidth', 2.5)
legend('CG nor', 'Strang PCG nor', 'FP adj wei', 'CG adj nor', 'sinc2 PCG adj nor', 'INUDFT','location', 'northwest')

