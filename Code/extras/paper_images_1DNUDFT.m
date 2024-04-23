%some figures for papers 

clear all
close all
n = 2^4;
%%
%roots of unity
w = exp(pi*1i/n); 
r = w.^(2*(1:n)); 
rang = angle(r);
%r = r(s); 
%%

%cluster boundaries:
y1 = w.^(2*(1:n) + 1) ; y1 = y1(:);


figure(1)
clf
set(gca, 'fontsize', 18)
set(gcf,'color','w')
%title('node locations')
hold on
%%
%start with a polar grid:
rr = linspace(0, 1, n/2); 
t = trigpts(101, [0, 2*pi]); t = [t; 0];

[TT, RR] = meshgrid(t,rr);
[TT1, RR1] = meshgrid(rang,rr);
%[TT1, RR1] = meshgrid([trigpts(n,dom(1:2)); dom(2)],r);
XX = RR.*cos(TT);
YY = RR.*sin(TT);
XX1 = RR1.*cos(TT1);
YY1 = RR1.*sin(TT1);
grdgry = [120 120 120]/256;
clrGrid = [0.542968750000000, 0, 0];
plot(XX1,YY1,'-','Color',grdgry);
%plot(XX1,YY1,'-');
hold on
plot(XX',YY','-','Color',grdgry);
hold on
plot(exp(1i*linspace(0, 2*pi, 100)), 'k-')
axis square
axis off
%%
%now plot the boundaries ticks: 
rdist = [0.95, 1.05]; 
ptset1 = y1.*rdist;  
ptset2 = circshift(ptset1,1); 
y2 = circshift(y1,1);


%%
for j = 1:n
    plot(ptset1(j,:), '-', 'color', 'k', 'linewidth', 2.2 )
    %plot(ptset2(j,:), '--', 'color', grdgry )
    %tt = rdist(2)*exp(1i*linspace(mod(angle(y1(j)), 2*pi), mod(angle(y2(j)), 2*pi), 100));
    %plot(tt, '--', 'color', 'k' )
end
%%
% add in last piece: 
   %plot(ptset1(1,:), '--', 'color', grdgry )
   %plot(ptset2(1,:), '--', 'color', grdgry )
   %tt = rdist(2)*exp(1i*linspace(angle(ptset2(1,2)), angle(ptset1(1,2)), 100));
   %plot(tt, '--', 'color', grdgry )
   
%%
% add real and imag axes: 
axis equal
hold on
plot([0, 0], [-1.1, 1.1], '-', 'color', 'k', 'linewidth', 1)
plot([-1.1, 1.1], [0 0], '-', 'color', 'k', 'linewidth', 1)
%%
%now plot some clusters: 
clr1 = [53, 100, 176]./255;
clr2 = [120, 24, 89]./255;

beta = 4/7;
nclust = []; 
hold on
for j = 1:n
    k = max(8,randi(20)); %number of nodes in cluster
    cl = 2*j + (2*rand(k,1)-1)*beta; %random locations in cluster
    cl = w.^(cl); 
    if mod(j,2)
        clr = clr1;
    else
        clr = clr2;
    end
    plot(cl, '.', 'markersize', 12, 'color', clr);
    shg
    nclust(j) = k; 
end
%%
% now plot roots
%clrt = [194, 124, 4]./255; 
%plot(r, 'sq', 'markersize', 6, 'color', clrt, ...
%   'MarkerFaceColor',clrt); 


ylim([-1.1, 1.1])
xlim([-1.1, 1.1])
axis equal
hold off
%%
export_fig -m2 -pdf 'nodes_ROU'
%%
%let's add a matrix picture: 
m = sum(nclust); 

xborders = [ 1 : n, n * ones(1, m), n:-1:1, zeros(1, m) ];
yborders = [ zeros(1,n), 1:m, m * ones(1, n), m:-1:1 ];
%%
st = 0; 
clf
for j = 1:n
    gap = nclust(j);
    x = [1, 5*n, 5*n, 1]; 
    y = [st, st, st-gap+1, st-gap+1]; 
    st = st-gap+1; 
    if mod(j,2)
        clr = clr1;
    else
        clr = clr2;
    end
    fill(x,y,clr)
    hold on
end
axis off
axis equal
%%
export_fig -m2 -pdf 'rows_clustered'
%%


%%
%x = exp(-1i*rand(n,1)*2*pi); %random perturbation
%x = exp(-1i*rand(floor(60*n/100),1)*pi/4); %random perturbation with clustering
%x = [x ;exp(-1i*rand(n-length(x),1)*2*pi)];
%x = exp(2*pi*1i*((0:n-1)/n + 1/(2*n))); %perturbed by 1/2 (sanity check)


% %%
% %bounds on singular values of submatrices: 
% clear all
% n = 2^8;
% %x = exp(-1i*rand(n,1)*2*pi); %random perturbation
% %perturbation bounded by p 
% p = .25; 
% signs = (2*(rand(n,1)>.5) - 1);
% x = exp(-2*pi*1i*((0:n-1).'/n + p/n*rand(n,1).*signs));
% 
% %%
% %build matrix
% %[~, s] = sort(angle(x)); 
% %x = x(s); 
% V = bsxfun(@(x,y) x.^y, x, 0:n-1); 
% V = (fft(V'))'; %apply fft on the right
% 
% %submatrix 1: HODLR type 
% I1 = 1:n/4;
% I2 = (n/4+1):n/2;
% S1 = V(I1, I2); 
% s1 = svd(S1); s1 = s1/s1(1); 
% %%
% %bound: 
% p = exp(pi^2/4/log(8*n)); 
% 
% k = [2:40];
% nn= length(k); 
% 
% %plot the bound
% figure(2)
% set(gca, 'fontsize', 18)
% set(gcf,'color','w')
% semilogy(k, s1(2:nn+1), '.k','markersize', 20)
% hold on
% semilogy(k, 4*p.^(-2*(k-1)), '-k', 'linewidth', 2)





 